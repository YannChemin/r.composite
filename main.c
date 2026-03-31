/*****************************************************************************
 *
 * MODULE:       r.composite
 * AUTHOR(S):    Glynn Clements - glynn.clements@virgin.net
 *               Modern GIS quality improvements
 * PURPOSE:      Combines red, green, and blue raster maps into a
 *               true-colour composite (default: three FCELL maps) or a
 *               palette-indexed CELL map (-p).
 *
 *               Improvements over original:
 *                 true-colour FCELL output (default, -p for legacy palette)
 *                 perceptual Oklab quantisation (-k, palette mode)
 *                 gamma-correct quantisation (-g, palette mode)
 *                 Bayer ordered dithering (-b, palette mode)
 *                 per-band NULL tracking and fill (null_value=)
 *                 adaptive palette levels (levels=auto)
 *                 OpenMP column-level parallelism (non-dithering paths)
 *
 * COPYRIGHT:    (C) 2001-2024 by the GRASS Development Team
 *
 *               This program is free software under the GNU General Public
 *               License (>=v2). Read the file COPYING that comes with GRASS
 *               for details.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>

/* ═══════════════════════════════════════════════════════════════════════════
 * Oklab colour-space matrices
 * Reference: https://bottosson.github.io/posts/oklab/
 * ═══════════════════════════════════════════════════════════════════════════ */

/* linear sRGB → LMS */
static const float M1[3][3] = {
    {0.4122214708f, 0.5363325363f, 0.0514459929f},
    {0.2119034982f, 0.6806995451f, 0.1073969566f},
    {0.0883024619f, 0.2817188376f, 0.6299787005f},
};
/* LMS^(1/3) → Lab */
static const float M2[3][3] = {
    { 0.2104542553f,  0.7936177850f, -0.0040720468f},
    { 1.9779984951f, -2.4285922050f,  0.4505937099f},
    { 0.0259040371f,  0.7827717662f, -0.8086757660f},
};
/* Lab → LMS^(1/3) */
static const float M2_inv[3][3] = {
    {1.0f,  0.3963377774f,  0.2158037573f},
    {1.0f, -0.1055613458f, -0.0638541728f},
    {1.0f, -0.0894841775f, -1.2914855480f},
};
/* LMS → linear sRGB */
static const float M1_inv[3][3] = {
    { 4.0767416621f, -3.3077115913f,  0.2309699292f},
    {-1.2684380046f,  2.6097574011f, -0.3413193965f},
    {-0.0041960863f, -0.7034186147f,  1.7076147010f},
};

/* Safe quantisation domain for Oklab a and b axes (covers the sRGB gamut
 * with margin; actual sRGB maps into roughly a∈[-0.23,+0.28], b∈[-0.31,+0.20]) */
#define OKLAB_A_MIN (-0.5f)
#define OKLAB_A_MAX  (0.5f)
#define OKLAB_B_MIN (-0.5f)
#define OKLAB_B_MAX  (0.5f)

/* ═══════════════════════════════════════════════════════════════════════════
 * 8×8 Bayer ordered-dither threshold matrix (values 0–63)
 * ═══════════════════════════════════════════════════════════════════════════ */
static const int bayer8[8][8] = {
    { 0, 32,  8, 40,  2, 34, 10, 42},
    {48, 16, 56, 24, 50, 18, 58, 26},
    {12, 44,  4, 36, 14, 46,  6, 38},
    {60, 28, 52, 20, 62, 30, 54, 22},
    { 3, 35, 11, 43,  1, 33,  9, 41},
    {51, 19, 59, 27, 49, 17, 57, 25},
    {15, 47,  7, 39, 13, 45,  5, 37},
    {63, 31, 55, 23, 61, 29, 53, 21},
};

/* ═══════════════════════════════════════════════════════════════════════════
 * Data structures
 * ═══════════════════════════════════════════════════════════════════════════ */

struct band {
    struct Option *opt_name;
    struct Option *opt_levels;
    char *name;
    int levels;
    int maxlev;
    int offset;
    int file;
    int type;
    int size;
    unsigned char *array[3]; /* channel arrays from colour-table lookup */
    short *floyd[2];         /* Floyd-Steinberg error rows */
    struct Colors colors;
};

static const char *const color_names[3] = {"red", "green", "blue"};

static struct band B[3];
static int closest;
static int use_gamma;
static int use_oklab;

static unsigned char *dummy; /* shared sink for unused colour channels */

/* ═══════════════════════════════════════════════════════════════════════════
 * Utility
 * ═══════════════════════════════════════════════════════════════════════════ */

static inline int clamp_byte(int v)
{
    return v < 0 ? 0 : v > 255 ? 255 : v;
}

static inline float clamp_f(float v)
{
    return v < 0.0f ? 0.0f : v > 1.0f ? 1.0f : v;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * sRGB ↔ linear-light (IEC 61966-2-1)
 * ═══════════════════════════════════════════════════════════════════════════ */

static inline float srgb_to_linear(int v)
{
    float f = v / 255.0f;
    return f <= 0.04045f ? f / 12.92f
                         : powf((f + 0.055f) / 1.055f, 2.4f);
}

static inline int linear_to_srgb(float f)
{
    f = clamp_f(f);
    float s = f <= 0.0031308f ? 12.92f * f
                              : 1.055f * powf(f, 1.0f / 2.4f) - 0.055f;
    return (int)(s * 255.0f + 0.5f);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Oklab conversions
 * ═══════════════════════════════════════════════════════════════════════════ */

static void rgb_to_oklab(float r, float g, float b,
                          float *oL, float *oa, float *ob)
{
    float l = M1[0][0]*r + M1[0][1]*g + M1[0][2]*b;
    float m = M1[1][0]*r + M1[1][1]*g + M1[1][2]*b;
    float s = M1[2][0]*r + M1[2][1]*g + M1[2][2]*b;

    l = cbrtf(l > 0.0f ? l : 0.0f);
    m = cbrtf(m > 0.0f ? m : 0.0f);
    s = cbrtf(s > 0.0f ? s : 0.0f);

    *oL = M2[0][0]*l + M2[0][1]*m + M2[0][2]*s;
    *oa = M2[1][0]*l + M2[1][1]*m + M2[1][2]*s;
    *ob = M2[2][0]*l + M2[2][1]*m + M2[2][2]*s;
}

static void oklab_to_rgb(float L, float a, float b,
                          float *or_, float *og, float *ob_out)
{
    float l = M2_inv[0][0]*L + M2_inv[0][1]*a + M2_inv[0][2]*b;
    float m = M2_inv[1][0]*L + M2_inv[1][1]*a + M2_inv[1][2]*b;
    float s = M2_inv[2][0]*L + M2_inv[2][1]*a + M2_inv[2][2]*b;

    l = l*l*l;
    m = m*m*m;
    s = s*s*s;

    *or_    = M1_inv[0][0]*l + M1_inv[0][1]*m + M1_inv[0][2]*s;
    *og     = M1_inv[1][0]*l + M1_inv[1][1]*m + M1_inv[1][2]*s;
    *ob_out = M1_inv[2][0]*l + M1_inv[2][1]*m + M1_inv[2][2]*s;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Quantisation helpers
 * ═══════════════════════════════════════════════════════════════════════════ */

/* Standard per-channel sRGB quantisation (original algorithm). */
static int quantize(int c, int x)
{
    return closest ? (x + B[c].offset) * B[c].maxlev / 256
                   : x * B[c].levels / 256;
}

/* Gamma-correct quantisation: linearise, quantise uniformly in linear light. */
static int quantize_gamma(int c, int x)
{
    float lin = srgb_to_linear(x);
    int q = (int)(lin * B[c].levels);
    if (q >= B[c].levels)
        q = B[c].maxlev;
    return q;
}

/* Whole-pixel Oklab quantisation.
 * Returns a single packed index encoding L, a, b quantisation levels.
 * Index layout: (b_idx * na + a_idx) * nL + L_idx, where nL = B[0].levels,
 * na = B[1].levels, nb = B[2].levels. */
static int quantize_oklab(int xr, int xg, int xb)
{
    int nL = B[0].levels;
    int na = B[1].levels;
    int nb = B[2].levels;

    float lr = srgb_to_linear(xr);
    float lg = srgb_to_linear(xg);
    float lb = srgb_to_linear(xb);

    float L, a, bv;
    rgb_to_oklab(lr, lg, lb, &L, &a, &bv);

    int l_idx = (int)(L * nL);
    if (l_idx < 0)    l_idx = 0;
    if (l_idx >= nL)  l_idx = nL - 1;

    int a_idx = (int)((a - OKLAB_A_MIN) / (OKLAB_A_MAX - OKLAB_A_MIN) * na);
    if (a_idx < 0)    a_idx = 0;
    if (a_idx >= na)  a_idx = na - 1;

    int b_idx = (int)((bv - OKLAB_B_MIN) / (OKLAB_B_MAX - OKLAB_B_MIN) * nb);
    if (b_idx < 0)    b_idx = 0;
    if (b_idx >= nb)  b_idx = nb - 1;

    return (b_idx * na + a_idx) * nL + l_idx;
}

/* Recover sRGB triple for a given Oklab palette index.
 * Used by Floyd-Steinberg to compute per-channel quantisation error. */
static void oklab_idx_to_srgb(int idx, int *r, int *g, int *b)
{
    int nL = B[0].levels;
    int na = B[1].levels;

    int l_idx = idx % nL;
    int a_idx = (idx / nL) % na;
    int b_idx = idx / (nL * na);

    float L  = (B[0].maxlev > 0) ? (float)l_idx / B[0].maxlev : 0.5f;
    float av = (B[1].maxlev > 0)
               ? (float)a_idx / B[1].maxlev * (OKLAB_A_MAX - OKLAB_A_MIN) + OKLAB_A_MIN
               : 0.0f;
    float bv = (B[2].maxlev > 0)
               ? (float)b_idx / B[2].maxlev * (OKLAB_B_MAX - OKLAB_B_MIN) + OKLAB_B_MIN
               : 0.0f;

    float lr, lg, lb;
    oklab_to_rgb(L, av, bv, &lr, &lg, &lb);

    *r = linear_to_srgb(clamp_f(lr));
    *g = linear_to_srgb(clamp_f(lg));
    *b = linear_to_srgb(clamp_f(lb));
}

/* Bayer ordered-dither quantisation.
 * Adds a deterministic per-position signed bias before quantising. */
static int quantize_bayer(int c, int x, int atrow, int atcol)
{
    int step  = (B[c].levels > 1) ? 256 / B[c].levels : 256;
    int threshold = bayer8[atrow & 7][atcol & 7]; /* 0–63 */
    int bias  = (threshold * step / 64) - step / 2;
    int dithered = clamp_byte(x + bias);
    return use_gamma ? quantize_gamma(c, dithered) : quantize(c, dithered);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Colour cube construction
 * ═══════════════════════════════════════════════════════════════════════════ */

/* Standard RGB uniform colour cube (original algorithm). */
static void make_color_cube_rgb(struct Colors *colors)
{
    int nr = B[0].levels;
    int ng = B[1].levels;
    int nb = B[2].levels;
    int mr = B[0].maxlev;
    int mg = B[1].maxlev;
    int mb = B[2].maxlev;
    int g, b;
    int i = 0;

    Rast_init_colors(colors);
    G_message(_("Creating RGB colour cube for output..."));

    for (b = 0; b < nb; b++) {
        G_percent(b, nb, 5);
        for (g = 0; g < ng; g++) {
            int blu = (mb > 0) ? b * 255 / mb : 0;
            int grn = (mg > 0) ? g * 255 / mg : 0;
            CELL i0 = i;
            CELL i1 = i + mr;
            Rast_add_c_color_rule(&i0, 0, grn, blu, &i1, 255, grn, blu, colors);
            i += nr;
        }
    }
    G_percent(nb, nb, 1);
}

/* Oklab perceptual colour cube.
 * Palette entries are equally spaced in Oklab space; each entry maps to the
 * nearest in-gamut sRGB colour after converting back from Oklab. */
static void make_color_cube_oklab(struct Colors *colors)
{
    int nL = B[0].levels;
    int na = B[1].levels;
    int nb = B[2].levels;
    int bi, ai;
    int i = 0;

    Rast_init_colors(colors);
    G_message(_("Creating Oklab perceptual colour cube for output..."));

    for (bi = 0; bi < nb; bi++) {
        G_percent(bi, nb, 5);
        for (ai = 0; ai < na; ai++) {
            float av = (na > 1)
                       ? (float)ai / (na - 1) * (OKLAB_A_MAX - OKLAB_A_MIN) + OKLAB_A_MIN
                       : 0.0f;
            float bv = (nb > 1)
                       ? (float)bi / (nb - 1) * (OKLAB_B_MAX - OKLAB_B_MIN) + OKLAB_B_MIN
                       : 0.0f;

            /* Compute sRGB at L=0 and L=1 to build a gradient rule.
             * The Oklab L axis is designed to be perceptually linear, so a
             * linear interpolation over L gives a visually smooth ramp. */
            float r0, g0, b0, r1, g1, b1;
            oklab_to_rgb(0.0f, av, bv, &r0, &g0, &b0);
            oklab_to_rgb(1.0f, av, bv, &r1, &g1, &b1);

            CELL i0 = i;
            CELL i1 = i + nL - 1;

            Rast_add_c_color_rule(
                &i0,
                linear_to_srgb(clamp_f(r0)),
                linear_to_srgb(clamp_f(g0)),
                linear_to_srgb(clamp_f(b0)),
                &i1,
                linear_to_srgb(clamp_f(r1)),
                linear_to_srgb(clamp_f(g1)),
                linear_to_srgb(clamp_f(b1)),
                colors);

            i += nL;
        }
    }
    G_percent(nb, nb, 1);
}

/* ═══════════════════════════════════════════════════════════════════════════
 * Adaptive palette-level estimation via histogram sampling
 * ═══════════════════════════════════════════════════════════════════════════ */

/* Fill hist[256] from band_idx by sampling every ~100th row. */
static void sample_histogram(int band_idx, int *hist, int nrows, int ncols)
{
    struct band *b = &B[band_idx];
    unsigned char *buf    = G_malloc(ncols);
    unsigned char *sink   = G_malloc(ncols);
    unsigned char *nullbuf = G_malloc(ncols);
    int step = (nrows > 100) ? nrows / 100 : 1;
    int row, col;

    memset(hist, 0, 256 * sizeof(int));

    for (row = 0; row < nrows; row += step) {
        unsigned char *arr[3];
        int j;
        for (j = 0; j < 3; j++)
            arr[j] = (j == band_idx) ? buf : sink;

        Rast_get_row_colors(b->file, row, &b->colors,
                            arr[0], arr[1], arr[2], nullbuf);

        for (col = 0; col < ncols; col++)
            if (!nullbuf[col])
                hist[buf[col]]++;
    }

    G_free(buf);
    G_free(sink);
    G_free(nullbuf);
}

/* Choose level count from histogram: ~4 source values per palette step,
 * clamped to [8, 256]. */
static int adaptive_levels(const int *hist)
{
    long total = 0;
    long cumsum;
    int p2 = 0, p98 = 255;
    long lo, hi;
    int i, range, levels;

    for (i = 0; i < 256; i++)
        total += hist[i];

    if (total == 0)
        return 32;

    lo = total * 2  / 100;
    hi = total * 98 / 100;
    cumsum = 0;

    for (i = 0; i < 256; i++) {
        cumsum += hist[i];
        if (cumsum <= lo)  p2  = i;
        if (cumsum <= hi)  p98 = i;
    }

    range  = p98 - p2 + 1;
    levels = (range + 3) / 4; /* ceil(range/4) */
    if (levels < 8)   levels = 8;
    if (levels > 256) levels = 256;
    return levels;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * main
 * ═══════════════════════════════════════════════════════════════════════════ */

int main(int argc, char **argv)
{
    struct GModule *module;
    struct Option *opt_out;
    struct Option *opt_lev;
    struct Option *opt_null;
    struct Flag   *flg_p;  /* palette (indexed) mode */
    struct Flag   *flg_d;  /* Floyd-Steinberg dither */
    struct Flag   *flg_b;  /* Bayer ordered dither */
    struct Flag   *flg_c;  /* closest-colour rounding */
    struct Flag   *flg_g;  /* gamma-correct quantisation */
    struct Flag   *flg_k;  /* Oklab perceptual quantisation */

    int dither, bayer, palette_mode;
    char *out_name;
    int atrow, atcol;
    struct Cell_head window;
    int i, j;
    struct History history;
    int do_null_fill;
    int null_fill;
    int auto_levels;
    int levels;

    /* Per-band null buffers and combined null mask */
    unsigned char *nulls[3];
    unsigned char *combined_nulls;

    /* Palette mode output */
    int    out_file_pal  = -1;
    CELL  *out_arr_pal   = NULL;
    struct Colors out_colors;

    /* True-colour mode output (three FCELL maps) */
    char   out_name_r[GNAME_MAX], out_name_g[GNAME_MAX], out_name_b[GNAME_MAX];
    int    out_file_r = -1, out_file_g = -1, out_file_b = -1;
    FCELL *out_arr_r  = NULL, *out_arr_g = NULL, *out_arr_b = NULL;

    G_gisinit(argv[0]);

    module = G_define_module();
    G_add_keyword(_("raster"));
    G_add_keyword(_("composite"));
    G_add_keyword("RGB");
    module->description =
        _("Combines red, green and blue raster maps into a composite. "
          "Default output: three true-colour FCELL maps (output.r/g/b). "
          "Use -p for a legacy palette-indexed CELL map.");

    for (i = 0; i < 3; i++) {
        struct Option *opt;
        char buff[80];

        B[i].opt_name = opt = G_define_standard_option(G_OPT_R_INPUT);
        snprintf(buff, sizeof(buff), "%s", color_names[i]);
        opt->key = G_store(buff);
        opt->answer = NULL;
        snprintf(buff, sizeof(buff),
                 _("Name of raster map to be used for <%s>"), color_names[i]);
        opt->description = G_store(buff);
    }

    opt_lev = G_define_option();
    opt_lev->key         = "levels";
    opt_lev->type        = TYPE_STRING; /* accepts integers and "auto" */
    opt_lev->required    = NO;
    opt_lev->answer      = "64";
    opt_lev->description =
        _("Number of levels per component for palette mode (1-256 or \"auto\")");
    opt_lev->guisection  = _("Levels");

    for (i = 0; i < 3; i++) {
        struct Option *opt;
        char buff[80];

        B[i].opt_levels = opt = G_define_option();
        snprintf(buff, sizeof(buff), "level_%s", color_names[i]);
        opt->key = G_store(buff);
        opt->type        = TYPE_INTEGER;
        opt->required    = NO;
        opt->options     = "1-256";
        snprintf(buff, sizeof(buff),
                 _("Number of levels for <%s> (palette mode, overrides levels=)"),
                 color_names[i]);
        opt->description = G_store(buff);
        opt->guisection  = _("Levels");
    }

    opt_out = G_define_standard_option(G_OPT_R_OUTPUT);

    opt_null = G_define_option();
    opt_null->key         = "null_value";
    opt_null->type        = TYPE_INTEGER;
    opt_null->required    = NO;
    opt_null->options     = "0-255";
    opt_null->description =
        _("Fill value (0-255) for NULL pixels; if omitted, NULLs are propagated");

    flg_p = G_define_flag();
    flg_p->key         = 'p';
    flg_p->description =
        _("Palette mode: produce a single palette-indexed CELL map (legacy behaviour)");

    flg_d = G_define_flag();
    flg_d->key         = 'd';
    flg_d->description = _("Floyd-Steinberg dithering (palette mode)");

    flg_b = G_define_flag();
    flg_b->key         = 'b';
    flg_b->description =
        _("Bayer ordered dithering (palette mode; parallelisable, no streak artefacts)");

    flg_c = G_define_flag();
    flg_c->key         = 'c';
    flg_c->description = _("Use closest-colour rounding (palette mode)");

    flg_g = G_define_flag();
    flg_g->key         = 'g';
    flg_g->description =
        _("Gamma-correct quantisation: linearise to linear light before quantising "
          "(palette mode)");

    flg_k = G_define_flag();
    flg_k->key         = 'k';
    flg_k->description =
        _("Oklab perceptual quantisation: uniform palette in Oklab colour space "
          "(palette mode; implies linear-light processing)");

    if (G_parser(argc, argv))
        exit(EXIT_FAILURE);

    /* ── Parse flags ───────────────────────────────────────────────────── */
    palette_mode = flg_p->answer;
    dither       = flg_d->answer;
    bayer        = flg_b->answer;
    closest      = flg_c->answer;
    use_gamma    = flg_g->answer;
    use_oklab    = flg_k->answer;

    if (dither && bayer)
        G_fatal_error(_("Flags -d (Floyd-Steinberg) and -b (Bayer) are mutually exclusive"));
    if (use_oklab && use_gamma)
        G_warning(_("Flag -g is redundant with -k: Oklab already operates in linear light"));
    if (!palette_mode && (dither || bayer || flg_c->answer || use_gamma || use_oklab))
        G_warning(_("Dithering and quantisation flags (-d/-b/-c/-g/-k) only affect "
                    "palette mode (-p)"));

    /* ── Parse null-fill ───────────────────────────────────────────────── */
    do_null_fill = (opt_null->answer != NULL);
    null_fill    = do_null_fill ? atoi(opt_null->answer) : 0;

    /* ── Parse levels (global default; may be "auto") ──────────────────── */
    if (strcmp(opt_lev->answer, "auto") == 0) {
        auto_levels = 1;
        levels      = 64; /* placeholder; replaced per-band below */
    }
    else {
        auto_levels = 0;
        levels = atoi(opt_lev->answer);
        if (levels < 1 || levels > 256)
            G_fatal_error(_("levels must be 1-256 or \"auto\""));
    }

    /* ── Validate output name length for band-suffix appending ─────────── */
    out_name = opt_out->answer;
    if (!palette_mode && strlen(out_name) > (size_t)(GNAME_MAX - 3))
        G_fatal_error(_("Output map name <%s> is too long; must leave room for "
                        ".r/.g/.b suffixes"), out_name);

    /* ── Read window ───────────────────────────────────────────────────── */
    G_get_window(&window);

    dummy         = G_malloc(window.cols);
    combined_nulls = G_malloc(window.cols);
    for (i = 0; i < 3; i++)
        nulls[i] = G_malloc(window.cols);

    /* ── Open input bands ──────────────────────────────────────────────── */
    for (i = 0; i < 3; i++) {
        struct band *b = &B[i];

        b->name = b->opt_name->answer;
        b->file = Rast_open_old(b->name, "");
        b->type = Rast_get_map_type(b->file);
        b->size = Rast_cell_size(b->type);

        if (Rast_read_colors(b->name, "", &b->colors) == -1)
            G_fatal_error(_("Unable to read colour file of raster map <%s>"),
                          b->name);

        for (j = 0; j < 3; j++)
            b->array[j] = (i == j) ? G_malloc(window.cols) : dummy;

        /* Per-band level assignment (band-specific option overrides global) */
        if (b->opt_levels->answer) {
            b->levels = atoi(b->opt_levels->answer);
        }
        else if (auto_levels && palette_mode) {
            b->levels = levels; /* placeholder; updated after histogram pass */
        }
        else {
            b->levels = levels;
        }
        b->maxlev = b->levels - 1;
        b->offset = (b->maxlev > 0) ? 128 / b->maxlev : 0;

        if (dither)
            for (j = 0; j < 2; j++)
                b->floyd[j] = G_calloc(window.cols + 2, sizeof(short));
    }

    /* ── Adaptive level estimation (palette mode + auto only) ──────────── */
    if (auto_levels && palette_mode) {
        G_message(_("Sampling histograms for adaptive level selection..."));
        for (i = 0; i < 3; i++) {
            if (!B[i].opt_levels->answer) {
                int hist[256];
                sample_histogram(i, hist, window.rows, window.cols);
                B[i].levels = adaptive_levels(hist);
                B[i].maxlev = B[i].levels - 1;
                B[i].offset = (B[i].maxlev > 0) ? 128 / B[i].maxlev : 0;
                G_verbose_message(_("Band %s: adaptive levels = %d"),
                                  color_names[i], B[i].levels);
            }
        }
    }

    /* ── Open output ───────────────────────────────────────────────────── */
    if (palette_mode) {
        out_file_pal = Rast_open_c_new(out_name);
        out_arr_pal  = Rast_allocate_c_buf();
        if (use_oklab)
            make_color_cube_oklab(&out_colors);
        else
            make_color_cube_rgb(&out_colors);
    }
    else {
        snprintf(out_name_r, sizeof(out_name_r), "%s.r", out_name);
        snprintf(out_name_g, sizeof(out_name_g), "%s.g", out_name);
        snprintf(out_name_b, sizeof(out_name_b), "%s.b", out_name);

        out_file_r = Rast_open_fp_new(out_name_r);
        out_file_g = Rast_open_fp_new(out_name_g);
        out_file_b = Rast_open_fp_new(out_name_b);

        out_arr_r = Rast_allocate_f_buf();
        out_arr_g = Rast_allocate_f_buf();
        out_arr_b = Rast_allocate_f_buf();
    }

    G_message(_("Writing raster map <%s>..."), out_name);

    /* ── Main processing loop ──────────────────────────────────────────── */
    for (atrow = 0; atrow < window.rows; atrow++) {
        G_percent(atrow, window.rows, 2);

        /* Read colour values from all three input bands (per-band nulls) */
        for (i = 0; i < 3; i++) {
            struct band *b = &B[i];
            Rast_get_row_colors(b->file, atrow, &b->colors,
                                b->array[0], b->array[1], b->array[2],
                                nulls[i]);
            if (dither) {
                short *tmp = b->floyd[0];
                b->floyd[0] = b->floyd[1];
                for (atcol = 0; atcol < window.cols + 2; atcol++)
                    tmp[atcol] = 0;
                b->floyd[1] = tmp;
            }
        }

        /* Build combined null mask */
        for (atcol = 0; atcol < window.cols; atcol++)
            combined_nulls[atcol] =
                nulls[0][atcol] | nulls[1][atcol] | nulls[2][atcol];

        /* ── True-colour (FCELL) path ──────────────────────────────────── */
        if (!palette_mode) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (atcol = 0; atcol < window.cols; atcol++) {
                /* Per-band null handling: each output map tracks its band */
                int n0 = nulls[0][atcol];
                int n1 = nulls[1][atcol];
                int n2 = nulls[2][atcol];

                if (n0 && !do_null_fill)
                    Rast_set_f_null_value(&out_arr_r[atcol], 1);
                else
                    out_arr_r[atcol] =
                        (FCELL)(n0 ? null_fill : (int)B[0].array[0][atcol]);

                if (n1 && !do_null_fill)
                    Rast_set_f_null_value(&out_arr_g[atcol], 1);
                else
                    out_arr_g[atcol] =
                        (FCELL)(n1 ? null_fill : (int)B[1].array[1][atcol]);

                if (n2 && !do_null_fill)
                    Rast_set_f_null_value(&out_arr_b[atcol], 1);
                else
                    out_arr_b[atcol] =
                        (FCELL)(n2 ? null_fill : (int)B[2].array[2][atcol]);
            }

            Rast_put_f_row(out_file_r, out_arr_r);
            Rast_put_f_row(out_file_g, out_arr_g);
            Rast_put_f_row(out_file_b, out_arr_b);
        }

        /* ── Palette (indexed CELL) path ───────────────────────────────── */
        else {
            /* Floyd-Steinberg is inherently serial (forward error diffusion) */
            if (dither) {
                for (atcol = 0; atcol < window.cols; atcol++) {
                    if (combined_nulls[atcol] && !do_null_fill) {
                        Rast_set_c_null_value(&out_arr_pal[atcol], 1);
                        continue;
                    }

                    if (use_oklab) {
                        /* Collect dither-adjusted sRGB values for all bands */
                        int vr = combined_nulls[atcol] ? null_fill
                             : clamp_byte((int)B[0].array[0][atcol]
                                          + B[0].floyd[0][atcol + 1] / 16);
                        int vg = combined_nulls[atcol] ? null_fill
                             : clamp_byte((int)B[1].array[1][atcol]
                                          + B[1].floyd[0][atcol + 1] / 16);
                        int vb = combined_nulls[atcol] ? null_fill
                             : clamp_byte((int)B[2].array[2][atcol]
                                          + B[2].floyd[0][atcol + 1] / 16);

                        int idx = quantize_oklab(vr, vg, vb);

                        /* Recover reconstructed sRGB to compute per-channel error */
                        int rv, rg, rb;
                        oklab_idx_to_srgb(idx, &rv, &rg, &rb);

                        int d0 = vr - rv;
                        int d1 = vg - rg;
                        int d2 = vb - rb;

                        B[0].floyd[0][atcol + 2] += 7 * d0;
                        B[0].floyd[1][atcol + 0] += 3 * d0;
                        B[0].floyd[1][atcol + 1] += 5 * d0;
                        B[0].floyd[1][atcol + 2] += 1 * d0;

                        B[1].floyd[0][atcol + 2] += 7 * d1;
                        B[1].floyd[1][atcol + 0] += 3 * d1;
                        B[1].floyd[1][atcol + 1] += 5 * d1;
                        B[1].floyd[1][atcol + 2] += 1 * d1;

                        B[2].floyd[0][atcol + 2] += 7 * d2;
                        B[2].floyd[1][atcol + 0] += 3 * d2;
                        B[2].floyd[1][atcol + 1] += 5 * d2;
                        B[2].floyd[1][atcol + 2] += 1 * d2;

                        out_arr_pal[atcol] = (CELL)idx;
                    }
                    else {
                        int val[3];
                        for (i = 0; i < 3; i++) {
                            int v = combined_nulls[atcol]
                                    ? null_fill
                                    : (int)B[i].array[i][atcol];
                            int r, w, d;
                            v = clamp_byte(v + B[i].floyd[0][atcol + 1] / 16);
                            r = use_gamma ? quantize_gamma(i, v) : quantize(i, v);
                            w = (B[i].maxlev > 0) ? r * 255 / B[i].maxlev : 0;
                            d = v - w;
                            B[i].floyd[0][atcol + 2] += 7 * d;
                            B[i].floyd[1][atcol + 0] += 3 * d;
                            B[i].floyd[1][atcol + 1] += 5 * d;
                            B[i].floyd[1][atcol + 2] += 1 * d;
                            val[i] = r;
                        }
                        out_arr_pal[atcol] =
                            (CELL)(val[2] * B[1].levels + val[1]) * B[0].levels + val[0];
                    }
                }
            }
            else {
                /* Non-dithering or Bayer: parallelisable */
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
                for (atcol = 0; atcol < window.cols; atcol++) {
                    if (combined_nulls[atcol] && !do_null_fill) {
                        Rast_set_c_null_value(&out_arr_pal[atcol], 1);
                        continue;
                    }

                    if (use_oklab) {
                        int vr = combined_nulls[atcol] ? null_fill : (int)B[0].array[0][atcol];
                        int vg = combined_nulls[atcol] ? null_fill : (int)B[1].array[1][atcol];
                        int vb = combined_nulls[atcol] ? null_fill : (int)B[2].array[2][atcol];
                        out_arr_pal[atcol] = (CELL)quantize_oklab(vr, vg, vb);
                    }
                    else {
                        int qv[3];
                        int ci;
                        for (ci = 0; ci < 3; ci++) {
                            int x = combined_nulls[atcol] ? null_fill
                                                          : (int)B[ci].array[ci][atcol];
                            qv[ci] = bayer  ? quantize_bayer(ci, x, atrow, atcol)
                                   : use_gamma ? quantize_gamma(ci, x)
                                   :             quantize(ci, x);
                        }
                        out_arr_pal[atcol] =
                            (CELL)(qv[2] * B[1].levels + qv[1]) * B[0].levels + qv[0];
                    }
                }
            }

            Rast_put_row(out_file_pal, out_arr_pal, CELL_TYPE);
        }
    }
    G_percent(window.rows, window.rows, 1);

    /* ── Close input files ─────────────────────────────────────────────── */
    for (i = 0; i < 3; i++)
        Rast_close(B[i].file);

    /* ── Write output and metadata ─────────────────────────────────────── */
    if (palette_mode) {
        Rast_close(out_file_pal);
        Rast_write_colors(out_name, G_mapset(), &out_colors);
        Rast_short_history(out_name, "raster", &history);
        Rast_command_history(&history);
        Rast_write_history(out_name, &history);
        G_done_msg(_("Raster map <%s> created."), out_name);
    }
    else {
        const char *tc_names[3] = {out_name_r, out_name_g, out_name_b};
        int        *tc_files[3] = {&out_file_r, &out_file_g, &out_file_b};

        for (i = 0; i < 3; i++) {
            struct Colors tcol;
            FCELL f0 = 0.0f, f255 = 255.0f;

            Rast_close(*tc_files[i]);

            /* Greyscale ramp 0–255 for each component map */
            Rast_init_colors(&tcol);
            Rast_add_f_color_rule(&f0, 0, 0, 0, &f255, 255, 255, 255, &tcol);
            Rast_write_colors(tc_names[i], G_mapset(), &tcol);
            Rast_free_colors(&tcol);

            Rast_short_history(tc_names[i], "raster", &history);
            Rast_command_history(&history);
            Rast_write_history(tc_names[i], &history);
        }
        G_done_msg(_("True-colour maps <%s.r>, <%s.g>, <%s.b> created."),
                   out_name, out_name, out_name);
    }

    exit(EXIT_SUCCESS);
}
