## DESCRIPTION

*r.composite* combines three separate raster maps (red, green, blue
channels) into a composite RGB output.

**Default (true-colour) mode** writes three FCELL raster maps named
`output.r`, `output.g`, and `output.b`, each holding the raw 0–255
component values from the input maps' colour tables. This is lossless
and directly compatible with downstream GRASS modules such as
*[i.cluster](i.cluster.md)*, *[r.mapcalc](r.mapcalc.md)*, and
*[d.rgb](d.rgb.md)*.

**Palette (legacy) mode** (`-p`) produces a single palette-indexed CELL
map using colour-space quantisation, identical in concept to the
original algorithm. All dithering and quantisation flags apply only in
this mode.

## OPTIONS

### Output mode

| Flag | Description |
|------|-------------|
| (default) | True-colour FCELL output: `output.r`, `output.g`, `output.b` |
| `-p` | Palette mode: single indexed CELL map |

### Quantisation quality (palette mode `-p`)

| Flag/Parameter | Description |
|----------------|-------------|
| `levels=N` | Intensity levels per component (1–256, default **64**, or `auto`) |
| `level_red=`, `level_green=`, `level_blue=` | Per-channel level overrides |
| `-g` | Gamma-correct: linearise to linear light before quantising |
| `-k` | Oklab: quantise in a perceptually uniform colour space |
| `-c` | Closest-colour rounding (alternative to floor rounding) |

### Dithering (palette mode `-p`)

| Flag | Description |
|------|-------------|
| `-d` | Floyd-Steinberg error-diffusion dithering (sequential, smooth gradients) |
| `-b` | Bayer 8×8 ordered dithering (deterministic, parallelisable, no streak artefacts) |

### NULL handling

| Parameter | Description |
|-----------|-------------|
| `null_value=N` | Fill value (0–255) for NULL pixels; omit to propagate NULLs |

In true-colour mode, NULLs are tracked per input band: if only the red
input has a NULL at a pixel, only `output.r` gets a NULL there — the
green and blue output maps retain their valid values.

## NOTES

### Quantisation levels (`levels=`)

The default **64 levels** per component yields 64³ = 262,144 palette
entries (equivalent to ~18 bits). Use `levels=auto` to let the module
sample each band's histogram and choose a level count proportional to
the effective data range (approximately 4 source values per palette
step, clamped to 8–256).

### Gamma-correct quantisation (`-g`)

Standard sRGB values are gamma-encoded. Quantising directly in sRGB
distributes palette steps unevenly: too many in highlights, too few in
dark tones where human vision is most sensitive. The `-g` flag removes
the sRGB gamma transfer function (IEC 61966-2-1) before quantising and
re-applies it to the palette representative colours, distributing
palette steps evenly in linear light.

### Oklab perceptual quantisation (`-k`)

Oklab is a perceptually uniform colour space designed so that equal
distances correspond to equal perceived colour differences. The `-k`
flag converts each pixel to Oklab before quantising, then builds the
colour table from the Oklab grid back-converted to sRGB. This reduces
visible banding at the same level count compared to sRGB or even
gamma-correct quantisation. Oklab always operates in linear light, so
`-g` is redundant when `-k` is active.

### Bayer vs Floyd-Steinberg dithering

Floyd-Steinberg (`-d`) diffuses quantisation error to neighbouring
pixels, producing smooth gradients but sequential (non-parallelisable)
processing and diagonal streak artefacts on edges. Bayer ordered
dithering (`-b`) uses a deterministic 8×8 threshold matrix: each pixel
is processed independently (fully parallelisable), artefacts are
spatially structured rather than streaky, and the result is stable
across tiled processing.

## EXAMPLES

### True-colour composite (default)

```sh
g.region raster=lsat7_2002_10
r.composite blue=lsat7_2002_10 green=lsat7_2002_20 red=lsat7_2002_30 \
            output=lsat7_2002_rgb
# Creates: lsat7_2002_rgb.r, lsat7_2002_rgb.g, lsat7_2002_rgb.b
```

### Display with d.rgb

```sh
d.rgb red=lsat7_2002_rgb.r green=lsat7_2002_rgb.g blue=lsat7_2002_rgb.b
```

### Legacy palette mode (backward-compatible)

```sh
r.composite -p blue=lsat7_2002_10 green=lsat7_2002_20 red=lsat7_2002_30 \
               output=lsat7_2002_pal
```

### Perceptual palette with Bayer dithering

```sh
r.composite -pkb red=elevation.r green=elevation.g blue=elevation.b \
               output=elev.composite
```

### Adaptive levels with Oklab quantisation

```sh
r.composite -pk levels=auto red=s2_b4 green=s2_b3 blue=s2_b2 \
               output=s2_composite
```

### Fill NULLs instead of propagating them

```sh
r.composite red=s2_b4 green=s2_b3 blue=s2_b2 \
            output=s2_composite null_value=0
```

## SEE ALSO

*[d.rast](d.rast.md), [d.rgb](d.rgb.md), [r.blend](r.blend.md),
[r.colors](r.colors.md), [r.rgb](r.rgb.md)*

*[Wikipedia: Floyd-Steinberg dithering](https://en.wikipedia.org/wiki/Floyd-Steinberg_dithering)*
*[Wikipedia: Ordered dithering](https://en.wikipedia.org/wiki/Ordered_dithering)*
*[Oklab colour space](https://bottosson.github.io/posts/oklab/)*

## AUTHORS

Glynn Clements
Modern GIS quality improvements (true-colour output, Oklab, Bayer
dithering, gamma correction, adaptive levels, per-band NULL handling,
OpenMP)
