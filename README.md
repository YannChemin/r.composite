# r.composite

GRASS GIS module that combines red, green, and blue raster maps into a
composite RGB output.

## Output modes

**Default — true-colour FCELL** writes three lossless floating-point maps:

```
output.r   output.g   output.b
```

Each map holds the raw 0–255 component values from the input maps'
colour tables, making them directly usable with `d.rgb`, `r.mapcalc`,
`i.cluster`, and any other GRASS module that works on single-band
rasters.

**Palette mode** (`-p`) produces a single palette-indexed CELL map
using colour-space quantisation — the original algorithm, preserved for
backward compatibility.

## Quick start

```sh
# True-colour composite (default)
g.region raster=lsat7_2002_10
r.composite red=lsat7_2002_30 green=lsat7_2002_20 blue=lsat7_2002_10 \
            output=lsat7_rgb
# → lsat7_rgb.r, lsat7_rgb.g, lsat7_rgb.b

# Display
d.rgb red=lsat7_rgb.r green=lsat7_rgb.g blue=lsat7_rgb.b

# Legacy single-map output
r.composite -p red=lsat7_2002_30 green=lsat7_2002_20 blue=lsat7_2002_10 \
               output=lsat7_pal
```

## Flags

| Flag | Description |
|------|-------------|
| `-p` | Palette mode: single indexed CELL map |
| `-k` | Oklab perceptual quantisation (palette mode) |
| `-g` | Gamma-correct quantisation in linear light (palette mode) |
| `-b` | Bayer 8×8 ordered dithering (palette mode) |
| `-d` | Floyd-Steinberg error-diffusion dithering (palette mode) |
| `-c` | Closest-colour rounding (palette mode) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `red=` | — | Input raster for red channel |
| `green=` | — | Input raster for green channel |
| `blue=` | — | Input raster for blue channel |
| `output=` | — | Output map name (base name in true-colour mode) |
| `levels=` | `64` | Levels per component for palette mode (1–256 or `auto`) |
| `level_red=` | — | Per-channel level override |
| `level_green=` | — | Per-channel level override |
| `level_blue=` | — | Per-channel level override |
| `null_value=` | — | Fill NULLs with this value (0–255) instead of propagating |

## Palette quality options

### `levels=auto`

Samples each band's histogram to choose a level count proportional to
the effective data range (~4 source values per palette step, clamped to
8–256). Useful for narrow-range inputs that don't need the full 64
levels.

### `-g` gamma-correct quantisation

sRGB values are gamma-encoded. Quantising directly in sRGB wastes
palette steps in highlights and crowds them in shadows where human
vision is most sensitive. `-g` linearises to linear light (IEC
61966-2-1) before quantising and re-applies gamma to the palette
representative colours.

### `-k` Oklab perceptual quantisation

[Oklab](https://bottosson.github.io/posts/oklab/) is a perceptually
uniform colour space where equal distances correspond to equal perceived
colour differences. `-k` converts each pixel to Oklab before quantising
and builds the colour table from the Oklab grid back-converted to sRGB.
Produces less visible banding than sRGB or gamma-correct quantisation at
the same level count. Implies linear-light processing, so `-g` is
redundant alongside `-k`.

### `-b` Bayer vs `-d` Floyd-Steinberg dithering

| | Floyd-Steinberg (`-d`) | Bayer (`-b`) |
|---|---|---|
| Error propagation | Forward-diffused to neighbours | None (threshold matrix) |
| Processing order | Sequential (serial) | Independent per pixel (parallelisable) |
| Artefact character | Smooth gradients, diagonal streaks | Structured grain, no streaks |
| Tiled processing | Seams at tile boundaries | Seamless |

## NULL handling

In **true-colour mode**, NULLs are tracked per input band independently.
If only the red band is NULL at a pixel, only `output.r` gets a NULL —
`output.g` and `output.b` keep their valid values.

In **palette mode**, any NULL in any band produces a NULL output cell
(same as the original behaviour).

`null_value=N` overrides both modes: NULLs are filled with the given
value rather than propagated.

## Building

r.composite is part of the GRASS GIS source tree:

```sh
cd grass
make MODULE_TOPDIR=$(pwd) -C raster/r.composite
```

Requires a GRASS build environment with `$(MATHLIB)` (`-lm`) for `powf`
and `cbrtf`, and optionally OpenMP for multi-threaded processing.

## Tests

```sh
cd raster/r.composite
python -m pytest tests/ -v
```

## See also

`d.rgb`, `d.rast`, `r.blend`, `r.colors`, `r.rgb`

## Authors

Glynn Clements (original)
Modern improvements: true-colour output, Oklab, Bayer dithering, gamma
correction, adaptive levels, per-band NULL handling, OpenMP
