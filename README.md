![mgart](static/icon.svg)

# Mgart

**M**achine **G**enerated **Art**, short Mgart and pronounced 
"em-gart" is a rust crate and CLI application for generating 
algorithmic art.

## Table of contents

<!--ts-->
<!--te-->

## Install

### Cargo

Note that you need to have the rust toolchain installed on your
computer if you want to install Mgart using `cargo`.

Mgart is distributed via [crates.io](https://crates.io) and can be 
installed with:

```bash
cargo install mgart
```

If you'd like to install a specific version of Mgart, use the
`--version` flag:

```bash
cargo install --version $VERSION
```

If you have an old version of Mgart already installed and wish to 
update it to the newest version, use the `--force` flag:

```bash
cargo install --force mgart
```


## Example

If you run `mgart -h`, you will find that Mgart takes a single file as
input argument.
The content of the file contains the configuration for the art
pieces you wish to generate.
Currently, Mgart supports Json and [Jsonnet](https://jsonnet.org/) 
input files.
The input file contains an array of `algorithm` objects where each one
describes an artwork you would like to create.
This is an example Json file that would generate a single artwork, a 
rendering of a [buddhabrot](https://en.wikipedia.org/wiki/Buddhabrot):

```json
[
   {
      "algorithm": "buddhabrot",
      "filename": "buddhabrot_20000.png",
      "width": 3800,
      "height": 2000,
      "center": {
         "im": -0.35,
         "re": 0
      },
      "zoom": 0.45,
      "exponent": 2,
      "iter": 20000,
      "rotation": 90,
      "sample_count": 500000000,
      "sampler": {
         "r": 3,
         "type": "uniform_polar"
      },
      "post_processing": [
         {
            "process": "normalize"
         },
         {
            "bin_count": 256,
            "contrast_limit": 500,
            "process": "clahe",
            "tile_size_x": 380,
            "tile_size_y": 200
         },
         {
            "h": 0.0005,
            "n": 7,
            "process": "smoothing",
            "type": "non_local_means",
            "window_size": 21
         }
      ],
      "color_map": {
         "gradient": {
            "factor": 1,
            "type": "linear"
         },
         "map": [
            {
               "b": 0,
               "g": 0,
               "r": 0,
               "type": "rgb"
            },
            {
               "b": 255,
               "g": 255,
               "r": 255,
               "type": "rgb"
            }
         ]
      }
   }
]
```

You simply run `mgart $FILE` and the resulting image looks like:

![Buddhabrot](examples/buddhabrot/greyscale/buddhabrot_20000.png)

You can find more example artworks and their configuration in the 
`examples/` folder.


## Supported Algorithms

Below you will find a list of algorithms either already supported by
Mgart or planned to be supported by a future release.

### Fractals

Algorithms for creating various types of fractal art.

#### Mandelbrot and Julia Sets

#### Buddhabrot

* [x] Buddhabrot

* [ ] Anti-Buddhabrot

* [ ] Nebulabrot

#### Other

* [ ] Fractal Flames

* [ ] Newton Fractals

* [ ] Strange Attractors

* [ ] L-Systems

### AI Art

* [ ] DeepDream-like filter

* [ ] Text-to-image

## Goals

* API-first + readability over performance

* aesthetic value over correctness

* fast algorithms over fast hardware
