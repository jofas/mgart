![mgart](static/icon.svg)

# Mgart

**M**achine **G**enerated **Art**, short Mgart and pronounced 
"em-gart" is a rust crate and CLI application for generating 
algorithmic art.

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
