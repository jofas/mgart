local constants = import "constants.libsonnet";
local types = import "types.libsonnet";

{
  buddhabrot(
    filename="buddhabrot.png",
    width=1920,
    height=1080,
    center=constants.numbers.origin,
    zoom=0.5,
    rotation=null,
    iter=20,
    sample_count=constants.numbers["100m"],
    exponent=2,
    sampler=types.sampler.uniform_polar(2),
    post_processing=[],
    color_map=constants.color_maps.greyscale,
  ): {
    width: width,
    height: height,
    command: 'buddhabrot',
    center: center,
    zoom: zoom,
    rotation: rotation,
    filename: filename,
    iter: iter,
    sample_count: sample_count,
    exponent: exponent,
    sampler: sampler,
    post_processing: post_processing,
    color_map: color_map,
  },
}
