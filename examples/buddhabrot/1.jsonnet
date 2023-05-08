local constants = import "jsonnet/constants.libsonnet";
local types = import "jsonnet/types.libsonnet";
local algorithms = import "jsonnet/algorithms.libsonnet";

[
  algorithms.buddhabrot(
    filename="1.webp",
    width=3840,
    height=3840,
    center=types.complex.cartesian(0, -0.55),
    zoom=4,
    rotation=90,
    sample_count=5 * constants.numbers["100m"],
    iter=20000,
    sampler=types.sampler.kernel_density_estimation(
      weighted=true,
      kernel=types.sampler.uniform(2e-2),
      population=10 * constants.numbers["100m"],
      p_min=0.01,
      pre_sampler=types.sampler.uniform_polar(4),
    ),
    post_processing=[
      types.post_processing.normalize,
      types.post_processing.clamp_and_normalize(min=0.05),
      types.post_processing.clahe(100, 256, 192, 192),
      types.post_processing.smoothing.non_local_means(7, 25, 5e-4),
      types.post_processing.gradient(types.gradient.exp(1 / 2)),
    ],
    color_map=types.color_map_1d(
      gradient=types.gradient.linear(),
      map=[
        constants.colors.black,
        constants.colors.white,
      ],
    ),
  ),
]
