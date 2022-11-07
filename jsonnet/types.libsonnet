{
  complex: {
    cartesian(re, im): {
      re: re,
      im: im,
    },
    polar(r, theta): {
      r: r,
      theta: theta
    },
  },
  color: {
    lch(l, c, h): {
      l: l,
      c: c,
      h: h,
      type: 'lch',
    },
    rgb(r, g, b): {
      r: r,
      g: g,
      b: b,
      type: 'rgb',
    },
  },
  sampler: {
    uniform_polar(r): {
      type: 'uniform_polar',
      r: r,
    },
    kde(population, h, p_min): {
      type: 'kde',
      population: population,
      h: h,
      p_min: p_min,
    },
  },
  gradient: {
    linear(factor=1): {
      type: 'linear',
      factor: factor,
    },
  },
  color_map_1d(gradient, map): {
    gradient: gradient,
    map: map,
  },
  post_processing: {
    normalize: {
      process: 'normalize',
    },
    clahe(contrast_limit, bin_count, tile_size_x, tile_size_y): {
      process: 'clahe',
      contrast_limit: contrast_limit,
      bin_count: bin_count,
      tile_size_x: tile_size_x,
      tile_size_y: tile_size_y,
    },
    smoothing: {
      non_local_means(n, window_size, h): {
        process: 'smoothing',
        type: 'non_local_means',
        n: n,
        window_size: window_size,
        h: h,
      },
    },
  },
}
