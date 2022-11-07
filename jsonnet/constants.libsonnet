local types = import "types.libsonnet";

local colors = {
  black: types.color.rgb(0, 0, 0),
  white: types.color.rgb(255, 255, 255),
  red: types.color.rgb(255, 0, 0),
  green: types.color.rgb(0, 255, 0),
  blue: types.color.rgb(0, 0, 255),
};

{
  numbers: {
    "1k": 1000,
    "1m": 1000000,
    "10m": 10000000,
    "100m": 100000000,
    "1bn": 1000000000,
    origin: types.complex.cartesian(0, 0),
  },
  colors: colors,
  color_maps: {
    greyscale: types.color_map_1d(
      gradient=types.gradient.linear(),
      map=[colors.black, colors.white],
    ),
  }
}
