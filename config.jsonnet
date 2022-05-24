local LCH(l, c, h) = {
  type: 'lch',
  l: l,
  c: c,
  h: h,
};

local RGBA(r, g, b, a) = {
  type: 'rgba',
  r: r,
  g: g,
  b: b,
  a: a,
};

[
  {
    width: 3400,
    height: 800,
    command: 'color-map-1d',
    filename: 'lch.png',
    color_map: {
      method: 'linear',
      color_space: 'lch',
      map: [
        //LCH(100, 0, 'NaN'),
        LCH(53.24, 104.55, 40),
        LCH(87.73, 119.78, 136.02),
        //LCH(32.3, 133.81, 306.28),
        //LCH(0, 0, 'NaN'),
      ],
    },
  },
  {
    width: 3400,
    height: 800,
    command: 'color-map-1d',
    filename: 'rgba.png',
    color_map: {
      method: 'linear',
      color_space: 'rgba',
      map: [
        //LCH(100, 0, 'NaN'),
        LCH(53.24, 104.55, 40),
        LCH(87.73, 119.78, 136.02),
        //LCH(32.3, 133.81, 306.28),
        //LCH(0, 0, 'NaN'),
      ],
    },
  },
  /*
  {
    width: 800,
    height: 900,
    filename: '1.png',
    zoom: 5.0e11,
    zpx: 0.3750001200618655,
    zpy: -0.2166393884377127,
    iter: 1000000,
    color_map: [
      std.parseHex('FFFFFFFF'),
      std.parseHex('FF0000FF'),
      std.parseHex('00FF00FF'),
      std.parseHex('0000FFFF'),
      std.parseHex('000000FF'),
    ],
  },
  */
]
