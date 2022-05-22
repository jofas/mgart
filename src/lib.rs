use cgmath::Vector3;

pub struct ColorMap1D {
  colors: Vec<[u8; 4]>,
}

impl ColorMap1D {
  pub fn new(colors: Vec<[u8; 4]>) -> Self {
    let colors = if colors.len() >= 2 {
      colors
    } else if colors.len() == 1 {
      vec![[0, 0, 0, 0], colors[0]]
    } else {
      vec![[0, 0, 0, 0], [255, 255, 255, 255]]
    };

    Self { colors }
  }

  pub fn value(&self, x: f32) -> [u8; 4] {
    let x = 0.0_f32.max(0.99_f32.min(x));

    let interval = x * (self.colors.len() - 1) as f32;
    let pos = interval.fract();

    let c1 = self.colors[interval as usize];
    let c2 = self.colors[interval as usize + 1];

    let v1: Vector3<f32> = self.color_to_vec3(&c1);
    let v2: Vector3<f32> = self.color_to_vec3(&c2);

    let res = (v2 - v1) * pos + v1;

    [res.x.abs() as u8, res.y.abs() as u8, res.z.abs() as u8, 255]
  }

  fn color_to_vec3(&self, c: &[u8; 4]) -> Vector3<f32> {
    Vector3::new(c[0] as f32, c[1] as f32, c[2] as f32)
  }
}
