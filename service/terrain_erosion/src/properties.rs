#[derive(Debug, Clone)]
pub struct DirectionalProperties {
    pub left: f64,
    pub right: f64,
    pub top: f64,
    pub bottom: f64
}

#[derive(Debug, Clone)]
pub struct CellProperties {
  pub terrain_height: f64,
  pub water_height: f64,
  pub suspended_sediment: f64,
  pub water_outflow_flux: DirectionalProperties,
  pub velocity: (f64, f64),
}

pub struct CellPropertiesBuilder {
  terrain_height: f64,
  water_height: f64,
  suspended_sediment: f64,
  water_outflow_flux: DirectionalProperties,
  velocity: (f64, f64),
}

impl CellPropertiesBuilder {
  pub fn new() -> Self {
      CellPropertiesBuilder {
          terrain_height: 0.0,
          water_height: 0.0,
          suspended_sediment: 0.0,
          water_outflow_flux: DirectionalProperties {
            left: f64::NAN,
            right: f64::NAN,
            top: f64::NAN,
            bottom: f64::NAN
          },
          velocity: (0.0, 0.0),
      }
  }

  pub fn terrain_height(mut self, value: f64) -> Self {
      self.terrain_height = value;
      self
  }

  pub fn water_height(mut self, value: f64) -> Self {
      self.water_height = value;
      self
  }

  pub fn suspended_sediment(mut self, value: f64) -> Self {
      self.suspended_sediment = value;
      self
  }

  pub fn water_outflow_flux(mut self, value: DirectionalProperties) -> Self {
      self.water_outflow_flux = value;
      self
  }

  pub fn velocity(mut self, value: (f64, f64)) -> Self {
      self.velocity = value;
      self
  }

  pub fn build(self) -> CellProperties {
      CellProperties {
          terrain_height: self.terrain_height,
          water_height: self.water_height,
          suspended_sediment: self.suspended_sediment,
          water_outflow_flux: self.water_outflow_flux,
          velocity: self.velocity,
      }
  }
}