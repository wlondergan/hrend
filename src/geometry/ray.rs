use super::vector::{Vector, Vector3f};

pub trait GenericRay {
    fn o(&self) -> Vector3f;
    fn d(&self) -> Vector3f;
    fn time(&self) -> f32;

    fn at(&self, t: f32) -> Vector3f {
        self.o() + self.d() * t
    }
}

#[derive(Clone, Copy)]
pub struct Ray {
    pub o: Vector3f,
    pub d: Vector3f,
    pub time: f32,
    // medium: Medium
}

pub struct DiffRay {
    ray: Ray,
    pub has_diff: bool,
    pub rx_origin: Vector3f,
    pub ry_origin: Vector3f,
    pub rx_direction: Vector3f,
    pub ry_direction: Vector3f
}

impl Ray {
    pub fn new(o: Vector3f, d: Vector3f) -> Ray {
        Ray {o, d, time: 0.0}
    }
}

impl GenericRay for Ray {
    fn o(&self) -> Vector3f {
        self.o
    }

    fn d(&self) -> Vector3f {
        self.d
    }

    fn time(&self) -> f32 {
        self.time
    }
}

impl DiffRay {
    
    pub fn new(ray: &Ray) -> DiffRay {
        DiffRay {
            ray: *ray,
            has_diff: false,
            rx_origin: Vector3f::all_x(0.0),
            ry_origin: Vector3f::all_x(0.0),
            rx_direction: Vector3f::all_x(0.0),
            ry_direction: Vector3f::all_x(0.0)
        }
    }

    pub fn scale_diffs(&mut self, s: f32) {
        self.rx_origin = self.o() + (self.rx_origin - self.o()) * s;
        self.ry_origin = self.o() + (self.ry_origin - self.o()) * s;
        self.rx_direction = self.d() + (self.rx_direction - self.d()) * s;
        self.ry_direction = self.d() + (self.ry_direction - self.d()) * s;
    }

}

impl GenericRay for DiffRay {
    fn o(&self) -> Vector3f {
        self.ray.o()
    }

    fn d(&self) -> Vector3f {
        self.ray.d()
    }

    fn time(&self) -> f32 {
        self.ray.time()
    }
}