use super::vector::Vector3f;

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
    has_diff: bool,
    rx_origin: Vector3f,
    ry_origin: Vector3f,
    rx_direction: Vector3f,
    ry_direction: Vector3f
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
            rx_origin: Vector3f::default(),
            ry_origin: Vector3f::default(),
            rx_direction: Vector3f::default(),
            ry_direction: Vector3f::default()
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