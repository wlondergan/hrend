use std::f64;
use std::ops;
use std::fmt;

/// Represents a vector in 2 dimensional space.
/// The fields of this struct are ``f64``, which is unlike `pbrt`'s approach. They have a runtime flag that lets
/// the user choose between `f64` and `f32`, which tends to complicate things unnecessarily for a toy renderer
/// like this one.
/// I have chosen to exclude the indexing operator included in `pbrt`. It doesn't seem important.
#[derive(Copy, Clone, PartialEq)]
pub struct Vector2 {
    pub x: f64,
    pub y: f64
}

impl Vector2 {
    /// makes a `Vector2`. Equivalent to just building it manually, but in more Java-y style.
    pub fn new(x: f64, y: f64) -> Vector2 {
        Vector2 {
            x,
            y
        }
    }

    /// Determines whether any of the fields of this struct are `NaN`.
    pub fn has_nan(&self) -> bool {
        f64::is_nan(self.x) || f64::is_nan(self.y)
    }

    // Returns the length of this vector, squared.
    pub fn length_sqr(&self) -> f64 {
        self.x*self.x + self.y*self.y
    }

    // Returns the length of this vector. Expensive, because it uses a `sqrt` call, which isn't ideal.
    pub fn length(&self) -> f64 {
        f64::sqrt(self.length_sqr())
    }

    pub fn mult(&self, by: f64) -> Vector2 {
        Vector2::new(self.x * by, self.y * by)
    }

    pub fn mult_mut(&mut self, by: f64) {
        self.x *= by;
        self.y *= by
    }

    pub fn div(&self, by: f64) -> Vector2 {
        Vector2::new(self.x / by, self.y / by)
    }

    pub fn div_mut(&mut self, by: f64) {
        self.x /= by;
        self.y /= by;
    }
}

impl ops::Add for Vector2 {
    type Output = Vector2;

    fn add(self, other: Vector2) -> Vector2 {
        assert!(!other.has_nan());

        Vector2 {
            x: self.x + other.x,
            y: self.y + other.y
        }
    }
}

impl ops::Sub for Vector2 {
    type Output = Vector2;

    fn sub(self, other: Vector2) -> Vector2 {
        assert!(!other.has_nan());

        Vector2 {
            x: self.x - other.x,
            y: self.y - other.y
        }
    }
}

impl ops::AddAssign for Vector2 {
    fn add_assign(&mut self, other: Self) {
        assert!(!other.has_nan());
        self.x += other.x;
        self.y += other.y;
    }
}

impl ops::SubAssign for Vector2 {
    fn sub_assign(&mut self, other: Self) {
        assert!(!other.has_nan());
        self.x -= other.x;
        self.y -= other.y;
    }
}

impl fmt::Display for Vector2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}, {}]", self.x, self.y)
    }
}

/// The zero vector in 2-dimensional space.
pub const ZERO_2: Vector2 = Vector2 {
    x: 0.,
    y: 0.
};

/// Represents a vector in 3 dimensional space. See notes on `Vector2` for more information about the
/// specific implementation.
#[derive(Copy, Clone, PartialEq)]
pub struct Vector3 {
    pub x: f64,
    pub y: f64,
    pub z: f64
}

impl Vector3 {
    /// makes a `Vector3`. Equivalent to just building it manually, but in more Java-y style.
    pub fn new(x: f64, y: f64, z: f64) -> Vector3 {
        Vector3 {
            x,
            y,
            z
        }
    }

    /// Determines whether any of the fields of this struct are `NaN`.
    pub fn has_nan(&self) -> bool {
        f64::is_nan(self.x) || f64::is_nan(self.y) || f64::is_nan(self.z)
    }

    // Returns the length of this vector, squared.
    pub fn length_sqr(&self) -> f64 {
        self.x*self.x + self.y*self.y + self.z * self.z
    }

    // Returns the length of this vector. Expensive, because it uses a `sqrt` call, which isn't ideal.
    pub fn length(&self) -> f64 {
        f64::sqrt(self.length_sqr())
    }

    pub fn mult(&self, by: f64) -> Vector3 {
        Vector3::new(self.x * by, self.y * by, self.z * by)
    }

    pub fn mult_mut(&mut self, by: f64) {
        self.x *= by;
        self.y *= by;
        self.z *= by;
    }

    pub fn div(&self, by: f64) -> Vector3 {
        Vector3::new(self.x / by, self.y / by, self.z / by)
    }

    pub fn div_mut(&mut self, by: f64) {
        self.x /= by;
        self.y /= by;
        self.z /= by;
    }
}

impl ops::Add for Vector3 {
    type Output = Vector3;

    fn add(self, other: Vector3) -> Vector3 {
        assert!(!other.has_nan());

        Vector3 {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z
        }
    }
}

impl ops::Sub for Vector3 {
    type Output = Vector3;

    fn sub(self, other: Vector3) -> Vector3 {
        assert!(!other.has_nan());

        Vector3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z
        }
    }
}

impl ops::AddAssign for Vector3 {
    fn add_assign(&mut self, other: Self) {
        assert!(!other.has_nan());
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl ops::SubAssign for Vector3 {
    fn sub_assign(&mut self, other: Self) {
        assert!(!other.has_nan());
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
    }
}

impl fmt::Display for Vector3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}, {}, {}]", self.x, self.y, self.z)
    }
}

pub const ZERO_3: Vector3 = Vector3 {
    x: 0.,
    y: 0.,
    z: 0.
};