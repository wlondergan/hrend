/* hrend - A Toy Rendering Engine
*
* Author: Hughes Londergan <whlondergan@gmail.com>
* 
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
use std::f64;
use std::fmt;
use super::vectors::*;

/// Represents a point in 3D (Euclidean) space, with double precision floating points for each
/// component.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct Point3 {
    pub x: f64,
    pub y: f64,
    pub z: f64
}

impl Point3 {
    pub fn new(x: f64, y: f64, z: f64) -> Point3 {
        Point3 {
            x,
            y,
            z
        }
    }

    /// Determines whether any of the fields of this struct are `NaN`.
    pub fn has_nan(&self) -> bool {
        f64::is_nan(self.x) || f64::is_nan(self.y) || f64::is_nan(self.z)
    }

    pub fn add_vec(&self, other: &Vector3) -> Point3 {
        Point3::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }

    pub fn add_vec_mut(&mut self, other: &Vector3) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }

    pub fn sub_vec(&self, other: &Vector3) -> Point3 {
        Point3::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }

    pub fn sub_vec_mut(&mut self, other: &Vector3) {
        self.x -= other.x;
        self.y -= other.y;
        self.z -= other.z;
    }

    pub fn dist(&self, b: &Point3) -> f64 {
        (*self - *b).length()
    }

    pub fn dist_sqr(&self, b: &Point3) -> f64 {
        (*self - *b).length_sqr()
    }

    pub fn mult(&self, other: f64) -> Point3 {
        Point3::new(self.x * other, self.y * other, self.z * other)
    }

    pub fn mult_mut(&mut self, other: f64) {
        self.x *= other;
        self.y *= other;
        self.z *= other;
    }

    pub fn floor(&self) -> Point3 {
        Point3 {
            x: self.x.floor(),
            y: self.y.floor(),
            z: self.z.floor()
        }
    }

    pub fn ceil(&self) -> Point3 {
        Point3 {
            x: self.x.ceil(),
            y: self.y.ceil(),
            z: self.z.ceil()
        }
    }

    pub fn abs(&self) -> Point3 {
        Point3 {
            x: self.x.abs(),
            y: self.y.abs(),
            z: self.z.abs()
        }
    }

    pub fn by_ind(&self, ind: usize) -> f64 {
        match ind {
            0 => self.x,
            1 => self.y,
            2 => self.z,
            _ => panic!("Index out of range")
        }
    }

    pub fn permute(&self, x: usize, y: usize, z: usize) -> Point3 {
        Point3 {
            x: self.by_ind(x),
            y: self.by_ind(y),
            z: self.by_ind(z)
        }
    }
}

impl std::ops::Sub for Point3 {
    type Output = Vector3;

    fn sub(self, other: Point3) -> Self::Output {
        Vector3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z
        }
    }
}

impl fmt::Display for Point3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}, {}, {}]", self.x, self.y, self.z)
    }
}

impl std::ops::Add for Point3 {
    type Output = Point3;
    fn add(self, other: Point3) -> Point3 {
        Point3 {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z
        }
    }
}

/// Represents a point in 2D (Cartesian) space, with double precision floating points for each
/// component.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct Point2 {
    pub x: f64,
    pub y: f64
}

impl Point2 {
    pub fn new(x: f64, y: f64) -> Point2 {
        Point2 {
            x,
            y
        }
    }

    /// Determines whether any of the fields of this struct are `NaN`.
    pub fn has_nan(&self) -> bool {
        f64::is_nan(self.x) || f64::is_nan(self.y)
    }
}

/// Linearly interpolates between the points.
pub fn lerp(t: f64, a: &Point3, b: &Point3) -> Point3 {
    a.mult(1. - t) + b.mult(t)
}

/// Gives the minimum of each component of the two points.
pub fn min(a: &Vector3, b: &Vector3) -> Vector3 {
    Vector3 {
        x: a.x.min(b.x),
        y: a.y.min(b.y),
        z: a.z.min(b.z)
    }
}

pub fn max(a: &Vector3, b: &Vector3) -> Vector3 {
    Vector3 {
        x: a.x.max(b.x),
        y: a.y.max(b.y),
        z: a.z.max(b.z)
    }
}

//TODO: implement all of these functions for Point2Ds as well.