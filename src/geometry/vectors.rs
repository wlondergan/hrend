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
use std::ops;
use std::fmt;
// TODO: add iteration to vector types

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
        let recip = 1./by;
        Vector2::new(self.x * recip, self.y * recip)
    }

    pub fn div_mut(&mut self, by: f64) {
        let recip = 1./by;
        self.x *= recip;
        self.y *= recip; // floating point division is expensive.
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

impl ops::Neg for Vector2 {
    type Output = Vector2;
    fn neg(self) -> Vector2 {
        Vector2::new(-self.x, -self.y)
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
        let recip = 1./by;
        Vector3::new(self.x * recip, self.y * recip, self.z * recip)
    }

    pub fn div_mut(&mut self, by: f64) {
        let recip = 1./by;
        self.x *= recip;
        self.y *= recip;
        self.z *= recip;
    }

    //FIXME: this probably doesn't work.
    pub fn min_component(&self) -> f64 {
        *[self.x, self.y, self.z].iter().enumerate().min_by(|x,y| x.partial_cmp(y).unwrap()).unwrap().1
    }

    pub fn max_component(&self) -> f64 {
        *[self.x, self.y, self.z].iter().enumerate().max_by(|x,y| x.partial_cmp(y).unwrap()).unwrap().1
    }

    pub fn min_ind(&self) -> usize {
        [self.x, self.y, self.z].iter().enumerate().min_by(|x,y| x.partial_cmp(y).unwrap()).unwrap().0
    }

    pub fn max_ind(&self) -> usize {
        [self.x, self.y, self.z].iter().enumerate().max_by(|x,y| x.partial_cmp(y).unwrap()).unwrap().0
    }

    pub fn by_ind(&self, ind: usize) -> f64 {
        match ind {
            0 => self.x,
            1 => self.y,
            2 => self.z,
            _ => panic!("Index out of range")
        }
    }

    pub fn permute(&self, x: usize, y: usize, z: usize) -> Vector3 {
        Vector3 {
            x: self.by_ind(x),
            y: self.by_ind(y),
            z: self.by_ind(z)
        }
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

impl ops::Neg for Vector3 {
    type Output = Vector3;
    fn neg(self) -> Vector3 {
        Vector3::new(-self.x, -self.y, -self.z)
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

pub fn dot_2(a: &Vector2, b: &Vector2) -> f64 {
    a.x * b.x + a.y * b.y
}

pub fn abs_dot_2(a: &Vector2, b: &Vector2) -> f64 {
    dot_2(a, b).abs()
}

/// Computes the dot product of two 3-vectors.
pub fn dot(a: &Vector3, b: &Vector3) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}

/// Takes the absolute value of the dot product of two vectors.
pub fn abs_dot(a: &Vector3, b: &Vector3) -> f64 {
    dot(a, b).abs()
}

/// Computes the cross product of two vectors.
pub fn cross(a: &Vector3, b: &Vector3) -> Vector3 {
    Vector3 {
        x: a.y * b.z - a.z * b.y,
        y: a.z * b.x - a.x * b.z,
        z: a.x * b.y - a.y * b.x
    }
}

pub fn norm_2(a: &Vector2) -> Vector2 {
    a.div(a.length())
}

pub fn norm(a: &Vector3) -> Vector3 {
    a.div(a.length())
}

///# Safety
/// 
/// There aren't really any issues with safety here, as long as you're okay with having these pointers rewritten afterwards.
pub unsafe fn coord_sys(a: &Vector3, b: *mut Vector3, c: *mut Vector3) {
    if a.x.abs() > a.y.abs() {
        *b = Vector3::new(-a.z, 0., a.x).div(a.length())
    }
    else {
        *b = Vector3::new(0., a.z, -a.y).div(a.length())
    }
    *c = cross(a, &*b);
}

/// A safer version of the previous example, except it allocates new vectors instead.
pub fn safe_coord_sys(a: &Vector3) -> (Vector3, Vector3) {
    let b = if a.x.abs() > a.y.abs() {
        Vector3::new(-a.z, 0., a.x).div(a.length())
    }
    else {
        Vector3::new(0., a.z, -a.y).div(a.length()) 
    };
    (b, cross(a, &b))
}


