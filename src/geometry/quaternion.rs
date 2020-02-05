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

use super::vectors::*;
use std::ops::{
    Add, Sub, SubAssign, AddAssign, Mul
};
use crate::geometry::transform::*;
use crate::math::*;

/// Represents a "Quaternion", or the 4D imaginary number. Here, we use it to
/// compute pure rotations.
#[derive(PartialEq, Copy, Clone)] 
pub struct Quaternion {
    pub v: Vector3, //the imaginary "vector" part of the quaternion
    pub w: f64 //the real part of the quaternion
}

impl Quaternion {

    /// Constructs a unit quaternion.
    pub fn default() -> Quaternion {
        Quaternion {
            v: Vector3::new(0., 0., 0.),
            w: 1.
        }
    }

    pub fn new(v: &Vector3, w: f64) -> Quaternion {
        Quaternion {
            v: *v,
            w
        } 
    }

    pub fn dot(&self, other: &Quaternion) -> f64 {
        dot(&self.v, &other.v) + self.w * other.w
    }

    /// Creates a Quaternion from a Transform, via the method suggested in
    /// `pbrt`. According to the `pbrt` documentation, we should consult
    /// Shoemake (1991).
    #[allow(clippy::many_single_char_names)]
    pub fn from_transform(t: &Transform) -> Quaternion {
        let trace = t.m.m[0][0] + t.m.m[1][1] + t.m.m[2][2];

        if trace > 0. {
            let mut s = f64::sqrt(trace + 1.);
            let w = s / 2.;
            s = 0.5 / s;
            Quaternion::new(
                &Vector3::new(
                    (t.m.m[2][1] - t.m.m[1][2]) * s,
                    (t.m.m[0][2] - t.m.m[2][0]) * s,
                    (t.m.m[1][0] - t.m.m[0][1]) * s),
                w,
            )
        } else {
            let nxt = [1, 2, 0];
            let mut q = [0., 0., 0.];
            let mut i = if t.m.m[1][1] > t.m.m[0][0] {
                1
            } else {
                0
            };
            if t.m.m[2][2] > t.m.m[i][i] {
                i = 2;
            }
            let j = nxt[i];
            let k = nxt[j];
            let mut s = f64::sqrt((t.m.m[i][i] - (t.m.m[j][j] + t.m.m[k][k])) + 1.);
            q[i] = s * 0.5;
            if s != 0. {
                s = 0.5/s;
            }
            q[j] = (t.m.m[j][i] + t.m.m[i][j]) * s;
            q[k] = (t.m.m[k][i] + t.m.m[i][k]) * s;
            Quaternion {
                w: (t.m.m[k][j] - t.m.m[j][k]) * s,
                v: Vector3 {
                    x: q[0],
                    y: q[1],
                    z: q[2]
                }
            } //can probably be optimized.
        }
    }

    pub fn norm(&self) -> Quaternion {
        self.div(f64::sqrt(self.dot(&self)))
    }

    pub fn mult(&self, by: f64) -> Quaternion {
        Quaternion::new(&self.v.mult(by), self.w * by)
    }

    pub fn mult_mut(&mut self, by: f64) {
        self.v.mult_mut(by);
        self.w *= by;
    }

    pub fn div(&self, by: f64) -> Quaternion {
        let recip = 1./by;
        self.mult(recip)
    }

    pub fn div_mut(&mut self, by: f64) {
        let recip = 1./by;
        self.mult_mut(recip);
    }

    pub fn to_transform(&self) -> Transform {
        // per "quaternion.cpp" from pbrt
        let xx = self.v.x * self.v.x;
        let yy = self.v.y * self.v.y;
        let zz = self.v.z * self.v.z;
        let xy = self.v.x * self.v.y;
        let xz = self.v.x * self.v.z;
        let yz = self.v.y * self.v.z;
        let wx = self.v.x * self.w;
        let wy = self.v.y * self.w;
        let wz = self.v.z * self.w;

        let mut m = IDENTITY; //the identity matrix

        m.m[0][0] = 1. - 2. * (yy + zz);
        m.m[0][1] = 2. * (xy + wz);
        m.m[0][2] = 2. * (xz - wy);
        m.m[1][0] = 2. * (xy - wz);
        m.m[1][1] = 1. - 2. * (xx + zz);
        m.m[1][2] = 2. * (yz + wx);
        m.m[2][0] = 2. * (xz + wy);
        m.m[2][1] = 2. * (yz - wx);
        m.m[2][2] = 1. - 2. * (xx + yy);

        Transform::from_mat(m.transpose().m) //TODO little gross to unwrap the matrix like this
    }

    /// Performs spherical linear interpolation between q1 and q2.
    pub fn slerp(t: f64, q1: &Quaternion, q2: &Quaternion) -> Quaternion {
        let cos = q1.dot(q2);
        if cos > 0.9995 { // if we're close enough simply use linear interpolation
            q1.mult(1. - t) + q2.mult(t)
        } else { // otherwise we have to use the expensive trig calculations
            let theta = f64::acos(clamp(cos, -1., 1.));
            let thetap = theta * t;

            let q_perp = (*q2 - *q1).mult(cos).norm();
            q1.mult(f64::cos(thetap)) + q_perp.mult(f64::sin(thetap))
        }
    }

}

//dot, normalize, to transform, from transform, slerp

impl Add for Quaternion {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Quaternion::new(&(self.v + other.v), self.w + other.w)
    }
}

impl AddAssign for Quaternion {
    fn add_assign(&mut self, other: Self) {
        self.v += other.v;
        self.w += other.w;
    }
}

impl Sub for Quaternion {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Quaternion::new(&(self.v - other.v), self.w - other.w)
    }
}

impl SubAssign for Quaternion {
    fn sub_assign(&mut self, other: Self) {
        self.v -= other.v;
        self.w -= other.w;
    }
}

impl Mul for Quaternion {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Quaternion {
            v: cross(&self.v, &other.v) + other.v.mult(self.w) + self.v.mult(other.w),
            w: self.w * other.w - dot(&self.v, &other.v)
        }
    }
}

/* find a better way to do this.
impl MulAssign for Quaternion {
    fn mul_assign(&mut self, other: Self) {
        self = &(*self * other);
    }
}
*/


