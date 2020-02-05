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

#[derive(Copy, Clone, PartialEq, Debug)]
pub struct Normal3 {
    pub x: f64,
    pub y: f64,
    pub z: f64
}

impl Normal3 {
    pub fn new(x: f64, y: f64, z: f64) -> Normal3 {
        Normal3 {x, y, z}
    }

    pub fn from_vec3(from: Vector3) -> Normal3 {
        Normal3::new(from.x, from.y, from.z)
    }

    pub fn default() -> Normal3 {
        Normal3::new(0., 0., 0.)
    }

    pub fn mult(&self, by: f64) -> Normal3 {
        Normal3::new(self.x * by, self.y * by, self.z * by)
    }

    pub fn mult_mut(&mut self, by: f64) {
        self.x *= by;
        self.y *= by;
        self.z *= by;
    }

    // TODO make this more efficient
    /// Makes a vector face in the same hemisphere as the given vector.
    pub fn faceforward(&self, other: &Vector3) -> Normal3 {
        if dot(&Vector3::from_norm3(&self), other) < 0. {
            -*self
        }
        else {
            *self
        }
    }
}

impl std::ops::Neg for Normal3 {
    type Output = Normal3;
    fn neg(self) -> Normal3 {
        Normal3::new(-self.x, -self.y, -self.z)
    }
}
//TODO implement other operators for this class