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

/// Provides implementations for both 2D and 3D vectors.
pub mod vectors;

/// Provides implementations for 2D and 3D points, similarly to vectors.
pub mod points;

/// Provides implementation for surface normals.
pub mod normal;

/// Provides implementation for rays.
pub mod ray;

/// Provides implementation for square bounds in 2 and 3D.
pub mod bounds;

/// Provides implementation for linear transformations on homogenous coordinates.
pub mod transform;

/// Linearly interpolates on a line given a point `t` along the line to find and points 
///`v1` and `v2` defining the line segment.
pub fn lerp(t: f64, v1: f64, v2: f64) -> f64 {
    (1. - t) * v1 + t * v2
}
