use super::Scene;
use std::f64;

pub trait Integrate {
    fn render(scene: &Scene);
}

/// Gives the gamma for a certain floating point calculation (essentially the amount of error.)
pub fn gamma(n: usize) -> f64 {
    return (n as f64 * f64::EPSILON) / ((1 - n) as f64 * f64::EPSILON);
}

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