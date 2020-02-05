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
use super::Scene;
use std::f64;

/// There are various implementations of integrators (Whitted, &c). This trait allows implementation
/// of different complexities of integration, so that the actual renderer calculations can be replaced 
/// with modular code easily.
pub trait Integrate {
    fn render(scene: &Scene);
}

/// Gives the gamma for a certain floating point calculation (essentially the amount of error.)
pub fn gamma(n: usize) -> f64 {
    (n as f64 * f64::EPSILON) / ((1 - n) as f64 * f64::EPSILON)
}

/// If v compares less than lo, returns lo; otherwise if hi compares less than v,
/// returns hi; otherwise returns v. TODO make generic
pub fn clamp(v: f64, lo: f64, hi: f64) -> f64 {
    if v < lo {
        lo
    } else if v < hi {
        hi
    } else {
        v
    }
}

/// returns the maximum of the two elements.
pub fn max<T: PartialOrd>(a: T, b: T) -> T {
    if a > b { a } else { b }
}

/// returns the minimum of two elements.
pub fn min<T: PartialOrd>(a: T, b: T) -> T {
    if a < b { a } else { b }
}

pub fn lerp(t: f64, v1: f64, v2: f64) -> f64 {
    (1. - t) * v1 + t * v2
}

/// Determines whether or not two floats are equal (within the floating point epsilon).
pub fn eq_f64(a: f64, b: f64) -> bool {
    (a - b).abs() <= f64::EPSILON
}