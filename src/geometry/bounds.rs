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
use super::points::*;
use super::vectors::*;
use std::f64;

#[derive(Copy, Clone, PartialEq, Debug)]
pub struct Bounds2 {
    pub p_min: Point2,
    pub p_max: Point2
}

impl Bounds2 {
    pub fn new(p1: &Point2, p2: &Point2) -> Bounds2 {
        Bounds2 {
            p_min: p1,
            p_max: p2
        }
    }

    pub fn default() -> Bounds2 {
        let min = f64::MIN;
        let max = f64::MAX;
        Bounds2 {
            p_min: Point2::new(max, max),
            p_max: Point2::new(min, min);
        }
    }

    pub fn from_point(p: &Point2) {
        Bounds2 {
            p_min: p,
            p_max: p
        }
    }

    pub fn diagonal(&self) -> Vector2 {
        self.p_max - self.p_min
    }

    pub fn area(&self) -> f64 {
        let d = self.diagonal();
        d.x * d.y
    }

    pub fn max_extent(&self) -> isize {
        let d = self.diagonal();
        if diag.x > diag.y {
            0
        } else {
            1
        }
    }

    pub fn by_ind(&self, ind: usize) {
        assert!(0 <= ind && 1 >= ind);
        if ind == 0 {
            self.p_min
        } else {
            self.p_max
        }
    }

    pub fn lerp(&self, t: &Point2) -> Point2 {
        Point2::new(
            super::lerp(t.x, self.p_min.x, self.p_max.x),
            super::lerp(t.y, self.p_min.y, self.p_max.y)
        )
    }

    pub fn offset(&self, p: &Point2) -> Vector2 {
        let mut o = p - self.p_min;
        if self.p_max.x > self.p_min.x {
            o.x /= self.p_max.x - self.p_min.x;
        }
        if self.p_max.y > self.p_min.y {
            o.y /= self.p_max.y - self.p_min.y;
        }
        o
    }

    pub fn mut_bound_sphere(&self, c: &mut Point2, rad: &mut f64) {
        c = (self.p_min + self.p_max).div(2.);
        rad = if inside(c, &mut self) {
            distance(c, self.p_max)
        } else { 0 }
    }

    /* TODO implement
    pub fn bound_sphere(&self) -> (Point2, f64) {

    }
    */
}

pub struct Bounds3 {
    pub pMin: Point3,
    pub pMax: Point3
}

impl Bounds3 {

}