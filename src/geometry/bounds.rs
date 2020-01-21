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
use super::{
    points::*,
    vectors::*
};
use std::f64;

//TODO implement iteration over bounds maximum extents 
//http://www.pbr-book.org/3ed-2018/Geometry_and_Transformations/Bounding_Boxes.html#Bounds3::Inside

/// Represents a 2 dimensional bounding box, represented by minimum and maximum extents.
/// This is the 2D equivalent of ``Bounds3`` and is somewhat less useful as a construct, but
/// it's useful to have a 2D equivalent to reduce the cost of some computations.
/// 
/// You can construct bounds directly if you want, or you can use the ``Bounds2::new`` function, 
/// which will make copies of the given bounding points.
/// ### Usage
/// ```
/// # use hrend::geometry::bounds::*;
/// # use hrend::geometry::points::*;
/// let min = Point2::new(1., 1.);
/// let max = Point2::new(2., 2.);
/// 
/// let bounds_copy = Bounds2::new(&min, &max);
/// let bounds_owned = Bounds2 {
///     p_min: min,
///     p_max: max
/// };
/// ```
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct Bounds2 {
    pub p_min: Point2,
    pub p_max: Point2
}

impl Bounds2 {

    /// Constructs a new Bounds from references to points. This will copy the points to create
    /// the new Bounds, so it's advisable to only use this if you actually want this behavior.
    /// Otherwise, it's more advisable to create the Bounds2 object directly.
    pub fn new(p_min: &Point2, p_max: &Point2) -> Bounds2 {
        Bounds2 {
            p_min: *p_min,
            p_max: *p_max
        }
    }

    pub fn default() -> Bounds2 {
        let min = f64::MIN;
        let max = f64::MAX;
        Bounds2 {
            p_min: Point2::new(max, max),
            p_max: Point2::new(min, min)
        }
    }

    pub fn from_point(p: &Point2) -> Bounds2 {
        Bounds2 {
            p_min: *p,
            p_max: *p
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
        if d.x > d.y {
            0
        } else {
            1
        }
    }

    pub fn by_ind(&self, ind: usize) -> Point2{
        assert!(1 >= ind); //TODO: figure out if this is expensive and unneeded
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
        let mut o = Vector2::new(p.x - self.p_min.x, p.y - self.p_min.y); // because normal subtraction is just too easy
        if self.p_max.x > self.p_min.x {
            o.x /= self.p_max.x - self.p_min.x;
        }
        if self.p_max.y > self.p_min.y {
            o.y /= self.p_max.y - self.p_min.y;
        }
        o
    }

    pub fn inside(&self, p: &Point2) -> bool {
        p.x >= self.p_min.x && p.x <= self.p_max.x &&
        p.y >= self.p_min.y && p.y <= self.p_max.y
    }

    /// Same thing as `inside`, except it doesn't consider things on the upper bound to be inside.
    /// Not that useful.
    pub fn inside_excl(&self, p: &Point2) -> bool {
        p.x >= self.p_min.x && p.x < self.p_max.x &&
        p.y >= self.p_min.y && p.y < self.p_max.y
    }

    /*
    pub fn mut_bound_sphere(&self, c: &mut Point2, rad: &mut f64) {
        c = &mut (self.p_min + self.p_max).div(2.);
        rad = if inside(c, &mut self) {
            distance(c, self.p_max)
        } else { 0 }
    }
    */

    pub fn bound_sphere(&self, c: &mut Point2, rad: &mut f64) -> (Point2, f64) {
        //TODO: figure out what's actually going on with this function, and what rad does
        (
            (self.p_min + self.p_max).div(2.), 
            if self.inside(c) {c.dist(&self.p_max)} else {0.}
        )
    }

    pub fn overlaps(&self, other: &Bounds2) -> bool {
        self.p_max.x >= other.p_min.x && self.p_min.x <= other.p_max.x &&
        self.p_max.y >= other.p_min.y && self.p_min.y <= other.p_min.y
    }
}

#[derive(Copy, Clone, PartialEq, Debug)]
pub struct Bounds3 {
    pub p_min: Point3,
    pub p_max: Point3
}

impl Bounds3 {
    pub fn new(p1: &Point3, p2: &Point3) -> Bounds3 {
        Bounds3 {
            p_min: *p1,
            p_max: *p2
        }
    }

    pub fn default() -> Bounds3 {
        let min = f64::MIN;
        let max = f64::MAX;
        Bounds3 {
            p_min: Point3::new(max, max, max),
            p_max: Point3::new(min, min, min)
        }
    }

    pub fn from_point(p: &Point3) -> Bounds3 {
        Bounds3 {
            p_min: *p,
            p_max: *p
        }
    }

    pub fn diagonal(&self) -> Vector3 {
        self.p_max - self.p_min
    }

    pub fn surface_area(&self) -> f64 {
        let d = self.diagonal();
        2. * (d.x * d.y + d.x * d.z + d.y * d.z)
    }

    pub fn volume(&self) -> f64 {
        let d = self.diagonal();
        d.x * d.y * d.z
    }

    pub fn corner(&self, corner: usize) -> Point3 {
        Point3 {
            x: self.by_ind(corner & 1).x,
            y: self.by_ind(if (corner & 2) != 0 {1} else {0}).y,
            z: self.by_ind(if (corner & 4) != 0 {1} else {0}).z
        }
    }

    pub fn max_extent(&self) -> isize {
        let d = self.diagonal();
        if d.x > d.y && d.x > d.z{
            0
        } else if d.y > d.z {
            1
        } else {
            2 
        }
    }

    pub fn by_ind(&self, ind: usize) -> Point3 {
        assert!(1 >= ind); // index must be one or two
        if ind == 0 {
            self.p_min
        } else {
            self.p_max
        }
    }

    pub fn lerp(&self, t: &Point3) -> Point3 {
        Point3::new(
            super::lerp(t.x, self.p_min.x, self.p_max.x),
            super::lerp(t.y, self.p_min.y, self.p_max.y),
            super::lerp(t.z, self.p_min.z, self.p_max.z)
        )
    }

    pub fn offset(&self, p: &Point3) -> Vector3 {
        let mut o = Vector3 {
            x: p.x - self.p_min.x,
            y: p.y - self.p_min.y,
            z: p.z - self.p_min.z
        };
        if self.p_max.x > self.p_min.x {
            o.x /= self.p_max.x - self.p_min.x;
        }
        if self.p_max.y > self.p_min.y {
            o.y /= self.p_max.y - self.p_min.y;
        }
        o
    }

    pub fn inside(&self, p: &Point3) -> bool {
        p.x >= self.p_min.x && p.x <= self.p_max.x &&
        p.y >= self.p_min.y && p.y <= self.p_max.y &&
        p.z >= self.p_min.z && p.z <= self.p_min.z
    }

    pub fn overlaps(&self, other: &Bounds3) -> bool {
        self.p_max.x >= other.p_min.x && self.p_min.x <= other.p_min.x &&
        self.p_max.y >= other.p_min.y && self.p_min.y <= other.p_min.y &&
        self.p_max.z >= other.p_min.z && self.p_min.z <= other.p_min.z
    }

    pub fn expand(&self, by: f64) -> Bounds3 {
        Bounds3 {
            p_min: self.p_min.sub_vec(&Vector3::new(by, by, by)),
            p_max: self.p_max.add_vec(&Vector3::new(by, by, by))
        }
    }

    pub fn expand_mut(&mut self, by: f64) {

    }

    /*
    pub fn mut_bound_sphere(&self, c: &mut Point2, rad: &mut f64) {
        c = (self.p_min + self.p_max).div(2.);
        rad = if inside(c, &mut self) {
            distance(c, self.p_max)
        } else { 0 }
    }
    */

    pub fn bound_sphere(&self, c: &Point3) -> (Point3, f64) {
        (
            (self.p_min + self.p_max).div(2.),
            if self.inside(c) {self.p_max.dist(c)} else {0.}
        )
    }
    
}