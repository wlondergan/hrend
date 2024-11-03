use std::ops::Index;
use std::marker::PhantomData;

use crate::math::{self, square, Num};
use super::vector::{Vector, Vector2, Vector3, Vector2i, Vector2f, Vector3i, Vector3f};

// TODO implement an iterator (or other algorithm somewhere) which allows iteration over integer bounds within a given bounding box

struct Bounds<T: Num, V: Vector<T>> {
    pub min_point: V,
    pub max_point: V,
    _t: PhantomData<T>
}

impl<T: Num, V: Vector<T>> Bounds<T, V> {

    fn default() -> Self {
        Self {
            min_point: V::default(),
            max_point: V::maximum(),
            _t: PhantomData
        }
    }

    fn new(p1: V, p2: V) -> Self {
        Self {
            min_point: Vector::min(p1, p2),
            max_point: Vector::max(p1, p2),
            _t: PhantomData
        }
    }

    fn from_point(p1: V) -> Self {
        Self {
            min_point: p1,
            max_point: p1,
            _t: PhantomData
        }
    }

    /// Adds the given point `p` to the existing bounds `b` to yield a new bounding box.
    fn add_point(b: &Self, p: V) -> Self {
        Bounds {
            min_point: Vector::min(b.min_point, p),
            max_point: Vector::max(b.max_point, p),
            _t: PhantomData
        }
    }

    fn union(b1: &Self, b2: &Self) -> Self {
        Bounds {
            min_point: Vector::min(b1.min_point, b2.min_point),
            max_point: Vector::max(b1.max_point, b2.max_point),
            _t: PhantomData
        }
    }

    fn intersection(b1: &Self, b2: &Self) -> Self {
        Bounds {
            min_point: Vector::max(b1.min_point, b2.min_point),
            max_point: Vector::min(b1.max_point, b2.max_point),
            _t: PhantomData
        }
    }

    fn expand(b: &Self, by: T) -> Self {
        Bounds {
            min_point: b.min_point - Vector::all_x(by),
            max_point: b.max_point + Vector::all_x(by),
            _t: PhantomData
        }
    }

    fn diagonal(b: &Self) -> V {
        b.max_point - b.min_point
    }

    fn volume(b: &Self) -> T {
        Vector::h_prod(Self::diagonal(b))
    }

}

impl<T: Num> Bounds<T, Vector2<T>> {

    fn corner(&self, c: usize) -> Vector2<T> {
        Vector2 {
            x: self[c & 1].x,
            y: self[c & 2].y
        }
    }

    fn overlap(b1: &Self, b2: &Self) -> bool {
        b1.max_point.x >= b2.min_point.x && b1.min_point.x <= b2.max_point.x &&
        b1.max_point.y >= b2.min_point.y && b1.min_point.y <= b2.max_point.y
    }

    fn inside(p: Vector2<T>, b: &Self) -> bool {
        p.x >= b.min_point.x && p.x <= b.max_point.x &&
        p.y >= b.min_point.y && p.y <= b.max_point.y
    }

    /// Same as `inside`, but excludes the top bound (generally only useful with integer bounds, but implemented
    /// generically anyways.)
    fn inside_exclusive(p: Vector2<T>, b: &Self) -> bool {
        p.x >= b.min_point.x && p.x < b.max_point.x &&
        p.y >= b.min_point.y && p.y < b.max_point.y
    }

    fn distance_squared(p: Vector2<T>, b: &Self) -> f32 {
        let minp = Vector::as_floats(b.min_point);
        let maxp = Vector::as_floats(b.max_point);
        let pf = Vector::as_floats(p);
        let dist_x = f32::max(0.0, f32::max(minp.x - pf.x, pf.x - maxp.x));
        let dist_y = f32::max(0.0, f32::max(minp.y - pf.y, pf.y - maxp.y));
        square(dist_x) + square(dist_y)
    }

    fn distance(p: Vector2<T>, b: &Self) -> f32 {
        Self::distance_squared(p, b)
    }

    fn empty(&self) -> bool {
        self.min_point.x >= self.max_point.x || self.min_point.y >= self.max_point.y
    }

    fn degenerate(&self) -> bool {
        self.min_point.x >  self.max_point.x || self.min_point.y >= self.max_point.y
    }

}

impl<T: Num> Bounds<T, Vector3<T>> {
    
    fn corner(&self, c: usize) -> Vector3<T> {
        Vector3 {
            x: self[c & 1].x,
            y: self[c & 2].y,
            z: self[c & 4].z
        }
    }

    fn overlap(b1: &Self, b2: &Self) -> bool {
        b1.max_point.x >= b2.min_point.x && b1.min_point.x <= b2.max_point.x &&
        b1.max_point.y >= b2.min_point.y && b1.min_point.y <= b2.max_point.y &&
        b1.max_point.z >= b2.min_point.z && b1.min_point.z <= b2.max_point.z 
    }

    fn inside(p: Vector3<T>, b: &Self) -> bool {
        p.x >= b.min_point.x && p.x <= b.max_point.x &&
        p.y >= b.min_point.y && p.y <= b.max_point.y &&
        p.z >= b.min_point.z && p.z <= b.max_point.z
    }

    fn inside_exclusive(p: Vector3<T>, b: &Self) -> bool {
        p.x >= b.min_point.x && p.x <= b.max_point.x &&
        p.y >= b.min_point.y && p.y <= b.max_point.y &&
        p.z >= b.min_point.z && p.z <= b.max_point.z
    }

    fn distance_squared(p: Vector3<T>, b: &Self) -> f32 {
        let minp = Vector::as_floats(b.min_point);
        let maxp = Vector::as_floats(b.max_point);
        let pf = Vector::as_floats(p);
        let dist_x = f32::max(0.0, f32::max(minp.x - pf.x, pf.x - maxp.x));
        let dist_y = f32::max(0.0, f32::max(minp.y - pf.y, pf.y - maxp.y));
        let dist_z = f32::max(0.0, f32::max(minp.z - pf.z, pf.z - maxp.z));
        square(dist_x) + square(dist_y) + square(dist_z)
    }

    fn distance(p: Vector3<T>, b: &Self) -> f32 {
        Num::sqrt(Self::distance_squared(p, b))
    }

    fn surface_area(b: &Self) -> T {
        let d = Self::diagonal(b);
        T::from_i32(2) * d.x * d.y + d.x * d.z + d.y * d.z
    }

    fn max_dimension(b: &Self) -> usize {
        let d = Self::diagonal(b);
        if d.x > d.y && d.x > d.z {
            0
        } else if d.y > d.z {
            1
        } else {
            2
        }
    }

    fn bounding_sphere(&self) -> (Vector3<T>, f32) {
        let center = (self.max_point + self.min_point) / T::from_i32(2);
        (
            center,
            if Self::inside(center, self) {Vector::distance(center, self.max_point)} else {0.0}
        )
    }

    fn empty(&self) -> bool {
        self.min_point.x >= self.max_point.x || self.min_point.y >= self.max_point.y || self.min_point.z >= self.max_point.z
    }

    fn degenerate(&self) -> bool {
        self.min_point.x >  self.max_point.x || self.min_point.y >= self.max_point.y || self.min_point.z > self.max_point.z
    }

}

impl Bounds<i32, Vector2<i32>> {
    fn distance_squared_int(p: Vector2<i32>, b: &Self) -> i32 {
        let dist_x = Ord::max(0, Ord::max(b.min_point.x - p.x, p.x - b.max_point.x));
        let dist_y = Ord::max(0, Ord::max(b.min_point.y - p.y, p.y - b.max_point.y));
        square(dist_x) + square(dist_y)
    }
}

impl Bounds<i32, Vector3<i32>> {
    fn distance_squared_int(p: Vector3<i32>, b: &Self) -> i32 {
        let dist_x = Ord::max(0, Ord::max(b.min_point.x - p.x, p.x - b.max_point.x));
        let dist_y = Ord::max(0, Ord::max(b.min_point.y - p.y, p.y - b.max_point.y));
        let dist_z = Ord::max(0, Ord::max(b.min_point.z - p.z, p.z - b.max_point.z));
        square(dist_x) + square(dist_y) + square(dist_z)
    }
}

impl Bounds<f32, Vector3f> {
    fn lerp(&self, t: Vector3f) -> Vector3f {
        Vector3f {
            x: math::lerp(t.x, self.min_point.x, self.max_point.x),
            y: math::lerp(t.y, self.min_point.y, self.max_point.y),
            z: math::lerp(t.z, self.min_point.z, self.max_point.z)
        }
    }

    /// Gives the [0, 1] ratio of how far the given point is across the bounding box (i.e. the inverse of lerp).
    fn offset(&self, p: Vector3f) -> Vector3f {
        let mut o = p - self.min_point;
        if self.max_point.x > self.min_point.x {
            o.x /= self.max_point.x - self.min_point.x;
        }
        if self.max_point.y > self.min_point.y {
            o.y /= self.max_point.y - self.min_point.y;
        }
        if self.max_point.z > self.min_point.z {
            o.z /= self.max_point.z - self.min_point.z;
        }
        o
    }
}

impl<T: Num, V: Vector<T>> Index<usize> for Bounds<T, V> {
    type Output = V;

    fn index(&self, index: usize) -> &Self::Output {
        debug_assert!(index <= 1);
        match index {
            0 => &self.min_point,
            1 => &self.max_point,
            _ => panic!("invalid index in bounds struct")
        }
    }
}


pub type Bounds2i = Bounds<i32, Vector2i>;
pub type Bounds2f = Bounds<f32, Vector2f>;
pub type Bounds3i = Bounds<i32, Vector3i>;
pub type BOunds3f = Bounds<f32, Vector3f>;