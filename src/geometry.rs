use crate::math::square;
use std::ops::{Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Rem, Sub, SubAssign};

trait Num<Rhs = Self, Output = Self>: 
    Add<Rhs, Output = Output> + Sub<Rhs, Output = Output> + 
    Mul<Rhs, Output = Output> + Div<Rhs, Output = Output> + 
    Rem<Rhs, Output = Output> + Copy + Neg<Output=Self> + PartialOrd { 
    fn is_nan(a: Self) -> bool;
    fn sqrt(a: Self) -> f32;
    fn ceil(a: Self) -> Self;
    fn floor(a: Self) -> Self;
    fn to_float(a: Self) -> f32;
    fn abs(a: Self) -> Self;
    fn fma(a: Self, b: Self, c: Self) -> Self;
    fn min(a: Self, b: Self) -> Self;
    fn max(a: Self, b: Self) -> Self;
}

trait VecOps<T: Num, Rhs = Self, Output = Self>: 
    Add<Rhs, Output = Output> + Sub<Rhs, Output = Output> + 
    AddAssign<Rhs> + SubAssign<Rhs> + Neg + 
    Mul<T, Output = Output> + Div<T, Output = Output> +
    MulAssign<T> + DivAssign<T> { }

/// Defines a vector type to be used in the library. All of these methods are defined over any Vect type, specifically
/// Vector3 and Vector2.

pub trait Vect<T: Num>: Index<isize> {
    type FloatType;
    fn length(a: &Self) -> f32;
    fn length_squared(a: &Self) -> T;
    fn normalize(a: &Self) -> Self::FloatType;
    fn abs(a: &Self) -> Self;
    fn ceil(a: &Self) -> Self;
    fn floor(a: &Self) -> Self;
    fn as_floats(a: &Self) -> Self::FloatType;
    fn lerp(t: f32, a: &Self, b: &Self) -> Self::FloatType;
    fn fma(a: &Self, b: &Self, c: &Self) -> Self;
    fn min(a: &Self, b: &Self) -> Self;
    fn max(a: &Self, b: &Self) -> Self;
    fn min_component(a: &Self) -> T;
    fn max_component(a: &Self) -> T;
    fn min_component_index(a: &Self) -> isize;
    fn max_component_index(a: &Self) -> isize;
    fn permute(a: &Self, pm: &[isize]) -> Self;
    fn h_prod(a: &Self) -> T;
}


impl Num for i32 {
    fn is_nan(_: Self) -> bool {
        false
    }

    fn sqrt(a: Self) -> f32 {
        (a as f32).sqrt()
    }

    fn ceil(a: Self) -> Self {
        a
    }

    fn floor(a: Self) -> Self {
        a
    }

    fn to_float(a: Self) -> f32 {
        a as f32
    }

    fn abs(a: Self) -> Self {
        i32::abs(a)
    }
    
    fn fma(a: Self, b: Self, c: Self) -> Self {
        a*b+c
    }
    
    fn min(a: Self, b: Self) -> Self {
        Ord::min(a, b)
    }
    
    fn max(a: Self, b: Self) -> Self {
        Ord::max(a, b)
    }
}

impl Num for f32 {
    fn to_float(a: Self) -> f32 {
        a
    }
    
    fn is_nan(a: Self) -> bool {
        f32::is_nan(a)
    }
    
    fn sqrt(a: Self) -> f32 {
        f32::sqrt(a)
    }
    
    fn ceil(a: Self) -> Self {
        f32::ceil(a)
    }
    
    fn floor(a: Self) -> Self {
        f32::floor(a)
    }

    fn abs(a: Self) -> Self {
        f32::abs(a)
    }
    
    fn fma(a: Self, b: Self, c: Self) -> Self {
        f32::mul_add(a, b, c)
    }
    
    fn min(a: Self, b: Self) -> Self {
        f32::min(a, b)
    }
    
    fn max(a: Self, b: Self) -> Self {
        f32::max(a, b)
    }
}

#[derive(Copy, Clone)]
struct Vector2<T: Num>{
    pub x: T,
    pub y: T
}

impl<T: Num> Vector2<T> {
    pub fn new(x: T, y: T) -> Self {
        debug_assert!(!Num::is_nan(x) && !Num::is_nan(y));
        Self { x, y }
    }

    pub fn has_nan(&self) -> bool {
        Num::is_nan(self.x) || Num::is_nan(self.y)
    }
}

impl<T: Num> Add<Self> for Vector2<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Vector2 {
            x: self.x + rhs.x,
            y: self.y + rhs.y
        }
    }
}

impl<T: Num> Sub<Self> for Vector2<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Vector2 {
            x: self.x - rhs.x,
            y: self.y - rhs.y
        }
    }
}

impl<T: Num> AddAssign<Self> for Vector2<T> {
    
    fn add_assign(&mut self, rhs: Self) {
        self.x = self.x + rhs.x;
        self.y = self.y + rhs.y;
    }
}

impl<T: Num> SubAssign<Self> for Vector2<T> {

    fn sub_assign(&mut self, rhs: Self) {
        self.x = self.x - rhs.x;
        self.y = self.y - rhs.y;
    }
}

impl<T: Num> Neg for Vector2<T> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Vector2 {
            x: -self.x,
            y: -self.y
        }
    }
}

impl<T:Num> Mul<T> for Vector2<T> {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        Vector2 {
            x: self.x * rhs,
            y: self.y * rhs
        }
    }
}

impl<T:Num> Div<T> for Vector2<T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        Vector2 {
            x: self.x / rhs,
            y: self.y / rhs
        }
    }
}

impl<T: Num> MulAssign<T> for Vector2<T> {
    fn mul_assign(&mut self, rhs: T) {
        self.x = self.x * rhs;
        self.y = self.y * rhs;
    }
}

impl<T: Num> DivAssign<T> for Vector2<T> {
    fn div_assign(&mut self, rhs: T) {
        self.x = self.x / rhs;
        self.y = self.y / rhs;
    }
}

impl<T: Num> Index<isize> for Vector2<T> {
    type Output = T;

    fn index(&self, index: isize) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            _ => panic!("Index out of bounds in Vector2")
        }
    }
}

impl<T: Num> Vect<T> for Vector2<T> {

    type FloatType = Vector2<f32>;

    fn length_squared(a: &Self) -> T {
        square(a.x) + square(a.y)
    }

    fn length(a: &Self) -> f32 {
        Num::sqrt(Vect::length_squared(a))
    }

    fn normalize(a: &Self) -> Vector2<f32> {
        Vect::as_floats(a) / Vect::length(a)
    }
    
    fn abs(a: &Self) -> Self {
        Vector2 { 
            x: Num::abs(a.x), 
            y: Num::abs(a.y) 
        }
    }
    
    fn ceil(a: &Self) -> Self {
        Vector2 { 
            x: Num::ceil(a.x), 
            y: Num::ceil(a.y)
        }
    }
    
    fn floor(a: &Self) -> Self {
        Vector2 { 
            x: Num::floor(a.x), 
            y: Num::floor(a.y)
        }
    }

    fn as_floats(a: &Self) -> Self::FloatType {
        Vector2 { 
            x: Num::to_float(a.x), 
            y: Num::to_float(a.y)
        }
    }
    
    fn lerp(t: f32, a: &Self, b: &Self) -> Self::FloatType {
        Vect::as_floats(a) * (1.0 - t) + Vect::as_floats(b) * t
    }
    
    fn fma(a: &Self, b: &Self, c: &Self) -> Self {
        Vector2 {
            x: Num::fma(a.x, b.x, c.x),
            y: Num::fma(a.y, b.y, c.y)
        }
    }
    
    fn min(a: &Self, b: &Self) -> Self {
        Vector2 {
            x: Num::min(a.x, b.x),
            y: Num::min(a.y, b.y)
        }
    }
    
    fn max(a: &Self, b: &Self) -> Self {
        Vector2 {
            x: Num::max(a.x, b.x),
            y: Num::max(a.y, b.y)
        }
    }
    
    fn min_component(a: &Self) -> T {
        Num::min(a.x, a.y)
    }
    
    fn max_component(a: &Self) -> T {
        Num::max(a.x, a.y)
    }
    
    fn min_component_index(a: &Self) -> isize {
        let m = Vect::min_component(a);
        if m == a.x { 0 } else { 1 }
    }
    
    fn max_component_index(a: &Self) -> isize {
        let m = Vect::max_component(a);
        if m == a.x { 0 } else { 1 }
    }
    
    fn permute(a: &Self, pm: &[isize]) -> Self {
        Vector2 {
            x: a[pm[0]],
            y: a[pm[1]]
        }
    }
    
    fn h_prod(a: &Self) -> T {
        a.x * a.y
    }


}


struct Vector3<T: Num> {
    pub x: T,
    pub y: T,
    pub z: T
}

impl<T: Num> Vector3<T> {

}

impl<T: Num> Add<Self> for Vector3<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Vector3 {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z
        }
    }
}

impl<T: Num> Sub<Self> for Vector3<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Vector3 {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z
        }
    }
}

impl<T: Num> AddAssign<Self> for Vector3<T> {
    
    fn add_assign(&mut self, rhs: Self) {
        self.x = self.x + rhs.x;
        self.y = self.y + rhs.y;
        self.z = self.z + rhs.z;
    }
}

impl<T: Num> SubAssign<Self> for Vector3<T> {

    fn sub_assign(&mut self, rhs: Self) {
        self.x = self.x - rhs.x;
        self.y = self.y - rhs.y;
        self.z = self.z - rhs.z;
    }
}

impl<T: Num> Neg for Vector3<T> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Vector3 {
            x: -self.x,
            y: -self.y,
            z: -self.z
        }
    }
}

impl<T:Num> Mul<T> for Vector3<T> {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        Vector3 {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs
        }
    }
}

impl<T:Num> Div<T> for Vector3<T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        Vector3 {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs
        }
    }
}

impl<T: Num> MulAssign<T> for Vector3<T> {
    fn mul_assign(&mut self, rhs: T) {
        self.x = self.x * rhs;
        self.y = self.y * rhs;
        self.z = self.z * rhs;
    }
}

impl<T: Num> DivAssign<T> for Vector3<T> {
    fn div_assign(&mut self, rhs: T) {
        self.x = self.x / rhs;
        self.y = self.y / rhs;
        self.z = self.z / rhs;
    }
}

impl<T: Num> Index<isize> for Vector3<T> {
    type Output = T;

    fn index(&self, index: isize) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Index out of bounds in Vector3")
        }
    }
}

impl<T: Num> Vect<T> for Vector3<T> {
    type FloatType = Vector3<f32>;

    fn length(a: &Self) -> f32 {
        Num::sqrt(Vect::length_squared(a))
    }

    fn length_squared(a: &Self) -> T {
        square(a.x) + square(a.y) + square(a.z)
    }

    fn normalize(a: &Self) -> Self::FloatType {
        Vect::as_floats(a) / Vect::length(a)
    }

    fn abs(a: &Self) -> Self {
        Vector3 {
            x: Num::abs(a.x),
            y: Num::abs(a.y),
            z: Num::abs(a.z)
        }
    }

    fn ceil(a: &Self) -> Self {
        Vector3 {
            x: Num::ceil(a.x),
            y: Num::ceil(a.y),
            z: Num::ceil(a.z)
        }
    }

    fn floor(a: &Self) -> Self {
        Vector3 {
            x: Num::floor(a.x),
            y: Num::floor(a.y),
            z: Num::floor(a.z)
        }
    }

    fn as_floats(a: &Self) -> Self::FloatType {
        Vector3 {
            x: Num::to_float(a.x),
            y: Num::to_float(a.y),
            z: Num::to_float(a.z)
        }
    }

    fn lerp(t: f32, a: &Self, b: &Self) -> Self::FloatType {
        Vect::as_floats(a) * (1.0 - t) + Vect::as_floats(b) * t
    }

    fn fma(a: &Self, b: &Self, c: &Self) -> Self {
        Vector3 {
            x: Num::fma(a.x, b.x, c.x),
            y: Num::fma(a.y, b.y, c.y),
            z: Num::fma(a.z, b.z, c.z)
        }
    }

    fn min(a: &Self, b: &Self) -> Self {
        Vector3 {
            x: Num::min(a.x, b.x),
            y: Num::min(a.y, b.y),
            z: Num::min(a.z, b.z)
        }
    }

    fn max(a: &Self, b: &Self) -> Self {
        Vector3 {
            x: Num::max(a.x, b.x),
            y: Num::max(a.y, b.y),
            z: Num::max(a.z, b.z)
        }
    }

    fn min_component(a: &Self) -> T {
        Num::min(Num::min(a.x, a.y), a.z)
    }

    fn max_component(a: &Self) -> T {
        Num::max(Num::max(a.x, a.y), a.z)
    }

    fn min_component_index(a: &Self) -> isize {
        let m = Vect::min_component(a);
        if m == a.x {0} else if m == a.y {1} else {2}
    }

    fn max_component_index(a: &Self) -> isize {
        let m = Vect::max_component(a);
        if m == a.x {0} else if m == a.y {1} else {2}
    }

    fn permute(a: &Self, pm: &[isize]) -> Self {
        Vector3 {
            x: a[pm[0]],
            y: a[pm[1]],
            z: a[pm[2]]
        }
    }

    fn h_prod(a: &Self) -> T {
        a.x * a.y * a.z
    }
}

pub type Vector2f = Vector2<f32>;
pub type Vector2i = Vector2<i32>;
pub type Vector3f = Vector3<f32>;
pub type Vector3i = Vector3<i32>;