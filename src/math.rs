use core::f32;
use std::{i32, ops::{Add, Div, Mul, Neg, Sub}};
use crate::geometry::vector::Vector2f;

/// Squares the given number.
#[inline]
pub fn square<T>(a: T) -> T 
    where T: std::ops::Mul<T, Output=T> + Copy 
    {
    a * a
}

/// Interpolates between `a` and `b` at the point `x` (where `x=0` is `a`, `x=1` is `b`).
#[inline]
pub const fn lerp(x: f32, a: f32, b: f32) -> f32 {
    (1.0 - x) * a + x * b
}

/// Evaluates the probability density function for the linear interpolation function between `a` and `b` at point `x`.
#[inline]
pub const fn linear_pdf(x: f32, a: f32, b: f32) -> f32 {
    if x > 1.0 || x < 0.0 {
        return 0.0;
    }
    2.0 * lerp(x, a, b) / (a + b)
}

/// Generates a sample from the PDF of the linear interpolation function, given `c` from the `xi` distribution.
/// Clamps the output to be strictly lower than `1.0`.
#[inline]
pub fn linear_sample(u: f32, a: f32, b: f32) -> f32 {
    if u == 0.0 && a == 0.0 {
        return 0.0;
    }
    let x = u * (a + b) / (a + f32::sqrt(lerp(u, square(a), square(b))));
    f32::min(x, ONE_MINUS_EPSILON)
}

/// Finds the corresponding `u` for the value `x` taken from the linear interpolation function's CDF
#[inline]
pub const fn invert_linear_sample(x: f32, a: f32, b: f32) -> f32 {
    x * (a * (2.0 - x) + (b * x)) / (a + b)
}

/// The bilinear interpolation function, where `w` are the corners of the function and `x, y` interpolate over them in the 
/// unit square.
#[inline]
const fn bilerp(p: Vector2f, w: &[f32; 4]) -> f32 {
    (1.0 - p.x) * (1.0 - p.y) * w[0] + 
    p.x * (1.0 - p.y) * w[1] + 
    p.y * (1.0 - p.x) * w[2] + 
    p.x * p.y * w[3]
} 

/// The bilinear interpolation function's PDF.
#[inline]
pub const fn bilinear_pdf(p: Vector2f, w: &[f32; 4]) -> f32 {
    if p.x < 0.0 || p.x > 1.0 || p.y < 0.0 || p.y > 1.0 {
        return 0.0;
    }
    if w[0] + w[1] + w[2] + w[3] == 0.0 {
        return 1.0;
    }
    4.0 * bilerp(p, w) / (w[0] + w[1] + w[2] + w[3])
}

/// Takes a sample from the bilinear distribution, given `u` in the `xi` x `xi` distribution.
#[inline]
pub fn bilinear_sample(u: Vector2f, w: &[f32; 4]) -> Vector2f {
    let y = linear_sample(u.y, w[0] + w[1], w[2] + w[3]);
    let x = linear_sample(u.x, lerp(y, w[0], w[2]), lerp(y, w[1], w[3]));
    Vector2f {x, y}
}

pub const ONE_MINUS_EPSILON: f32 = 1.0 - f32::EPSILON;

pub fn arcsin(f: f32) -> f32 {
    f32::clamp(f32::asin(f), -1.0, 1.0)
}

pub fn arccos(f: f32) -> f32 {
    f32::clamp(f32::acos(f), -1.0, 1.0)
}

pub trait Num: 
    Add<Self, Output = Self> + Sub<Self, Output = Self> + 
    Mul<Self, Output = Self> + Div<Self, Output = Self> + 
    Copy + Neg<Output=Self> + PartialOrd { 
    fn is_nan(a: Self) -> bool;
    fn sqrt(a: Self) -> f32;
    fn ceil(a: Self) -> Self;
    fn floor(a: Self) -> Self;
    fn to_float(a: Self) -> f32;
    fn abs(a: Self) -> Self;
    fn fma(a: Self, b: Self, c: Self) -> Self;
    fn min(a: Self, b: Self) -> Self;
    fn max(a: Self, b: Self) -> Self;
    fn default() -> Self;
    fn maxiumum() -> Self;
    fn from_f32(x: f32) -> Self;
    fn from_i32(x: i32) -> Self;

    fn diff_products(a: Self, b: Self, c: Self, d: Self) -> Self {
        let cd = c * d;
        let diff_prod = Num::fma(a, b, -cd);
        let err = Num::fma(-c, d, cd);
        diff_prod + err
    }

    fn sum_products(a: Self, b: Self, c: Self, d: Self) -> Self {
        let cd = c * d;
        let sum_prod = Num::fma(a, b, cd);
        let err = Num::fma(c, d, -cd);
        sum_prod + err
    }
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

    fn default() -> Self {
        0
    }

    fn maxiumum() -> Self {
        i32::MAX
    }
    
    fn from_f32(x: f32) -> Self {
        x as i32
    }
    
    fn from_i32(x: i32) -> Self {
        x
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

    fn default() -> Self {
        0.0
    }

    fn maxiumum() -> Self {
        f32::MAX
    }
    
    fn from_f32(x: f32) -> Self {
        x
    }
    
    fn from_i32(x: i32) -> Self {
        x as f32
    }
}

// The following two functions were ripped from https://rust-lang.github.io/rfcs/3173-float-next-up-down.html.
// They're used for floating point error correction.

/// Returns the least number greater than `x`.
///
/// Let `TINY` be the smallest representable positive `f32`. Then,
///  - if `x.is_nan()`, this returns `x`;
///  - if `x` is `NEG_INFINITY`, this returns `-MAX`;
///  - if `x` is `-TINY`, this returns -0.0;
///  - if `x` is -0.0 or +0.0, this returns `TINY`;
///  - if `x` is `MAX` or `INFINITY`, this returns `INFINITY`;
///  - otherwise the unique least value greater than `self` is returned.
///
/// The identity `x.next_up() == -(-x).next_down()` holds for all `x`. When `x`
/// is finite `x == x.next_up().next_down()` also holds.
pub const fn next_up(x: f32) -> f32 {
    const TINY_BITS: u32 = 0x1; // Smallest positive f32.
    const CLEAR_SIGN_MASK: u32 = 0x7fff_ffff;

    let bits = x.to_bits();
    if x.is_nan() || bits == f32::INFINITY.to_bits() {
        return x;
    }
    
    let abs = bits & CLEAR_SIGN_MASK;
    let next_bits = if abs == 0 {
        TINY_BITS
    } else if bits == abs {
        bits + 1
    } else {
        bits - 1
    };
    f32::from_bits(next_bits)
}

/// Returns the greatest number less than `self`.
///
/// Let `TINY` be the smallest representable positive `f32`. Then,
///  - if `self.is_nan()`, this returns `self`;
///  - if `self` is `INFINITY`, this returns `MAX`;
///  - if `self` is `TINY`, this returns 0.0;
///  - if `self` is -0.0 or +0.0, this returns `-TINY`;
///  - if `self` is `-MAX` or `NEG_INFINITY`, this returns `NEG_INFINITY`;
///  - otherwise the unique greatest value less than `self` is returned.
///
/// The identity `x.next_down() == -(-x).next_up()` holds for all `x`. When `x`
/// is finite `x == x.next_down().next_up()` also holds.
pub const fn next_down(x: f32) -> f32 {
    const NEG_TINY_BITS: u32 = 0x8000_0001; // Smallest (in magnitude) negative f32.
    const CLEAR_SIGN_MASK: u32 = 0x7fff_ffff;

    let bits = x.to_bits();
    if x.is_nan() || bits == f32::NEG_INFINITY.to_bits() {
        return x;
    }
    
    let abs = bits & CLEAR_SIGN_MASK;
    let next_bits = if abs == 0 {
        NEG_TINY_BITS
    } else if bits == abs {
        bits - 1
    } else {
        bits + 1
    };
    f32::from_bits(next_bits)
}
