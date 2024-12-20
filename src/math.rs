use core::f32;
use std::{f32::consts::PI, i32, ops::{Add, Div, Mul, Neg, Sub}};
use crate::geometry::vector::{Vector, Vector2f, Vector3, Vector3f};

/// Squares the given number.
#[inline]
pub fn square<T: Num>(a: T) -> T {
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

/// The standard `f32::asin` function, but clamped between -1 and 1 for float purposes.
/// I'm not entirely sure if this is necessary, but pbrt uses something similar for its
/// arcsin implementation.
pub fn arcsin(f: f32) -> f32 {
    f32::clamp(f32::asin(f), -1.0, 1.0)
}

/// see `arcsin` for justtification on what this function does
pub fn arccos(f: f32) -> f32 {
    f32::clamp(f32::acos(f), -1.0, 1.0)
}

pub fn safe_sqrt(f: f32) -> f32 {
    f32::sqrt(f32::max(0.0, f))
}

/// Provides all of the traits that a number should have for the purposes of this library. 
/// Allows for split implementation of types across traits (see: `Vector`) where guaranteeing that
/// all of these operations are available is desirable.
/// 
/// Currently implemented for the types
/// - `f32`
/// - `i32`
pub trait Num: 
    Add<Self, Output = Self> + Sub<Self, Output = Self> + 
    Mul<Self, Output = Self> + Div<Self, Output = Self> + 
    Copy + Neg<Output=Self> + PartialOrd + 
    Mul<i32, Output = Self> + Mul<f32, Output = Self> { 
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

impl Mul<f32> for i32 {
    type Output = i32;

    fn mul(self, rhs: f32) -> Self::Output {
        todo!()
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

/// Computes the product of the two given numbers, where the first return value
/// is the return value and the second is the error incurred.
fn compensated_prod(a: f32, b: f32) -> (f32, f32) {
    (a * b, f32::fma(a, b, -(a * b)))
}

fn compensated_sum(a: f32, b: f32) -> (f32, f32) {
    let s = a + b;
    let delta = s - a;
    (s, (a - (s - delta)) + (b - delta))
}

pub fn compensated_inner_prod(a: &[f32]) -> (f32, f32) {
    debug_assert!(a.len() >= 2 && a.len() % 2 == 0);
    if a.len() == 2 {
        return compensated_prod(a[0], a[1]);
    } else if a.len() == 1 {
        return (a[0], 0.0);
    }

    let (ab_v, ab_err) = compensated_prod(a[0], a[1]);
    let (tp_v, tp_err) = compensated_inner_prod(&a[2..]);
    let (sum_v, sum_err) = compensated_sum(ab_v, tp_v);
    (sum_v, ab_err + (tp_err + sum_err))
}

pub fn inner_prod(a: &[f32]) -> f32 {
    let (res, err) = compensated_inner_prod(a);
    res + err
}

/// Computes the (solid angle) area of the spherical triangle with vertices given by `a`, `b`, and `c`.
pub fn spherical_triangle_area(a: Vector3f, b: Vector3f, c: Vector3f) -> f32 {
    f32::abs(2.0 * f32::atan2(Vector3f::dot(a, Vector3f::cross(b, c)), 
        1.0 + Vector3f::dot(a, b) + Vector3f::dot(a, c) + Vector3f::dot(b, c)))
}

/// Computes the (solid angle) area of the spherical quadrangle with vertices given by `a, b, c, d`.
pub fn spherical_quad_area(a: Vector3f, b: Vector3f, c: Vector3f, d: Vector3f) -> f32 {
    let mut a_cross_b = Vector3f::cross(a, b);
    let mut b_cross_c = Vector3f::cross(b, c);
    let mut c_cross_d = Vector3f::cross(c, d);
    let mut d_cross_a = Vector3f::cross(d, a);
    if Vector3::length_squared(a_cross_b) == 0.0 || Vector3::length_squared(b_cross_c) == 0.0 || 
    Vector3::length_squared(c_cross_d) == 0.0 || Vector3::length_squared(d_cross_a) == 0.0 {
        return 0.0;
    }
    a_cross_b = Vector3f::normalize(a_cross_b);
    b_cross_c = Vector3f::normalize(b_cross_c);
    c_cross_d = Vector3f::normalize(c_cross_d);
    d_cross_a = Vector3f::normalize(d_cross_a);

    let alpha = Vector3f::angle_between(d_cross_a, -a_cross_b);
    let beta = Vector3f::angle_between(a_cross_b, -b_cross_c);
    let gamma = Vector3f::angle_between(b_cross_c, -c_cross_d);
    let delta = Vector3f::angle_between(c_cross_d, -d_cross_a);

    f32::abs(alpha + beta + gamma + delta - 2.0 * PI)
}

/// Converts a set of spherical coordinates into a unit vector in (x, y, z).
/// Takes in an already computed `sin(Theta)`, `cos(Theta)`, but not `phi`.
pub fn spherical_direction(sin_theta: f32, cos_theta: f32, phi: f32) -> Vector3f {
    Vector3f::new(
        f32::clamp(sin_theta, -1.0, 1.0) * f32::cos(phi),
        f32::clamp(sin_theta, -1.0, 1.0) * f32::sin(phi),
        f32::clamp(cos_theta, -1.0, 1.0)
    )
}

/// Computes `theta` in spherical coordinates, given a unit vector in `x, y, z`.
pub fn spherical_theta(v: Vector3f) -> f32 {
    arccos(v.z)
}

/// Computes `phi` in spherical coordinates, given a unit vector in `x, y, z`.
pub fn spherical_phi(v: Vector3f) -> f32 {
    let p = f32::atan2(v.y, v.x);
    if p < 0.0 {
        p + 2.0 * PI
    } else {
        p
    }
}

/// Gets the value of `cos(theta)` from spherical coordinates, given the vector representation of it.
pub fn cos_theta(w: Vector3f) -> f32 {
    w.z
}

/// Gets the value of `cos^2(theta)` from spherical coordinates, given the vector representation of it.
pub fn cos2_theta(w: Vector3f) -> f32 {
    square(w.z)
}

/// Gets the value of `abs(cos(theta))` from spherical coordinates, given the vector representation of it.
pub fn abs_cos_theta(w: Vector3f) -> f32 {
    f32::abs(w.z)
}

/// Gets the value of `sin^2(theta)` from spherical coordinates, given the vector representation of it.
pub fn sin2_theta(w: Vector3f) -> f32 {
    f32::max(0.0, 1.0 - cos2_theta(w))
}

/// Gets the value of `sin(theta)` from spherical coordinates, given the vector representation of it.
pub fn sin_theta(w: Vector3f) -> f32 {
    f32::sqrt(sin2_theta(w))
}

pub fn tan_theta(w: Vector3f) -> f32 {
    sin_theta(w) / cos_theta(w)
}

pub fn tan2_theta(w: Vector3f) -> f32 {
    sin2_theta(w) / cos2_theta(w)
}

pub fn cos_phi(w: Vector3f) -> f32 {
    let sin_theta = sin_theta(w);
    if sin_theta == 0.0 {
        1.0
    } else {
        f32::clamp(w.x / sin_theta, -1.0, 1.0)
    }
}

pub fn sin_phi(w: Vector3f) -> f32 {
    let sin_theta = sin_theta(w);
    if sin_theta == 0.0 {
        0.0
    } else {
        f32::clamp(w.y / sin_theta, -1.0, 1.0)
    }
}

/// Finds the cosine of the angle (delta phi) between two different vectors' phi values.
pub fn cos_d_phi(wa: Vector3f, wb: Vector3f) -> f32 {
    let waxy = square(wa.x) + square(wa.y);
    let wbxy = square(wb.x) + square(wb.y);
    if waxy == 0.0 || wbxy == 0.0 {
        1.0
    } else {
        f32::clamp((wa.x + wb.x + wa.y * wb.y) / f32::sqrt(waxy * wbxy), -1.0, 1.0)
    }
}

pub fn eval_polynomial(t: f32, c: &[f32]) -> f32 {
    if c.len() == 1 {
        return c[0];
    }
    Num::fma(t, eval_polynomial(t, &c[1..]), c[0])
}

pub fn radians(deg: f32) -> f32 {
    (PI / 180.0) * deg
}

pub fn degrees(rad: f32) -> f32 {
    (180.0 / PI) * rad
}