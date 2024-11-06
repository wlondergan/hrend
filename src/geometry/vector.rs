use crate::math::{arccos, arcsin, degrees, eval_polynomial, safe_sqrt, square, Num};
use std::{f32::consts::PI, ops::{Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Sub, SubAssign}};
use std::mem;
use super::bounds::Bounds3f;

trait VecOps<T: Num, Rhs = Self, Output = Self>: 
    Add<Rhs, Output = Output> + Sub<Rhs, Output = Output> + 
    AddAssign<Rhs> + SubAssign<Rhs> + Neg<Output = Output> + 
    Mul<T, Output = Output> + Div<T, Output = Output> +
    MulAssign<T> + DivAssign<T> { }

/// Defines a vector type to be used in the library. All of these methods are defined over any V type, specifically
/// Vector3 and Vector2.
pub trait Vector<T: Num>: Index<usize> + VecOps<T> + Sized + Copy {
    
    type FloatType: Vector<f32>;

    fn default() -> Self;

    fn maximum() -> Self;

    fn all_x(x: T) -> Self;
    
    fn length(a: Self) -> f32;
    
    fn length_squared(a: Self) -> T;
    
    fn normalize(a: Self) -> Self::FloatType;
    
    fn abs(a: Self) -> Self;
    
    fn ceil(a: Self) -> Self;
    
    fn floor(a: Self) -> Self;
    
    fn as_floats(a: Self) -> Self::FloatType;

    fn lerp(t: f32, a: Self, b: Self) -> Self::FloatType;

    fn fma(a: Self, b: Self, c: Self) -> Self;

    fn min(a: Self, b: Self) -> Self;

    fn max(a: Self, b: Self) -> Self;

    fn min_component(a: Self) -> T;

    fn max_component(a: Self) -> T;

    fn min_component_index(a: Self) -> usize;

    fn max_component_index(a: Self) -> usize;

    fn permute(a: Self, pm: &[usize]) -> Self;

    fn h_prod(a: Self) -> T;

    fn dot(a: Self, b: Self) -> T;

    fn angle_between(a: Self, b: Self) -> f32 {
        if Self::FloatType::dot(Vector::as_floats(a), Vector::as_floats(b)) < 0.0 {
            return std::f32::consts::PI - 2.0 * arcsin(Vector::length(a + b));
        } else {
            return 2.0 * arcsin(Vector::length(a + b));
        }
    }

    fn abs_dot(a: Self, b: Self) -> T {
        Num::abs(Vector::dot(a, b))
    }

    fn gram_schmidt(v: Self, w: Self) -> Self {
        v - w * Vector::dot(v, w)
    }

    fn distance(p1: Self, p2: Self) -> f32 {
        Vector::length(p1 - p2)
    }

    fn distance_squared(p1: Self, p2: Self) -> f32 {
        Num::to_float(Vector::length_squared(p1 - p2))
    }

    fn face_forward(n: Self, v: Self) -> Self {
        if Num::to_float(Vector::dot(n, v)) < 0.0 {
            n
        } else {
            -n
        }
    }

}


#[derive(Copy, Clone)]
pub struct Vector2<T: Num>{
    pub x: T,
    pub y: T
}

impl<T: Num> VecOps<T, Self, Self> for Vector2<T> { }

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

impl<T: Num> Index<usize> for Vector2<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            _ => panic!("Index out of bounds in Vector2")
        }
    }
}

impl<T: Num> Vector<T> for Vector2<T> {

    type FloatType = Vector2<f32>;

    fn default() -> Self {
        Vector2 {
            x: T::default(),
            y: T::default()
        }
    }

    fn maximum() -> Self {
        Vector2 {
            x: T::default(),
            y: T::default()
        }
    }

    fn length_squared(a: Self) -> T {
        square(a.x) + square(a.y)
    }

    fn length(a: Self) -> f32 {
        Num::sqrt(Vector::length_squared(a))
    }

    fn normalize(a: Self) -> Vector2<f32> {
        Vector::as_floats(a) / Vector::length(a)
    }
    
    fn abs(a: Self) -> Self {
        Vector2 { 
            x: Num::abs(a.x), 
            y: Num::abs(a.y) 
        }
    }
    
    fn ceil(a: Self) -> Self {
        Vector2 { 
            x: Num::ceil(a.x), 
            y: Num::ceil(a.y)
        }
    }
    
    fn floor(a: Self) -> Self {
        Vector2 { 
            x: Num::floor(a.x), 
            y: Num::floor(a.y)
        }
    }

    fn as_floats(a: Self) -> Self::FloatType {
        Vector2 { 
            x: Num::to_float(a.x), 
            y: Num::to_float(a.y)
        }
    }
    
    fn lerp(t: f32, a: Self, b: Self) -> Self::FloatType {
        Vector::as_floats(a) * (1.0 - t) + Vector::as_floats(b) * t
    }
    
    fn fma(a: Self, b: Self, c: Self) -> Self {
        Vector2 {
            x: Num::fma(a.x, b.x, c.x),
            y: Num::fma(a.y, b.y, c.y)
        }
    }
    
    fn min(a: Self, b: Self) -> Self {
        Vector2 {
            x: Num::min(a.x, b.x),
            y: Num::min(a.y, b.y)
        }
    }
    
    fn max(a: Self, b: Self) -> Self {
        Vector2 {
            x: Num::max(a.x, b.x),
            y: Num::max(a.y, b.y)
        }
    }
    
    fn min_component(a: Self) -> T {
        Num::min(a.x, a.y)
    }
    
    fn max_component(a: Self) -> T {
        Num::max(a.x, a.y)
    }
    
    fn min_component_index(a: Self) -> usize {
        let m = Vector::min_component(a);
        if m == a.x { 0 } else { 1 }
    }
    
    fn max_component_index(a: Self) -> usize {
        let m = Vector::max_component(a);
        if m == a.x { 0 } else { 1 }
    }
    
    fn permute(a: Self, pm: &[usize]) -> Self {
        Vector2 {
            x: a[pm[0]],
            y: a[pm[1]]
        }
    }
    
    fn h_prod(a: Self) -> T {
        a.x * a.y
    }

    fn dot(a: Self, b: Self) -> T {
        a.x * b.x + a.y * b.y
    }
    
    fn all_x(x: T) -> Self {
        Self {
            x,
            y: x
        }
    }


}

#[derive(Copy, Clone)]
pub struct Vector3<T: Num> {
    pub x: T,
    pub y: T,
    pub z: T
}

impl<T: Num> Vector3<T> {

    pub fn new(x: T, y: T, z: T) -> Self {
        Vector3 {x, y, z}
    }

    pub fn cross(v: Self, w: Self) -> Self {
        Vector3 {
            x: Num::diff_products(v.y, v.z, w.z, w.y),
            y: Num::diff_products(v.z, w.x, v.x, w.z),
            z: Num::diff_products(v.x, w.y, v.y, w.x)
        }
    }

    pub fn coord_system(v: Self) -> (Vector3f, Vector3f) {
        let vf = Vector::as_floats(v);
        let sign = f32::copysign(1.0, vf.z);
        let a = -1.0 / (sign + vf.z);
        let b = vf.x * vf.y * a;
        (
            Vector3::new(1.0 + sign * square(vf.x) * a, sign * b, -sign * vf.x),
            Vector3::new(b, sign + square(vf.y) * a, -vf.y)
        )
    }
}

impl<T: Num> VecOps<T, Self, Self> for Vector3<T> { }

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

impl<T: Num> Index<usize> for Vector3<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Index out of bounds in Vector3")
        }
    }
}

impl<T: Num> Vector<T> for Vector3<T> {
    type FloatType = Vector3<f32>;

    fn length(a: Self) -> f32 {
        Num::sqrt(Vector::length_squared(a))
    }

    fn length_squared(a: Self) -> T {
        square(a.x) + square(a.y) + square(a.z)
    }

    fn normalize(a: Self) -> Self::FloatType {
        Vector::as_floats(a) / Vector::length(a)
    }

    fn abs(a: Self) -> Self {
        Vector3 {
            x: Num::abs(a.x),
            y: Num::abs(a.y),
            z: Num::abs(a.z)
        }
    }

    fn ceil(a: Self) -> Self {
        Vector3 {
            x: Num::ceil(a.x),
            y: Num::ceil(a.y),
            z: Num::ceil(a.z)
        }
    }

    fn floor(a: Self) -> Self {
        Vector3 {
            x: Num::floor(a.x),
            y: Num::floor(a.y),
            z: Num::floor(a.z)
        }
    }

    fn as_floats(a: Self) -> Self::FloatType {
        Vector3 {
            x: Num::to_float(a.x),
            y: Num::to_float(a.y),
            z: Num::to_float(a.z)
        }
    }

    fn lerp(t: f32, a: Self, b: Self) -> Self::FloatType {
        Vector::as_floats(a) * (1.0 - t) + Vector::as_floats(b) * t
    }

    fn fma(a: Self, b: Self, c: Self) -> Self {
        Vector3 {
            x: Num::fma(a.x, b.x, c.x),
            y: Num::fma(a.y, b.y, c.y),
            z: Num::fma(a.z, b.z, c.z)
        }
    }

    fn min(a: Self, b: Self) -> Self {
        Vector3 {
            x: Num::min(a.x, b.x),
            y: Num::min(a.y, b.y),
            z: Num::min(a.z, b.z)
        }
    }

    fn max(a: Self, b: Self) -> Self {
        Vector3 {
            x: Num::max(a.x, b.x),
            y: Num::max(a.y, b.y),
            z: Num::max(a.z, b.z)
        }
    }

    fn min_component(a: Self) -> T {
        Num::min(Num::min(a.x, a.y), a.z)
    }

    fn max_component(a: Self) -> T {
        Num::max(Num::max(a.x, a.y), a.z)
    }

    fn min_component_index(a: Self) -> usize {
        let m = Vector::min_component(a);
        if m == a.x {0} else if m == a.y {1} else {2}
    }

    fn max_component_index(a: Self) -> usize {
        let m = Vector::max_component(a);
        if m == a.x {0} else if m == a.y {1} else {2}
    }

    fn permute(a: Self, pm: &[usize]) -> Self {
        Vector3 {
            x: a[pm[0]],
            y: a[pm[1]],
            z: a[pm[2]]
        }
    }

    fn h_prod(a: Self) -> T {
        a.x * a.y * a.z
    }

    fn dot(a: Self, b: Self) -> T {
        a.x * b.x + a.y * b.y + a.z * b.z
    }
    
    fn default() -> Self {
        Vector3 {
            x: T::default(),
            y: T::default(),
            z: T::default()
        }
    }
    
    fn maximum() -> Self {
        Vector3 {
            x: T::maxiumum(),
            y: T::maxiumum(),
            z: T::maxiumum()
        }
    }
    
    fn all_x(x: T) -> Self {
        Self {
            x,
            y: x,
            z: x
        }
    }
}

pub type Vector2f = Vector2<f32>;
pub type Vector2i = Vector2<i32>;
pub type Vector3f = Vector3<f32>;
pub type Vector3i = Vector3<i32>;

pub type Normal3f = Vector3f;

/// Encodes the rotation of a Vector3f projected onto the unit sphere using two `u16`, with a mapping from the unit 
/// sphere projected onto an octahedron. Use `Vector3f as OctahedralVector` or `OctahedralVector as Vector3f` to achieve
/// this conversion.
pub struct OctahedralVector {
    x: u16,
    y: u16
}

const U16_MAX_FLOAT: f32 = 65535.0;

/// Takes a float in the range of [-1, 1] to an unsigned integer,
/// corresponding to a discretization of that range into 2^16 discrete segments.
fn ovec_float_encode(f: f32) -> u16 {
    f32::round(f32::clamp((f + 1.0) / 2.0, 0.0, 1.0) * U16_MAX_FLOAT) as u16
}

impl From<Vector3f> for OctahedralVector {
    fn from(value: Vector3f) -> Self {
        let mut v = value;
        v /= f32::abs(v.x) + f32::abs(v.y) + f32::abs(v.z);
        if v.z >= 0.0 {
            OctahedralVector {
                x: ovec_float_encode(v.x),
                y: ovec_float_encode(v.y)
            }
        } else {
         OctahedralVector {
                x: ovec_float_encode((1.0 - f32::abs(v.y)) * f32::signum(v.x)),
                y: ovec_float_encode((1.0 - f32::abs(v.x)) * f32::signum(v.y))
            }
        }
    }
}

impl Into<Vector3f> for OctahedralVector {
    fn into(self) -> Vector3f {
        let x = -1.0 + 2.0 * (self.x as f32 / U16_MAX_FLOAT);
        let y = -1.0 + 2.0 * (self.y as f32 / U16_MAX_FLOAT);
        let mut v = Vector3f {
            x,
            y,
            z: 1.0 - (f32::abs(x) + f32::abs(y))
        };
        if v.z < 0.0 {
            let vx = v.x;
            v.x = 1.0 - f32::abs(v.y) * f32::signum(vx);
            v.y = 1.0 - f32::abs(vx) * f32::signum(v.y);
        }
        Vector::normalize(v)
    }
}

/// Takes a point from the equal-area mapping onto a Vector3 on the unit sphere.
/// See: pbrt 4th edition, 3.9 for algorithm
/// From pbrt source code: -> Clarberg ~ Fast Equal-Area Mapping of the (Hemi)Sphere using SIMD
pub fn equal_area_square_to_sphere(p: Vector2f) -> Vector3f {
    debug_assert!(p.x >= 0.0 && p.x <= 1.0 && p.y >= 0.0 && p.y <= 1.0);
    let u = 2.0 * p.x - 1.0;
    let v = 2.0 * p.y - 1.0;
    let u_p = f32::abs(u);
    let v_p = f32::abs(v);

    let signed_dist = 1.0 - (u_p + v_p);
    let d = f32::abs(signed_dist);
    let r = 1.0 - d;

    let phi = (if r == 0.0 {1.0} else {(v_p - u_p) / r + 1.0}) * PI / 4.0;

    let z = f32::copysign(1.0 - square(r), signed_dist);

    let cos_phi = f32::copysign(f32::cos(phi), u);
    let sin_phi = f32::copysign(f32::sin(phi), v);
    Vector3f {
        x: cos_phi * r * safe_sqrt(2.0 - square(r)),
        y: sin_phi * r * safe_sqrt(2.0 - square(r)),
        z
    }
}

/// Takes a point from the unit sphere onto the equal-area mapped disk. See: 
/// `equal_area_square_to_sphere`
pub fn equal_area_sphere_to_square(d: Vector3f) -> Vector2f {
    debug_assert!(Vector::length_squared(d) > 0.999 && Vector::length_squared(d) < 1.001);
    let x = f32::abs(d.x);
    let y = f32::abs(d.y);
    let z = f32::abs(d.z);

    let r = safe_sqrt(1.0 - z);
    let a = f32::max(x, y);
    let mut b = f32::min(x, y);
    b = if a == 0.0 {0.0} else {b/a};

    // coefficients for a 6th degree approximation of the function `atan(x) * 2 / pi`
    const COEFFS: [f32;7] = 
    [0.406758566246788489601959989e-5,
    0.636226545274016134946890922156,
    0.61572017898280213493197203466e-2,
    -0.247333733281268944196501420480,
    0.881770664775316294736387951347e-1,
    0.419038818029165735901852432784e-1,
    -0.251390972343483509333252996350e-1];

    let mut phi = eval_polynomial(b, &COEFFS);

    if x < y {
        phi = 1.0 - phi;
    }

    let mut v = phi * r;
    let mut u = r - v;

    if d.z < 0.0 {
        mem::swap(&mut u, &mut v);
        u = 1.0 - u;
        v = 1.0 - v;
    }

    u = f32::copysign(u, d.x);
    v = f32::copysign(v, d.y);

    Vector2 {
        x: 0.5 * (u + 1.0),
        y: 0.5 * (v + 1.0)
    }
}

/// Wraps the given square for any `u, v` values that might fall outside of the regular `[-1, 1]^2` mapping.
pub fn wrap_equal_area_square(uv: Vector2f) -> Vector2f {
    let mut uv = uv;
    if uv.x < 0.0 {
        uv.x = -uv.x;
        uv.y = 1.0 - uv.y;
    } else if uv.x > 1.0 {
        uv.x = 2.0 - uv.x;
        uv.y = 1.0 - uv.y;
    }

    if uv.y < 0.0 {
        uv.x = 1.0 - uv.x;
        uv.y = -uv.y;
    } else if uv.y > 1.0 {
        uv.x = 1.0 - uv.x;
        uv.y = 2.0 - uv.y;
    }
    uv
}

#[derive(Clone, Copy)]
pub struct DirectionCone {
    pub w: Vector3f,
    pub cos_theta: f32
}

impl DirectionCone {

    pub fn default() -> Self {
        DirectionCone{w: Vector3f::default(), cos_theta: f32::INFINITY}
    }

    pub fn new(w: Vector3f, cos_theta: f32) -> Self {
        Self {w, cos_theta}
    }

    pub fn from_direction(w: Vector3f) -> Self {
        Self {w, cos_theta: 1.0}
    }

    pub fn empty(&self) -> bool {
        self.cos_theta == f32::INFINITY
    }

    pub fn all_directions() -> Self {
        Self::new(Vector3f::new(0.0, 0.0, 1.0), -1.0)
    }

    pub fn inside(d: &Self, w: Vector3f) -> bool {
        !d.empty() && Vector::dot(d.w, Vector::normalize(w)) >= d.cos_theta
    }

    /// Gives a DirectionCone which bounds the subtended angle by the given bounding box,
    /// specifically from the reference point of the given point.
    pub fn bound_directions(b: &Bounds3f, p: Vector3f) -> Self {
        let (center, radius) = b.bounding_sphere();
        if Vector::distance_squared(p, center) < square(radius) {
            return Self::all_directions();
        }
        let w = Vector::normalize(center - p);
        let sin2_theta_extent = square(radius) / Vector::distance_squared(center, p);
        let cos_theta_extent = safe_sqrt(1.0 - sin2_theta_extent);
        Self::new(w, cos_theta_extent)
    }

    pub fn union(a: &Self, b: &Self) -> Self {
        if a.empty() {
            return *b;
        } else if b.empty() {
            return *a;
        }

        let theta_a = arccos(a.cos_theta);
        let theta_b = arccos(b.cos_theta);
        let theta_d = Vector::angle_between(a.w, b.w);
        if f32::min(theta_d + theta_b, PI) <= theta_a {
            return *a;
        }
        if f32::min(theta_d + theta_a, PI) <= theta_b {
            return *b;
        }

        let theta_o = (theta_a + theta_d + theta_b) / 2.0;
        if theta_o >= PI {
            return Self::all_directions();
        }
        let theta_r = theta_o - theta_a;
        let w_r = Vector3::cross(a.w, b.w);
        if Vector::length_squared(w_r) == 0.0 {
            return Self::all_directions();
        }
        let w = Transform::rotate(degrees(theta_r), w_r).app(a.w);
        Self::new(w, f32::cos(theta_o))
    }
}

