
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

pub const ONE_MINUS_EPSILON: f32 = 1.0 - f32::EPSILON;

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
        return self;
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
