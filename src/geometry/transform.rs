use std::ops::{Add, AddAssign, Index, IndexMut, Mul, Div};
use crate::math::{inner_prod, Num, square};
use crate::geometry::vector::Vector3f;

use super::vector::{Vector, Vector3};

#[derive(PartialEq, Clone, Copy)]
pub struct SquareMatrix<const N: usize> {
    m: [[f32; N]; N],
}

impl<const N: usize> SquareMatrix<N> {
    pub fn zero() -> Self {
        Self::new([[0.0; N]; N])
    }

    pub fn ident() -> Self {
        let mut mat = Self::zero();
        for i in 0..N {
            for j in 0..N {
                if i == j {mat.m[i][j] = 1.0;}
            }
        }
        mat
    }

    pub fn new(m: [[f32; N]; N]) -> Self {
        Self {
            m
        }
    }

    pub fn diagonal(d: [f32; N]) -> Self {
        let mut mat = Self::zero();
        for i in 0..N {
            for j in 0..N {
                if i == j {mat.m[i][j] = d[i];}
            }
        }
        mat
    }

    pub fn is_ident(&self) -> bool {
        for i in 0..N {
            for j in 0..N {
                if i == j && self.m[i][j] != 1.0 {
                    return false;
                } else if self.m[i][j] != 0.0 {
                    return false;
                }
            }
        }
        true
    }
    
    //TODO: this might be way too many trait bounds for my purposes, and I might need to cool off. Unsure.
    pub fn mul<Result, Input, ResultIndex, InputIndex>(m: &Self, v: &Input) -> Result where 
        ResultIndex: From<f32> + AddAssign<f32>,
        InputIndex: Mul<f32, Output=f32> + Copy,
        Result: IndexMut<usize, Output=ResultIndex> + Default,
        Input: Index<usize, Output=InputIndex>,
        f32: Mul<InputIndex, Output=f32>
        {
        let mut res = Result::default();
        for i in 0..N {
            res[i] = ResultIndex::from(0.0);
            for j in 0..N {
                res[i] += m[i][j] * v[j];
            }
        }
        res
    }

    pub fn determinant(m: &Self) -> f32 {
        match N {
            1 => m[0][0],
            2 => Self::det2(m),
            3 => Self::det3(m),
            4 => Self::det4(m),
            _ => panic!("matrix of dimension > 4 attempted to take determinant, currently unsupported")
        }
    }

    fn det2(m: &Self) -> f32 {
        f32::diff_products(m[0][0], m[1][1], m[0][1], m[1][0])
    }
    
    fn det3(m: &Self) -> f32 {
        let m12 = f32::diff_products(m[1][1], m[2][2], m[1][2], m[2][1]);
        let m02 = f32::diff_products(m[1][0], m[2][2], m[1][2], m[2][0]);
        let m01 = f32::diff_products(m[1][0], m[2][1], m[1][1], m[2][0]);
        f32::fma(m[0][2], m01, f32::diff_products(m[0][0], m12, m[0][1], m02))
    }
    
    fn det4(m: &Self) -> f32 {
        let s0 = f32::diff_products(m[0][0], m[1][1], m[1][0], m[0][1]);
        let s1 = f32::diff_products(m[0][0], m[1][2], m[1][0], m[0][2]);
        let s2 = f32::diff_products(m[0][0], m[1][3], m[1][0], m[0][3]);
        let s3 = f32::diff_products(m[0][1], m[1][2], m[1][1], m[0][2]);
        let s4 = f32::diff_products(m[0][1], m[1][3], m[1][1], m[0][3]);
        let s5 = f32::diff_products(m[0][2], m[1][3], m[1][2], m[0][3]);

        let c0 = f32::diff_products(m[2][0], m[3][1], m[3][0], m[2][1]);
        let c1 = f32::diff_products(m[2][0], m[3][2], m[3][0], m[2][2]);
        let c2 = f32::diff_products(m[2][0], m[3][3], m[3][0], m[2][3]);
        let c3 = f32::diff_products(m[2][1], m[3][2], m[3][1], m[2][2]);
        let c4 = f32::diff_products(m[2][1], m[3][3], m[3][1], m[2][3]);
        let c5 = f32::diff_products(m[2][2], m[3][3], m[3][2], m[2][3]);

        f32::diff_products(s0, c5, s1, c4) + f32::diff_products(s2, c3, -s3, c2) +
            f32::diff_products(s5, c0, s4, c1)
    }

    pub fn inverse(m: &Self) -> Option<Self> {
        match N {
            1 => Some(*m),
            2 => Self::inv2(m),
            3 => Self::inv3(m),
            4 => Self::inv4(m),
            _ => panic!("matrix of dimensions > 4 attempted to take inverse, currently unsupported")
        }
    }

    fn inv2(m: &Self) -> Option<Self> {
        let det = Self::determinant(m);
        if det == 0.0 {
            return None;
        }
        let inv_det = 1.0 / det;
        let mut r = Self::zero();
        r[0][0] = inv_det * m[1][1];
        r[0][1] = inv_det * -m[0][1];
        r[1][0] = inv_det * -m[1][0];
        r[1][1] = inv_det * m[0][0];
        Some(r)
    }

    fn inv3(m: &Self) -> Option<Self> {
        let det = Self::determinant(m);
        if det == 0.0 {
            return None;
        }
        let inv_det = 1.0 / det;
        let mut r = Self::zero();

        r[0][0] = inv_det * f32::diff_products(m[1][1], m[2][2], m[1][2], m[2][1]);
        r[1][0] = inv_det * f32::diff_products(m[1][2], m[2][0], m[1][0], m[2][2]);
        r[2][0] = inv_det * f32::diff_products(m[1][0], m[2][1], m[1][1], m[2][0]);
        r[0][1] = inv_det * f32::diff_products(m[0][2], m[2][1], m[0][1], m[2][2]);
        r[1][1] = inv_det * f32::diff_products(m[0][0], m[2][2], m[0][2], m[2][0]);
        r[2][1] = inv_det * f32::diff_products(m[0][1], m[2][0], m[0][0], m[2][1]);
        r[0][2] = inv_det * f32::diff_products(m[0][1], m[1][2], m[0][2], m[1][1]);
        r[1][2] = inv_det * f32::diff_products(m[0][2], m[1][0], m[0][0], m[1][2]);
        r[2][2] = inv_det * f32::diff_products(m[0][0], m[1][1], m[0][1], m[1][0]);
        Some(r)
    }

    fn inv4(m: &Self) -> Option<Self> {
        let s0 = f32::diff_products(m[0][0], m[1][1], m[1][0], m[0][1]);
        let s1 = f32::diff_products(m[0][0], m[1][2], m[1][0], m[0][2]);
        let s2 = f32::diff_products(m[0][0], m[1][3], m[1][0], m[0][3]);
        let s3 = f32::diff_products(m[0][1], m[1][2], m[1][1], m[0][2]);
        let s4 = f32::diff_products(m[0][1], m[1][3], m[1][1], m[0][3]);
        let s5 = f32::diff_products(m[0][2], m[1][3], m[1][2], m[0][3]);

        let c0 = f32::diff_products(m[2][0], m[3][1], m[3][0], m[2][1]);
        let c1 = f32::diff_products(m[2][0], m[3][2], m[3][0], m[2][2]);
        let c2 = f32::diff_products(m[2][0], m[3][3], m[3][0], m[2][3]);
        let c3 = f32::diff_products(m[2][1], m[3][2], m[3][1], m[2][2]);
        let c4 = f32::diff_products(m[2][1], m[3][3], m[3][1], m[2][3]);
        let c5 = f32::diff_products(m[2][2], m[3][3], m[3][2], m[2][3]);


        let det = Self::determinant(m);
        if det == 0.0 {
            return None;
        }
        let s = 1.0 / det;
        let mut r = Self::zero();

        r[0][0] = s * inner_prod(&[m[1][1], c5, m[1][3], c3, -m[1][2], c4]);
        r[0][1] = s * inner_prod(&[-m[0][1], c5, m[0][2], c4, -m[0][3], c3]);
        r[0][2] = s * inner_prod(&[m[3][1], s5, m[3][3], s3, -m[3][2], s4]);
        r[0][3] = s * inner_prod(&[-m[2][1], s5, m[2][2], s4, -m[2][3], s3]);
        
        r[1][0] = s * inner_prod(&[-m[1][0], c5, m[1][2], c2, -m[1][3], c1]);
        r[1][1] = s * inner_prod(&[m[0][0], c5, m[0][3], c1, -m[0][2], c2]);
        r[1][2] = s * inner_prod(&[-m[3][0], s5, m[3][2], s2, -m[3][3], s1]);
        r[1][3] = s * inner_prod(&[m[2][0], s5, m[2][3], s1, -m[2][2], s2]);
        
        r[2][0] = s * inner_prod(&[m[1][0], c4, m[1][3], c0, -m[1][1], c2]);
        r[2][1] = s * inner_prod(&[-m[0][0], c4, m[0][1], c2, -m[0][3], c0]);
        r[2][2] = s * inner_prod(&[m[3][0], s4, m[3][3], s0, -m[3][1], s2]);
        r[2][3] = s * inner_prod(&[-m[2][0], s4, m[2][1], s2, -m[2][3], s0]);

        r[3][0] = s * inner_prod(&[-m[1][0], c3, m[1][1], c1, -m[1][2], c0]);
        r[3][1] = s * inner_prod(&[m[0][0], c3, m[0][2], c0, -m[0][1], c1]);
        r[3][2] = s * inner_prod(&[-m[3][0], s3, m[3][1], s1, -m[3][2], s0]);
        r[3][3] = s * inner_prod(&[m[2][0], s3, m[2][2], s0, -m[2][1], s1]);
        
        Some(r)
    }

    pub fn transpose(m: &Self) -> Self {
        let mut r = Self::zero();
        for i in 0..N {
            for j in 0..N {
                r[i][j] = m[j][i];
            }
        }
        r
    }

    pub fn linear_least_squares(a: &[[f32; N]], b: &[[f32; N]], rows: usize) -> Option<Self> {
        let mut ata = Self::zero();
        let mut atb = Self::zero();

        for i in 0..N {
            for j in 0..N {
                for r in 0..rows {
                    ata[i][j] += a[r][i] * a[r][j];
                    atb[i][j] += a[r][i] * b[r][j];
                }
            }
        }

        let ata_i = Self::inverse(&ata);
        match ata_i {
            None => None,
            Some(atai) => Some(Self::transpose(&(atai * &atb)))
        }
    }

    fn mul3(m1: &Self, m2: &Self) -> Self {
        let mut r = Self::zero();
        for i in 0..3 {
            for j in 0..3 {
                r[i][j] = inner_prod(&[m1[i][0], m2[0][j], m1[1][j], m2[1][j], m1[i][2], m2[2][j]]);
            }
        }
        r
    }

    fn mul4(m1: &Self, m2: &Self) -> Self {
        let mut r = Self::zero();
        for i in 0..4 {
            for j in 0..4 {
                r[i][j] = inner_prod(&[m1[i][0], m2[0][j], m1[i][1], m2[1][j], m1[i][2], m2[2][j], m1[i][3], m2[3][j]]);
            }
        }
        r
    }

    fn muln(m1: &Self, m2: &Self) -> Self {
        let mut r = Self::zero();
        for i in 0..N {
            for j in 0..N {
                for k in 0..N {
                    r[i][j] = f32::fma(m1[i][k], m2[k][j], r[i][j]);
                }
            }
        }
        r
    }

    pub fn nan() -> Self {
        let mut r = Self::zero();
        for i in 0..N {
            for j in 0..N {
                r[i][j] = f32::NAN;
            }
        }
        r
    }
}

impl<const N: usize> Index<usize> for SquareMatrix<N> {
    type Output = [f32; N];

    fn index<'a>(&'a self, index: usize) -> &'a Self::Output {
        &self.m[index]
    }
}

impl<const N: usize> IndexMut<usize> for SquareMatrix<N> {
    fn index_mut<'a>(&'a mut self, index: usize) -> &'a mut Self::Output {
        &mut self.m[index]
    }
}

impl<const N: usize> Eq for SquareMatrix<N> { }

impl<const N: usize> Add<&SquareMatrix<N>> for SquareMatrix<N> {
    type Output = Self;
    
    fn add(self, rhs: &Self) -> Self::Output {
        let mut r = self;
        for i in 0..N {
            for j in 0..N {
                r.m[i][j] += rhs.m[i][j]
            }
        }
        r
    }
}

impl<const N: usize> Mul<f32> for SquareMatrix<N> {
    type Output = Self;

    fn mul(self, rhs: f32) -> Self::Output {
        let mut r = self;
        for i in 0..N {
            for j in 0..N {
                r.m[i][j] *= rhs;
            }
        }
        r
    }
}

impl<const N: usize> Mul<&Self> for SquareMatrix<N> {
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self::Output {
        match N {
            3 => Self::mul3(&self, rhs),
            4 => Self::mul4(&self, rhs),
            _ => Self::muln(&self, rhs)
        }
    }
}

impl<const N: usize> Mul<Vector3f> for SquareMatrix<N> {
    type Output = Vector3f;
    
    fn mul(self, rhs: Vector3f) -> Self::Output {
        let mut res = Vector3f::new(0.0, 0.0, 0.0);
        for i in 0..N {
            res[i] = 0.0;
            for j in 0..N {
                res[i] += self.m[i][j] * rhs[j];
            }
        }
        res
    }
}

impl<const N: usize> Div<f32> for SquareMatrix<N> {
    type Output = Self;

    fn div(self, rhs: f32) -> Self::Output {
        let recip = 1.0 / rhs;
        self * recip
    }
}

impl<const N: usize> SquareMatrix<N> {

}

#[derive(PartialEq)]
pub struct Transform {
    pub m: SquareMatrix<4>,
    pub m_inv: SquareMatrix<4>
}

impl Transform {
    pub fn default() -> Self {
        Self {
            m: SquareMatrix::ident(),
            m_inv: SquareMatrix::ident()
        }
    }

    pub fn new(m: &SquareMatrix<4>) -> Self {
        match SquareMatrix::inverse(m) {
            Some(m_inv) => Self {m: m.clone(), m_inv},
            None => Self {m: m.clone(), m_inv: SquareMatrix::nan()}
        }
    }

    pub fn from_array(m: [[f32;4];4]) -> Self {
        Self::new(&SquareMatrix::new(m))
    }

    pub fn with_inv(m: &SquareMatrix<4>, m_inv: &SquareMatrix<4>) -> Self {
        Self {
            m: *m,
            m_inv: *m_inv
        }
    }

    pub fn invert(t: &Self) -> Self {
        Self::with_inv(&t.m_inv, &t.m)
    }

    pub fn transpose(t: &Self) -> Self {
        Self::with_inv(&SquareMatrix::transpose(&t.m), &SquareMatrix::transpose(&t.m_inv))
    }

    pub fn is_ident(t: &Self) -> bool {
        t.m.is_ident()
    }

    pub fn translate(by: Vector3f) -> Self {
        Self {
            m: SquareMatrix::new(
                [[1.0, 0.0, 0.0, by.x],
                [0.0, 1.0, 0.0, by.y],
                [0.0, 0.0, 1.0, by.z],
                [0.0, 0.0, 0.0, 1.0]]),
            m_inv: SquareMatrix::new(
                [[1.0, 0.0, 0.0, by.x],
                [0.0, 1.0, 0.0, by.y],
                [0.0, 0.0, 1.0, by.z],
                [0.0, 0.0, 0.0, 1.0]])
        }
    }

    pub fn scale(by: Vector3f) -> Self {
        Self {
            m: SquareMatrix::new(
                [[1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0]]),
            m_inv: SquareMatrix::new(
                [[1.0, 0.0, 0.0, by.x],
                [0.0, 1.0, 0.0, by.y],
                [0.0, 0.0, 1.0, by.z],
                [0.0, 0.0, 0.0, 1.0]])
        }
    }

    pub fn has_scaling(&self, tolerance: Option<f32>) -> bool {
        const DEFAULT_TOL: f32 = 1.0e-3;

        let la = Vector::length_squared(self.m * Vector3f::new(1.0, 0.0, 0.0));
        let lb = Vector::length_squared(self.m * Vector3f::new(0.0, 1.0, 0.0));
        let lc = Vector::length_squared(self.m * Vector3f::new(0.0, 0.0, 1.0));
        let tol = match tolerance {
            Some(t) => t,
            None => DEFAULT_TOL
        };
        f32::abs(la - 1.0) > tol && f32::abs(lb - 1.0) > tol && f32::abs(lc - 1.0) > tol
    }

    pub fn rotate_x(theta: f32) -> Self {
        let sin_theta = f32::sin(theta);
        let cos_theta = f32::cos(theta);
        let m = SquareMatrix::new(
            [[1.0, 0.0, 0.0, 0.0],
             [0.0, cos_theta, -sin_theta, 0.0],
             [0.0, sin_theta, cos_theta, 0.0],
             [0.0, 0.0, 0.0, 1.0]]);
        Self::with_inv(&m, &SquareMatrix::transpose(&m))
    }

    pub fn rotate_y(theta: f32) -> Self {
        let sin_theta = f32::sin(theta);
        let cos_theta = f32::cos(theta);
        let m = SquareMatrix::new(
            [[cos_theta, 0.0, sin_theta, 0.0],
             [0.0, 1.0, 0.0, 0.0],
             [-sin_theta, 0.0, cos_theta, 0.0],
             [0.0, 0.0, 0.0, 1.0]]);
        Self::with_inv(&m, &SquareMatrix::transpose(&m))
    }

    pub fn rotate_z(theta: f32) -> Self {
        let sin_theta = f32::sin(theta);
        let cos_theta = f32::cos(theta);
        let m = SquareMatrix::new(
            [[cos_theta, -sin_theta, 0.0, 0.0],
             [sin_theta, cos_theta, 0.0, 0.0],
             [0.0, 0.0, 0.0, 0.0],
             [0.0, 0.0, 0.0, 1.0]]);
        Self::with_inv(&m, &SquareMatrix::transpose(&m))
    }

    /// Provides a rotation matrix which rotates about the given axis by the specified angle,
    /// with the sin and cos of that angle already precalculated.
    pub fn rotate_precalc(sin_theta: f32, cos_theta: f32, axis: Vector3f) -> Self {
        let a = Vector::normalize(axis);
        let mut m = SquareMatrix4::ident();
        m[0][0] = square(a.x) + (1.0 - square(a.x)) * cos_theta;
        m[0][1] = square(a.y) * (1.0 - cos_theta) - a.z * sin_theta;
        m[0][2] = a.x * a.z * (1.0 - cos_theta) + a.y * sin_theta;
        m[0][3] = 0.0;

        m[1][0] = a.x * a.y * (1.0 - cos_theta) + a.z * sin_theta;
        m[1][1] = square(a.y) + (1.0 - square(a.y)) * cos_theta;
        m[1][2] = a.y * a.z * (1.0 - cos_theta) - a.x * sin_theta;
        m[1][3] = 0.0;

        m[2][0] = a.x * a.z * (1.0 - cos_theta) - a.y * sin_theta;
        m[2][1] = a.y * a.z * (1.0 - cos_theta) + a.x * sin_theta;
        m[2][2] = square(a.z) + (1.0 - square(a.z)) * cos_theta;
        m[2][3] = 0.0;

        Self::with_inv(&m, &SquareMatrix::transpose(&m))
    }

    pub fn rotate(theta: f32, axis: Vector3f) -> Self {
        Self::rotate_precalc(f32::sin(theta), f32::cos(theta), axis)
    }

    pub fn rotate_from_to(from: Vector3f, to: Vector3f) -> Self {
        const GT_SQRT2_ON_2: f32 = 0.72;

        let refl = if f32::abs(from.x) < GT_SQRT2_ON_2 && f32::abs(to.x) < GT_SQRT2_ON_2 {
            Vector3f::new(1.0, 0.0, 0.0)
        } else if f32::abs(from.y) < GT_SQRT2_ON_2 && f32::abs(to.y) < GT_SQRT2_ON_2 {
            Vector3f::new(0.0, 1.0, 0.0)
        } else {
            Vector3f::new(0.0, 0.0, 1.0)
        };
        let u = refl - from;
        let v = refl - to;
        let mut r = SquareMatrix4::ident();
        for i in 0..3 {
            for j in 0..3 {
                r[i][j] = if i == j {1.0} else {0.0} -
                    2.0 / Vector::dot(u, u) * u[i] * u[j] - 
                    2.0 / Vector::dot(v, v) * v[i] * v[j] +
                    4.0 * Vector::dot(u, v) / Vector::dot(u, u) * Vector::dot(v, v) * v[i] * u[j];
            }
        }
        Self::with_inv(&r, &SquareMatrix::transpose(&r))
    }

    pub fn look_at(pos: Vector3f, look: Vector3f, up: Vector3f) -> Self {
        let dir = Vector::normalize(look - pos);
        let right = Vector::normalize(Vector3::cross(Vector::normalize(up), dir));
        let new_up = Vector3::cross(dir, right);
        let m_inv = SquareMatrix::new(
            [[right.x, new_up.x, dir.x, pos.x],
             [right.y, new_up.y, dir.y, pos.y],
             [right.z, new_up.z, dir.z, pos.z],
             [0.0, 0.0, 0.0, 1.0]]
        )
        let m = SquareMatrix::inverse(&m_inv).unwrap();

        Self::with_inv(&m, &m_inv)
    }

}

impl Eq for Transform {}

impl Mul<Vector3f> for Transform {
    type Output = Vector3f;
    
    fn mul(self, rhs: Vector3f) -> Self::Output {
        unimplemented!()
    }
}

type SquareMatrix4 = SquareMatrix<4>;
