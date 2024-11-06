use std::ops::{Add, AddAssign, Index, IndexMut, Mul, Div};
use crate::math::Num;

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

    // TODO implement compensated inner product
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
        let inv_det = 1.0 / det;
        let mut r = Self::zero();



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

impl<const N: usize> Div<f32> for SquareMatrix<N> {
    type Output = Self;

    fn div(self, rhs: f32) -> Self::Output {
        let recip = 1.0 / rhs;
        self * recip
    }
}

impl<const N: usize> SquareMatrix<N> {

}



pub struct Transform {
    m: SquareMatrix<4>
}