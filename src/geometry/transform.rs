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

use super::vectors::*;
use super::points::*;
use super::ray::*;
use super::normal::*;
use super::bounds::*;
use crate::math::gamma;

#[derive(PartialEq, Copy, Clone, Debug)]
pub struct Transform {
    m: Matrix4x4,
    m_inv: Matrix4x4
}

#[allow(dead_code)]
impl Transform {
    pub fn new() -> Transform {
        Transform {
            m: IDENTITY,
            m_inv: IDENTITY
        }
    }

    pub fn from_mat(m: [[f64; 4]; 4]) -> Transform {
        let mat = Matrix4x4::new(m);
        Transform {
            m: mat,
            m_inv: mat.inverse()
        }
    }

    pub fn from_parts(m: [[f64; 4]; 4], m_inv: [[f64; 4]; 4]) -> Transform {
        Transform {m: Matrix4x4::new(m), m_inv: Matrix4x4::new(m_inv)}
    }

    fn from_matrices(m: &Matrix4x4, m_inv: &Matrix4x4) -> Transform {
        Transform {m: *m, m_inv: *m_inv}
    }

    pub fn inverse(&self) -> Transform {
        Transform::from_matrices(&self.m_inv, &self.m)
    }

    pub fn transpose(&self) -> Transform {
        Transform::from_matrices(&self.m.transpose(), &self.m_inv.transpose())
    }

    pub fn is_identity(&self) -> bool {
        self.m == IDENTITY
    }
    
    pub fn trans_point(&self, p: &Point3) -> Point3 {
        let xp = self.m.m[0][0] * p.x + self.m.m[0][1] * p.y + self.m.m[0][2]  * p.z + self.m.m[0][3];
        let yp = self.m.m[1][0] * p.x + self.m.m[1][1] * p.y + self.m.m[1][2]  * p.z + self.m.m[1][3];
        let zp = self.m.m[2][0] * p.x + self.m.m[2][1] * p.y + self.m.m[2][2]  * p.z + self.m.m[2][3];
        let wp = self.m.m[3][0] * p.x + self.m.m[3][1] * p.y + self.m.m[3][2]  * p.z + self.m.m[3][3];
        assert!(wp != 0.); // if wp is zero then the transform will explode a little bit

        if wp == 1. { //this check allows us to avoid unnecessary fp division
            Point3::new(xp, yp, zp)
        }
        else {
            Point3::new(xp, yp, zp).div(wp)
        }

    }

    pub fn trans_vec(&self, v: &Vector3) -> Vector3 {
        Vector3 {
            x: self.m.m[0][0] * v.x + self.m.m[1][0] * v.y + self.m.m[2][0] * v.z,
            y: self.m.m[0][1] * v.x + self.m.m[1][1] * v.y + self.m.m[2][1] * v.z,
            z: self.m.m[0][2] * v.x + self.m.m[1][2] * v.y + self.m.m[2][2] * v.z
        }
    }

    /// Transforms a vector and returns the error of the transformation.
    pub fn trans_vec_with_err(&self, v: &Vector3) -> (Vector3, Vector3) {
        let v_err = Vector3 {
            x: gamma(3) * (self.m.m[0][0] * v.x) + (self.m.m[0][1] * v.y) + (self.m.m[0][2] * v.z),
            y: gamma(3) * (self.m.m[1][0] * v.x) + (self.m.m[1][1] * v.y) + (self.m.m[1][2] * v.z),
            z: gamma(3) * (self.m.m[2][0] * v.x) + (self.m.m[2][1] * v.y) + (self.m.m[2][2] * v.z),
        };

        let t = Vector3 {
            x: self.m.m[0][0] * v.x + self.m.m[0][1] * v.y + self.m.m[0][2] * v.z,
            y: self.m.m[1][0] * v.x + self.m.m[1][1] * v.y + self.m.m[1][2] * v.z,
            z: self.m.m[2][0] * v.x + self.m.m[2][1] * v.y + self.m.m[2][2] * v.z,
        };

        (t, v_err)
    }

    /// Performs a transformation on a vector, given the vector and error already existing in the vector.
    pub fn trans_vec_from_err(self, v: &Vector3, v_err: &Vector3) -> (Vector3, Vector3) {
        let abs_err = Vector3 {
            x: (gamma(3) + 1.) * (self.m.m[0][0].abs() * v_err.x + self.m.m[0][1] * v_err.y + self.m.m[0][1] * v_err.z) + 
                gamma(3) * self.m.m[0][0] * v.x + self.m.m[0][1].abs() * v.y + self.m.m[0][2].abs() * v.z,

            y: (gamma(3) + 1.) * (self.m.m[1][0].abs() * v_err.x + self.m.m[1][1] * v_err.y + self.m.m[1][1] * v_err.z) + 
            gamma(3) * self.m.m[0][0] * v.x + self.m.m[0][1].abs() * v.y + self.m.m[0][2].abs() * v.z,
            
            z: (gamma(3) + 1.) * (self.m.m[2][0].abs() * v_err.x + self.m.m[2][1] * v_err.y + self.m.m[2][1] * v_err.z) + 
            gamma(3) * self.m.m[2][0] * v.x + self.m.m[2][1].abs() * v.y + self.m.m[2][2].abs() * v.z,
        };

        let t = Vector3 {
            x: self.m.m[0][0] * v.x + self.m.m[0][1] * v.y + self.m.m[0][2] * v.z,
            y: self.m.m[1][0] * v.x + self.m.m[1][1] * v.y + self.m.m[1][2] * v.z,
            z: self.m.m[2][0] * v.x + self.m.m[2][1] * v.y + self.m.m[2][2] * v.z,
        };

        (t, abs_err)
    }

    pub fn trans_norm(&self, n: &Normal3) -> Normal3 {
        Normal3 {
            x: self.m.m[0][0] * n.x + self.m.m[1][0] * n.y + self.m.m[2][0] * n.z,
            y: self.m.m[0][1] * n.x + self.m.m[1][1] * n.y + self.m.m[2][1] * n.z,
            z: self.m.m[0][2] * n.x + self.m.m[1][2] * n.y + self.m.m[2][2] * n.z
        }
    }

    /// Performs a transformation on the given point and returns it and the computed error (in vec form).
    pub fn trans_pt_with_err(&self, p: &Point3) -> (Point3, Vector3) {
        // compute trans on point
        let xp = (self.m.m[0][0] * p.x + self.m.m[0][1] * p.y) + (self.m.m[0][2] * p.z + self.m.m[0][3]);
        let yp = (self.m.m[1][0] * p.x + self.m.m[1][1] * p.y) + (self.m.m[1][2] * p.z + self.m.m[1][3]);
        let zp = (self.m.m[2][0] * p.x + self.m.m[2][1] * p.y) + (self.m.m[2][2] * p.z + self.m.m[2][3]);
        let wp = (self.m.m[3][0] * p.x + self.m.m[3][1] * p.y) + (self.m.m[3][2] * p.z + self.m.m[3][3]);

        // compute err for point
        let x_sum = (self.m.m[0][0] * p.x).abs() + (self.m.m[0][1] * p.x).abs() + (self.m.m[0][2] * p.x).abs() + self.m.m[0][3].abs();
        let y_sum = (self.m.m[1][0] * p.x).abs() + (self.m.m[1][1] * p.x).abs() + (self.m.m[1][2] * p.x).abs() + self.m.m[1][3].abs();
        let z_sum = (self.m.m[2][0] * p.x).abs() + (self.m.m[2][1] * p.x).abs() + (self.m.m[2][2] * p.x).abs() + self.m.m[2][3].abs();

        assert!(wp != 0.);
        (if wp == 1. {Point3::new(xp, yp, zp)} 
        else {Point3::new(xp, yp, zp).div(wp)}, 
        Vector3::new(x_sum, y_sum, z_sum).mult(gamma(3)))
    }

    /// Performs a transformation on the given point, with inputted error and returns the transformation plus the absolute error.
    pub fn trans_pt_from_err(&self, p: &Point3, p_err: &Vector3) -> (Point3, Vector3) {
        // compute trans on point
        let xp = (self.m.m[0][0] * p.x + self.m.m[0][1] * p.y) + (self.m.m[0][2] * p.z + self.m.m[0][3]);
        let yp = (self.m.m[1][0] * p.x + self.m.m[1][1] * p.y) + (self.m.m[1][2] * p.z + self.m.m[1][3]);
        let zp = (self.m.m[2][0] * p.x + self.m.m[2][1] * p.y) + (self.m.m[2][2] * p.z + self.m.m[2][3]);
        let wp = (self.m.m[3][0] * p.x + self.m.m[3][1] * p.y) + (self.m.m[3][2] * p.z + self.m.m[3][3]);

        let abs_err = Vector3 {
            x: (gamma(3) + 1.) * ((self.m.m[0][0]).abs() * p_err.x + self.m.m[0][1].abs() * p_err.y + self.m.m[0][0] * p_err.z) + 
               (gamma(3) * (self.m.m[0][0].abs() * p.x + self.m.m[0][1] * p.y + self.m.m[0][2] * p.z + self.m.m[0][3].abs())),

            y: (gamma(3) + 1.) * ((self.m.m[1][0]).abs() * p_err.x + self.m.m[1][1].abs() * p_err.y + self.m.m[1][0] * p_err.z) + 
               (gamma(3) * (self.m.m[1][0].abs() * p.x + self.m.m[1][1] * p.y + self.m.m[1][2] * p.z + self.m.m[1][3].abs())),
            
            z: (gamma(3) + 1.) * ((self.m.m[2][0]).abs() * p_err.x + self.m.m[2][1].abs() * p_err.y + self.m.m[2][0] * p_err.z) + 
               (gamma(3) * (self.m.m[2][0].abs() * p.x + self.m.m[2][1] * p.y + self.m.m[2][2] * p.z + self.m.m[2][3].abs())),
        };

        assert!(wp != 0.); //TODO write a function for this that checks within a range
        
        (
            if wp == 1. {
                Point3::new(xp, yp, zp)
            } else {
                Point3::new(xp, yp, zp).div(wp)
            },
            abs_err
        )

    }

    pub fn trans_ray(&self, r: &Ray) -> Ray {
        // If the ray in question is differential, then we handle it here
        if r.is_differential() {
            return self.trans_ray_diff(r); //must explicitly return because the function continues
        }
        let (o, o_err) = self.trans_pt_with_err(&r.o);
        let d = self.trans_vec(&r.d);

        let len_sqr = d.length_sqr();
        let t_max = r.t_max.get(); //make a copy of the t_max because it's only like 64 bytes and we'll be copying it anyways

        if len_sqr > 0. {
            let dt = dot(&d.abs(), &o_err) / len_sqr;
        }

        Ray::new(&o, &d, t_max, r.time /* ,TODO implement mediums r.medium*/)
    }

    /// Performs a transformation on a ray with a guaranteed differential.
    /// This guarantee allows this function to straight up force unwrap the differential and perform operations on it.
    fn trans_ray_diff(&self, r: &Ray) -> Ray {

        //apply a transformation to everything except the diff element
        let mut t = self.trans_ray(&Ray::new(&r.o, &r.d, r.t_max.get(), r.time));

        //TODO make sure we actually need the clone call
        let rd = r.diff.unwrap();

        let d = Differential { 
            rx_origin: self.trans_point(&rd.rx_origin),
            ry_origin: self.trans_point(&rd.ry_origin),
            rx_dir: self.trans_vec(&rd.rx_dir),
            ry_dir: self.trans_vec(&rd.ry_dir)
        };

        t.diff = Some(d);
        t
    }

    /// Computes the transformation on the given ray, computing the error for
    /// the origin transform and vector transforms, respectively.
    pub fn trans_ray_with_err(&self, r: &Ray) -> (Ray, Vector3, Vector3) {
        let (mut o, o_err) = self.trans_pt_with_err(&r.o);
        let (d, d_err) = self.trans_vec_with_err(&r.d);
        let t_max = r.t_max.get();
        let len_sqr = d.length_sqr();

        if len_sqr > 0. {
            let dt = dot(&d.abs(), &o_err) / len_sqr;
            o.add_vec_mut(&d.mult(dt));
            //t_max -= dt;
        }

        (Ray::new(&o, &d, t_max, r.time/*, r.medium*/), o_err, d_err)
    }

    /// Computes the transformation on the given ray, computing the error given already 
    /// existing error components.
    pub fn trans_ray_from_err
    (&self, r: &Ray, d_err: &Vector3, o_err: &Vector3)
    -> (Ray, Vector3, Vector3) {
        let (mut o, o_err) = self.trans_pt_from_err(&r.o, &o_err);
        let (d, d_err) = self.trans_vec_from_err(&r.d, &d_err);
        let t_max = r.t_max.get();
        let len_sqr = d.length_sqr();
        if len_sqr > 0. {
            let dt = dot(&d.abs(), &o_err) / len_sqr;
            o.add_vec_mut(&d.mult(dt));
        }

        (Ray::new(&o, &d, t_max, r.time/*, r.medium*/), o_err, d_err)
    }

    /*
    pub fn trans_bounds(&self, b: &Bounds3) -> Bounds3 {

    }
    */

    pub fn translate(delta: &Vector3) -> Transform {
        let m = [
                [1., 0., 0., delta.x],
                [0., 1., 0., delta.y], 
                [0., 0., 1., delta.z],
                [0., 0., 0., 1.]
                ];
        let m_inv = [
                    [1., 0., 0., -delta.x],
                    [0., 1., 0., -delta.y], 
                    [0., 0., 1., -delta.z],
                    [0., 0., 0., 1.]
                    ];
        Transform::from_parts(m, m_inv)
    }

    pub fn scale(x: f64, y: f64, z: f64) -> Transform {
        let m = [
                [x, 0., 0., 0.],
                [0., y, 0., 0.], 
                [0., 0., z, 0.],
                [0., 0., 0., 1.]
                ];
        let m_inv = [
                    [1./x, 0., 0., 0.],
                    [0., 1./y, 0., 0.], 
                    [0., 0., 1./z, 0.],
                    [0., 0., 0., 1.]
                    ];
        Transform::from_parts(m, m_inv)
    }

    /// Determines whether or not a given transformation has a scaling term in it.
    pub fn has_scale(&self) -> bool {
        let a = self.trans_vec(&Vector3::new(1., 0., 0.)).length_sqr(); //transform the standard basis vectors to determine.
        let b = self.trans_vec(&Vector3::new(0., 1., 0.)).length_sqr();
        let c = self.trans_vec(&Vector3::new(0., 0., 1.)).length_sqr();
        not_one(a) || not_one(b) || not_one(c)
    }

}

fn not_one(x: f64) -> bool {
    x < 0.999 && x > 1.001
}

#[derive(Copy, Clone, PartialEq, Debug)]
/// Represents a 4x4 matrix and some of the useful associated operations that can be performed on one.
/// 
/// TODO decide whether or not this gets to be public
struct Matrix4x4 {
    pub m: [[f64; 4]; 4],
}

///The identity matrix.
#[allow(dead_code)]
const IDENTITY: Matrix4x4 = Matrix4x4 {
    m: [[1., 0., 0., 0.],
        [0., 1., 0., 0.],
        [0., 0., 1., 0.],
        [0., 0., 0., 1.]]
};

#[allow(dead_code)]
impl Matrix4x4 {
    pub fn new(m: [[f64; 4]; 4]) -> Matrix4x4 {
        Matrix4x4 {
            m,
        }
    }

    pub fn empty() -> Matrix4x4 {
        Matrix4x4 {
            m: [[0., 0., 0., 0.],
                [0., 0., 0., 0.],
                [0., 0., 0., 0.],
                [0., 0., 0., 0.]]
        }
    }

    pub fn transpose(&self) -> Matrix4x4 {
        Matrix4x4 {
            m: [[self.m[0][0], self.m[1][0], self.m[2][0], self.m[3][0]],
                [self.m[0][1], self.m[1][1], self.m[2][1], self.m[3][1]],
                [self.m[0][2], self.m[1][2], self.m[2][2], self.m[3][2]],
                [self.m[0][3], self.m[1][3], self.m[2][3], self.m[3][3]]]
        }
    }

    /// Computes the inverse of the matrix. Currently uses stable Gauss-Jordan elimination, via the method used in `pbrt`.
    /// 
    /// TODO: reimplement this in a better manner
    /// 
    /// TODO: fix this function
    /// #### Panics 
    /// * if the matrix is not invertible
    /// * if the matrix is singular (these are equivalent)
    pub fn inverse(&self) -> Matrix4x4 {
        let mut indxc: [usize; 4] = [0; 4];
        let mut indxr: [usize; 4] = [0; 4];

        let mut ipiv: [usize; 4] = [0, 0, 0, 0];
        let mut minv = [[0.; 4]; 4];

        minv.copy_from_slice(&self.m);

        for i in 0..4 {
            let mut i_row = 0;
            let mut i_col = 0;

            let mut big = 0.;
            // Choose pivot
            for j in 0..4 {
                if ipiv[j] != 1 {
                    for k in 0..4 {
                        if ipiv[k] == 0 {
                            if (minv[j][k]).abs() >= big {
                                big = minv[j][k];
                                i_row = j;
                                i_col = k;
                            }
                        } else if ipiv[k] > 1 {
                            panic!("Singular matrix found while trying to invert");
                        }
                    }
                }
            }

            ipiv[i_col] += 1;
            // Swap rows _irow_ and _icol_ for pivot
            if i_row != i_col {
                for k in 0..4 {
                    // swap minv at [i_row] [k] and [i_col] [k]
                    let swap = minv[i_row][k];
                    minv[i_row][k] = minv[i_col][k];
                    minv[i_col][k] = swap;
                }
            }
            indxr[i] = i_row;
            indxc[i] = i_col;
            if minv[i_col][i_col] == 0. {
                panic!("Singular matrix found while trying to invert");
            }

            // Put this row's leading term in RREF
            let pivinv = 1. / minv[i_col][i_col];
            minv[i_col][i_col] = 1.;
            for j in 0..4 {
                minv[i_col][j] *= pivinv;
            } 

            // Take current row and subtract from other rows to create zeroes
            for j in 0..4 {
                if j != i_col {
                    let save = minv[j][i_col];
                    minv[j][i_col] = 0.;
                    for k in 0..4 {
                        minv[j][k] -= minv[i_col][k] * save;
                    } 
                }
            }
        }
        // Swap columns to reflect permutation
        for j in (0..3).rev() {
            if indxr[j] != indxc[j] {
                for k in 0..4 {
                    // swap minv[k][indxr[j]] and minv[k][indxc[j]]
                    let swap = minv[k][indxr[j]];
                    minv[j][indxr[j]] = minv[k][indxc[j]];
                    minv[k][indxc[j]] = swap;
                }
            }
        }
        Matrix4x4::new(minv)
    }

    /* TODO implement?
    pub fn write_to_file
    */
}

impl std::ops::Mul for Matrix4x4 {
    type Output = Self;
    fn mul(self, other: Matrix4x4) -> Matrix4x4 {
        let mut out = Matrix4x4::empty();
        for i in 0..4 {
            for j in 0..4 {
                out.m[i][j] = 
                    self.m[i][0] * other.m[0][j] + 
                    self.m[i][1] * other.m[1][j] + 
                    self.m[i][2] * other.m[2][j] + 
                    self.m[i][3] * other.m[3][j];
            }
        }
        out
    }
}

#[cfg(test)]
#[allow(dead_code)]
mod test {
    use super::*;
    
    #[test]
    fn test_inv_identity() {
        assert_eq!(IDENTITY, IDENTITY.inverse());
    }

    //#[test] This doesn't work, TODO figure out why
    fn test_inv_complicated() { // Got this one by calculating it by hand.
        let a = Matrix4x4 {
            m: [[1., 2., 2., 1.], 
                [1., 1., 1., 2.], 
                [-1., 1., -1., 4.], 
                [0., 1., 1., 1.]]
        }.inverse();

        let a_inv = Matrix4x4 {
            m: [[0.5, 0.5, 0., -1.5], 
                [1.5, -1., 0.5, -1.5], 
                [-1., 0.5, -0.5, 2.], 
                [-0.5, 0.5, 0., 0.5]]
        };

        let margin = 0.0000001;

        let mut inv_correct = true;
        for i in 0..4 {
            for j in 0..4 {
                if a.m[i][j] - a_inv.m[i][j] > margin {
                    inv_correct = false;
                }
            }
        }
        assert!(inv_correct);
    }
    
}