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
use std::cell::Cell;
use std::f64;

pub struct Ray {
    pub o: Point3,
    pub d: Vector3,
    // We'll always need our Rays to be mutable so we can mutate this
    pub t_max: Cell<f64>, 
    pub time: f64,
    //TODO: implement Mediums (chapter 11.3)
    //pub medium: Option<dyn Medium>,
    pub diff: Option<Differential>
}

impl Ray {
    pub fn new(o: &Point3, d: &Vector3, t_max: f64, time: f64 /*, medium: Option<dyn Medium>*/) -> Ray {
        Ray {
            o: *o,
            d: *d,
            t_max: Cell::from(t_max),
            time,
            /*
            medium,
            */
            diff: None
        }
    }

    pub fn from(origin: &Point3, dest: &Vector3) -> Ray {
        Ray {
            o: *origin,
            d: *dest,
            t_max: Cell::from(f64::INFINITY),
            time: 0.,
            /*
            medium: None
            */
            diff: None
        }
    }

    pub fn at(&self, time: f64) -> Point3{
        (self.o.add_vec(&self.d)).mult(time)
    }

    pub fn scale_diffs(&mut self, by: f64) {
        if let Some(diff) = &mut self.diff {
            diff.rx_origin = diff.rx_origin.add_vec(&(diff.rx_origin - self.o)).mult(by);
            diff.ry_origin = diff.ry_origin.add_vec(&(diff.rx_origin - self.o)).mult(by);
            diff.rx_dir = self.d + (diff.rx_dir - self.d).mult(by);
            diff.ry_dir = self.d + (diff.ry_dir - self.d).mult(by);
        }
    }

    pub fn is_differential(&self) -> bool {
        self.diff.is_some()
    }
}

//TODO: make default ray
//TODO: add NAN checking functions

pub struct Differential {
    pub rx_origin: Point3,
    pub ry_origin: Point3,
    pub rx_dir: Vector3,
    pub ry_dir: Vector3
}