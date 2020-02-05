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

/*
The pbrt book defines SurfaceInteraction and MediumInteraction objects as subclasses of Interaction.
Unless it becomes obvious that this approach is counterproductive, I'm going to make Interaction a
field of the SurfaceInteraction and MediumInteraction structs.
*/
use super::points::*;
use super::vectors::*;
use super::normal::*;
use super::ray::*;
use crate::media::mediuminterface::*;

/// Defines the common elements of an interaction. `dyn Inter` allows use of a
/// generic interaction.
pub trait Inter {
    fn is_surface_interaction(&self) -> bool;
    fn is_medium_interaction(&self) -> bool;
    fn interaction(&self) -> Interaction;
}

/// The common elements of a surface interaction.
#[derive(Copy, Clone, PartialEq)]
pub struct Interaction {
    pub p: Point3,
    pub time: f64,
    pub p_err: Vector3,
    /// the negative ray direction, w_0, which is the outgoing direction of light.
    pub wo: Vector3,
    pub n: Normal3,
    pub medium_interface: MediumInterface
}

impl Inter for Interaction {
    fn is_surface_interaction(&self) -> bool {
        false
    }

    fn is_medium_interaction(&self) -> bool {
        false
    }

    fn interaction(&self) -> Interaction {
        *self
    }
}

impl Interaction {

    /// Initializes an Interaction with all default values.
    pub fn default() -> Interaction {
        Interaction::new(Point3::default(), 0., ZERO_3, ZERO_3, Normal3::default(), MediumInterface::default())
    }

    /// Initializes an Interaction.
    pub fn new(p: Point3, time: f64, p_err: Vector3, wo: Vector3, n: Normal3, medium_interface: MediumInterface) -> Interaction {
        Interaction {
            p,
            time,
            p_err,
            wo,
            n,
            medium_interface
        }
    }
}

