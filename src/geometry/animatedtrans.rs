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

use super::transform::*;
use super::vectors::*;
use super::quaternion::*;
use super::points::*;
use super::bounds::*;
use super::ray::*;
use crate::math::*;
use std::ops::{Add, Sub, Mul};
use std::f64;
use crate::*;

/// Contains all of the needed information to compute the animated transformation, per `pbrt`.
pub struct AnimatedTransform {
    animated: bool, //is this transform actually an animation at all?
    rotating: bool, //does this transform have any rotation?
    start_trans: Transform, //TODO figure out if these can be owned or not
    end_trans: Transform,
    start_time: f64,
    end_time: f64,
    t: [Vector3; 2], //the start and end locations
    r: [Quaternion; 2], //the start and end rotations
    s: [Matrix4x4; 2], //the start and end scales
    c1: [DerivativeTerm; 3],
    c2: [DerivativeTerm; 3],
    c3: [DerivativeTerm; 3],
    c4: [DerivativeTerm; 3],
    c5: [DerivativeTerm; 3],
}

impl AnimatedTransform {
    pub fn new(start_trans: Transform, end_trans: Transform, start_time: f64, end_time: f64) -> AnimatedTransform {
        if start_trans == end_trans {
            return AnimatedTransform {
                animated: false,
                rotating: false,
                start_trans,
                end_trans,
                start_time,
                end_time,
                t: [Vector3::new(0., 0., 0.), Vector3::new(0., 0., 0.)],
                r: [Quaternion::default(), Quaternion::default()],
                s: [IDENTITY, IDENTITY],
                c1: [DerivativeTerm::default(), DerivativeTerm::default(), DerivativeTerm::default()],
                c2: [DerivativeTerm::default(), DerivativeTerm::default(), DerivativeTerm::default()],
                c3: [DerivativeTerm::default(), DerivativeTerm::default(), DerivativeTerm::default()],
                c4: [DerivativeTerm::default(), DerivativeTerm::default(), DerivativeTerm::default()],
                c5: [DerivativeTerm::default(), DerivativeTerm::default(), DerivativeTerm::default()]
            };
        }
        let (mut t0, mut r0, mut s0) = Self::decompose(&start_trans.m);
        let (mut t1, mut r1, mut s1) = Self::decompose(&end_trans.m);
        if r0.dot(&r1) < 0. {
            r1.mult_mut(-1.);
        }
        let cos = r0.dot(&r1);
        let rotating = cos < 0.9995;
        if rotating {
            let theta = f64::acos(clamp(cos, -1., 1.));
            let q_perp = ((r1-r0).mult(cos)).norm();

            let qperpx = q_perp.v.x;
            let qperpy = q_perp.v.y;
            let qperpz = q_perp.v.z;
            let qperpw = q_perp.w;
            let t0x = t0.x;
            let t0y = t0.y;
            let t0z = t0.z;
            let t1x = t1.x;
            let t1y = t1.y;
            let t1z = t1.z;
            let q0x = r0.v.x;
            let q0y = r0.v.y;
            let q0z = r0.v.z;
            let q0w = r0.w;
            let s000 = s0.m[0][0];
            let s001 = s0.m[0][1];
            let s002 = s0.m[0][2];
            let s010 = s0.m[1][0];
            let s011 = s0.m[1][1];
            let s012 = s0.m[1][2];
            let s020 = s0.m[2][0];
            let s021 = s0.m[2][1];
            let s022 = s0.m[2][2];
            let s100 = s1.m[0][0];
            let s101 = s1.m[0][1];
            let s102 = s1.m[0][2];
            let s110 = s1.m[1][0];
            let s111 = s1.m[1][1];
            let s112 = s1.m[1][2];
            let s120 = s1.m[2][0];
            let s121 = s1.m[2][1];
            let s122 = s1.m[2][2];

            let c10 = DerivativeTerm::new(
                -t0x + t1x,
            (-1. + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) *
                    s000 +
                q0w * q0z * s010 - qperpx * qperpy * s010 +
                qperpw * qperpz * s010 - q0w * q0y * s020 -
                qperpw * qperpy * s020 - qperpx * qperpz * s020 + s100 -
                q0y * q0y * s100 - q0z * q0z * s100 - qperpy * qperpy * s100 -
                qperpz * qperpz * s100 - q0w * q0z * s110 +
                qperpx * qperpy * s110 - qperpw * qperpz * s110 +
                q0w * q0y * s120 + qperpw * qperpy * s120 +
                qperpx * qperpz * s120 +
                q0x * (-(q0y * s010) - q0z * s020 + q0y * s110 + q0z * s120),
            (-1. + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) *
                    s001 +
                q0w * q0z * s011 - qperpx * qperpy * s011 +
                qperpw * qperpz * s011 - q0w * q0y * s021 -
                qperpw * qperpy * s021 - qperpx * qperpz * s021 + s101 -
                q0y * q0y * s101 - q0z * q0z * s101 - qperpy * qperpy * s101 -
                qperpz * qperpz * s101 - q0w * q0z * s111 +
                qperpx * qperpy * s111 - qperpw * qperpz * s111 +
                q0w * q0y * s121 + qperpw * qperpy * s121 +
                qperpx * qperpz * s121 +
                q0x * (-(q0y * s011) - q0z * s021 + q0y * s111 + q0z * s121),
            (-1. + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) *
                    s002 +
                q0w * q0z * s012 - qperpx * qperpy * s012 +
                qperpw * qperpz * s012 - q0w * q0y * s022 -
                qperpw * qperpy * s022 - qperpx * qperpz * s022 + s102 -
                q0y * q0y * s102 - q0z * q0z * s102 - qperpy * qperpy * s102 -
                qperpz * qperpz * s102 - q0w * q0z * s112 +
                qperpx * qperpy * s112 - qperpw * qperpz * s112 +
                q0w * q0y * s122 + qperpw * qperpy * s122 +
                qperpx * qperpz * s122 +
                q0x * (-(q0y * s012) - q0z * s022 + q0y * s112 + q0z * s122)
            );

            let c20 = DerivativeTerm::new(
                0.,
            -(qperpy * qperpy * s000) - qperpz * qperpz * s000 +
                qperpx * qperpy * s010 - qperpw * qperpz * s010 +
                qperpw * qperpy * s020 + qperpx * qperpz * s020 +
                q0y * q0y * (s000 - s100) + q0z * q0z * (s000 - s100) +
                qperpy * qperpy * s100 + qperpz * qperpz * s100 -
                qperpx * qperpy * s110 + qperpw * qperpz * s110 -
                qperpw * qperpy * s120 - qperpx * qperpz * s120 +
                2. * q0x * qperpy * s010 * theta -
                2. * q0w * qperpz * s010 * theta +
                2. * q0w * qperpy * s020 * theta +
                2. * q0x * qperpz * s020 * theta +
                q0y *
                    (q0x * (-s010 + s110) + q0w * (-s020 + s120) +
                     2. * (-2. * qperpy * s000 + qperpx * s010 + qperpw * s020) *
                         theta) +
                q0z * (q0w * (s010 - s110) + q0x * (-s020 + s120) -
                       2. * (2. * qperpz * s000 + qperpw * s010 - qperpx * s020) *
                           theta),
            -(qperpy * qperpy * s001) - qperpz * qperpz * s001 +
                qperpx * qperpy * s011 - qperpw * qperpz * s011 +
                qperpw * qperpy * s021 + qperpx * qperpz * s021 +
                q0y * q0y * (s001 - s101) + q0z * q0z * (s001 - s101) +
                qperpy * qperpy * s101 + qperpz * qperpz * s101 -
                qperpx * qperpy * s111 + qperpw * qperpz * s111 -
                qperpw * qperpy * s121 - qperpx * qperpz * s121 +
                2. * q0x * qperpy * s011 * theta -
                2. * q0w * qperpz * s011 * theta +
                2. * q0w * qperpy * s021 * theta +
                2. * q0x * qperpz * s021 * theta +
                q0y *
                    (q0x * (-s011 + s111) + q0w * (-s021 + s121) +
                     2. * (-2. * qperpy * s001 + qperpx * s011 + qperpw * s021) *
                         theta) +
                q0z * (q0w * (s011 - s111) + q0x * (-s021 + s121) -
                       2. * (2. * qperpz * s001 + qperpw * s011 - qperpx * s021) *
                           theta),
            -(qperpy * qperpy * s002) - qperpz * qperpz * s002 +
                qperpx * qperpy * s012 - qperpw * qperpz * s012 +
                qperpw * qperpy * s022 + qperpx * qperpz * s022 +
                q0y * q0y * (s002 - s102) + q0z * q0z * (s002 - s102) +
                qperpy * qperpy * s102 + qperpz * qperpz * s102 -
                qperpx * qperpy * s112 + qperpw * qperpz * s112 -
                qperpw * qperpy * s122 - qperpx * qperpz * s122 +
                2. * q0x * qperpy * s012 * theta -
                2. * q0w * qperpz * s012 * theta +
                2. * q0w * qperpy * s022 * theta +
                2. * q0x * qperpz * s022 * theta +
                q0y *
                    (q0x * (-s012 + s112) + q0w * (-s022 + s122) +
                     2. * (-2. * qperpy * s002 + qperpx * s012 + qperpw * s022) *
                         theta) +
                q0z * (q0w * (s012 - s112) + q0x * (-s022 + s122) -
                       2. * (2. * qperpz * s002 + qperpw * s012 - qperpx * s022) *
                           theta)
            );

            let c30 = DerivativeTerm::new(
                0.,
                -2. * (q0x * qperpy * s010 - q0w * qperpz * s010 +
                      q0w * qperpy * s020 + q0x * qperpz * s020 -
                      q0x * qperpy * s110 + q0w * qperpz * s110 -
                      q0w * qperpy * s120 - q0x * qperpz * s120 +
                      q0y * (-2. * qperpy * s000 + qperpx * s010 + qperpw * s020 +
                             2. * qperpy * s100 - qperpx * s110 - qperpw * s120) +
                      q0z * (-2. * qperpz * s000 - qperpw * s010 + qperpx * s020 +
                             2. * qperpz * s100 + qperpw * s110 - qperpx * s120)) *
                    theta,
                -2. * (q0x * qperpy * s011 - q0w * qperpz * s011 +
                      q0w * qperpy * s021 + q0x * qperpz * s021 -
                      q0x * qperpy * s111 + q0w * qperpz * s111 -
                      q0w * qperpy * s121 - q0x * qperpz * s121 +
                      q0y * (-2. * qperpy * s001 + qperpx * s011 + qperpw * s021 +
                             2. * qperpy * s101 - qperpx * s111 - qperpw * s121) +
                      q0z * (-2. * qperpz * s001 - qperpw * s011 + qperpx * s021 +
                             2. * qperpz * s101 + qperpw * s111 - qperpx * s121)) *
                    theta,
                -2. * (q0x * qperpy * s012 - q0w * qperpz * s012 +
                      q0w * qperpy * s022 + q0x * qperpz * s022 -
                      q0x * qperpy * s112 + q0w * qperpz * s112 -
                      q0w * qperpy * s122 - q0x * qperpz * s122 +
                      q0y * (-2. * qperpy * s002 + qperpx * s012 + qperpw * s022 +
                             2. * qperpy * s102 - qperpx * s112 - qperpw * s122) +
                      q0z * (-2. * qperpz * s002 - qperpw * s012 + qperpx * s022 +
                             2. * qperpz * s102 + qperpw * s112 - qperpx * s122)) *
                    theta
            );

            let c40 = DerivativeTerm::new(
                0.,
            -(q0x * qperpy * s010) + q0w * qperpz * s010 - q0w * qperpy * s020 -
                q0x * qperpz * s020 + q0x * qperpy * s110 -
                q0w * qperpz * s110 + q0w * qperpy * s120 +
                q0x * qperpz * s120 + 2. * q0y * q0y * s000 * theta +
                2. * q0z * q0z * s000 * theta -
                2. * qperpy * qperpy * s000 * theta -
                2. * qperpz * qperpz * s000 * theta +
                2. * qperpx * qperpy * s010 * theta -
                2. * qperpw * qperpz * s010 * theta +
                2. * qperpw * qperpy * s020 * theta +
                2. * qperpx * qperpz * s020 * theta +
                q0y * (-(qperpx * s010) - qperpw * s020 +
                       2. * qperpy * (s000 - s100) + qperpx * s110 +
                       qperpw * s120 - 2. * q0x * s010 * theta -
                       2. * q0w * s020 * theta) +
                q0z * (2. * qperpz * s000 + qperpw * s010 - qperpx * s020 -
                       2. * qperpz * s100 - qperpw * s110 + qperpx * s120 +
                       2. * q0w * s010 * theta - 2. * q0x * s020 * theta),
            -(q0x * qperpy * s011) + q0w * qperpz * s011 - q0w * qperpy * s021 -
                q0x * qperpz * s021 + q0x * qperpy * s111 -
                q0w * qperpz * s111 + q0w * qperpy * s121 +
                q0x * qperpz * s121 + 2. * q0y * q0y * s001 * theta +
                2. * q0z * q0z * s001 * theta -
                2. * qperpy * qperpy * s001 * theta -
                2. * qperpz * qperpz * s001 * theta +
                2. * qperpx * qperpy * s011 * theta -
                2. * qperpw * qperpz * s011 * theta +
                2. * qperpw * qperpy * s021 * theta +
                2. * qperpx * qperpz * s021 * theta +
                q0y * (-(qperpx * s011) - qperpw * s021 +
                       2. * qperpy * (s001 - s101) + qperpx * s111 +
                       qperpw * s121 - 2. * q0x * s011 * theta -
                       2. * q0w * s021 * theta) +
                q0z * (2. * qperpz * s001 + qperpw * s011 - qperpx * s021 -
                       2. * qperpz * s101 - qperpw * s111 + qperpx * s121 +
                       2. * q0w * s011 * theta - 2. * q0x * s021 * theta),
            -(q0x * qperpy * s012) + q0w * qperpz * s012 - q0w * qperpy * s022 -
                q0x * qperpz * s022 + q0x * qperpy * s112 -
                q0w * qperpz * s112 + q0w * qperpy * s122 +
                q0x * qperpz * s122 + 2. * q0y * q0y * s002 * theta +
                2. * q0z * q0z * s002 * theta -
                2. * qperpy * qperpy * s002 * theta -
                2. * qperpz * qperpz * s002 * theta +
                2. * qperpx * qperpy * s012 * theta -
                2. * qperpw * qperpz * s012 * theta +
                2. * qperpw * qperpy * s022 * theta +
                2. * qperpx * qperpz * s022 * theta +
                q0y * (-(qperpx * s012) - qperpw * s022 +
                       2. * qperpy * (s002 - s102) + qperpx * s112 +
                       qperpw * s122 - 2. * q0x * s012 * theta -
                       2. * q0w * s022 * theta) +
                q0z * (2. * qperpz * s002 + qperpw * s012 - qperpx * s022 -
                       2. * qperpz * s102 - qperpw * s112 + qperpx * s122 +
                       2. * q0w * s012 * theta - 2. * q0x * s022 * theta)
            );


            let c50 = DerivativeTerm::new(
                0.,
                2. * (qperpy * qperpy * s000 + qperpz * qperpz * s000 -
                    qperpx * qperpy * s010 + qperpw * qperpz * s010 -
                    qperpw * qperpy * s020 - qperpx * qperpz * s020 -
                    qperpy * qperpy * s100 - qperpz * qperpz * s100 +
                    q0y * q0y * (-s000 + s100) + q0z * q0z * (-s000 + s100) +
                    qperpx * qperpy * s110 - qperpw * qperpz * s110 +
                    q0y * (q0x * (s010 - s110) + q0w * (s020 - s120)) +
                    qperpw * qperpy * s120 + qperpx * qperpz * s120 +
                    q0z * (-(q0w * s010) + q0x * s020 + q0w * s110 - q0x * s120)) *
                    theta,
                2. * (qperpy * qperpy * s001 + qperpz * qperpz * s001 -
                    qperpx * qperpy * s011 + qperpw * qperpz * s011 -
                    qperpw * qperpy * s021 - qperpx * qperpz * s021 -
                    qperpy * qperpy * s101 - qperpz * qperpz * s101 +
                    q0y * q0y * (-s001 + s101) + q0z * q0z * (-s001 + s101) +
                    qperpx * qperpy * s111 - qperpw * qperpz * s111 +
                    q0y * (q0x * (s011 - s111) + q0w * (s021 - s121)) +
                    qperpw * qperpy * s121 + qperpx * qperpz * s121 +
                    q0z * (-(q0w * s011) + q0x * s021 + q0w * s111 - q0x * s121)) *
                    theta,
                2. * (qperpy * qperpy * s002 + qperpz * qperpz * s002 -
                    qperpx * qperpy * s012 + qperpw * qperpz * s012 -
                    qperpw * qperpy * s022 - qperpx * qperpz * s022 -
                    qperpy * qperpy * s102 - qperpz * qperpz * s102 +
                    q0y * q0y * (-s002 + s102) + q0z * q0z * (-s002 + s102) +
                    qperpx * qperpy * s112 - qperpw * qperpz * s112 +
                    q0y * (q0x * (s012 - s112) + q0w * (s022 - s122)) +
                    qperpw * qperpy * s122 + qperpx * qperpz * s122 +
                    q0z * (-(q0w * s012) + q0x * s022 + q0w * s112 - q0x * s122)) *
                    theta
            );

            let c11 = DerivativeTerm::new(
                -t0y + t1y,
                -(qperpx * qperpy * s000) - qperpw * qperpz * s000 - s010 +
                    q0z * q0z * s010 + qperpx * qperpx * s010 +
                    qperpz * qperpz * s010 - q0y * q0z * s020 +
                    qperpw * qperpx * s020 - qperpy * qperpz * s020 +
                    qperpx * qperpy * s100 + qperpw * qperpz * s100 +
                    q0w * q0z * (-s000 + s100) + q0x * q0x * (s010 - s110) + s110 -
                    q0z * q0z * s110 - qperpx * qperpx * s110 -
                    qperpz * qperpz * s110 +
                    q0x * (q0y * (-s000 + s100) + q0w * (s020 - s120)) +
                    q0y * q0z * s120 - qperpw * qperpx * s120 +
                    qperpy * qperpz * s120,
                -(qperpx * qperpy * s001) - qperpw * qperpz * s001 - s011 +
                    q0z * q0z * s011 + qperpx * qperpx * s011 +
                    qperpz * qperpz * s011 - q0y * q0z * s021 +
                    qperpw * qperpx * s021 - qperpy * qperpz * s021 +
                    qperpx * qperpy * s101 + qperpw * qperpz * s101 +
                    q0w * q0z * (-s001 + s101) + q0x * q0x * (s011 - s111) + s111 -
                    q0z * q0z * s111 - qperpx * qperpx * s111 -
                    qperpz * qperpz * s111 +
                    q0x * (q0y * (-s001 + s101) + q0w * (s021 - s121)) +
                    q0y * q0z * s121 - qperpw * qperpx * s121 +
                    qperpy * qperpz * s121,
                -(qperpx * qperpy * s002) - qperpw * qperpz * s002 - s012 +
                    q0z * q0z * s012 + qperpx * qperpx * s012 +
                    qperpz * qperpz * s012 - q0y * q0z * s022 +
                    qperpw * qperpx * s022 - qperpy * qperpz * s022 +
                    qperpx * qperpy * s102 + qperpw * qperpz * s102 +
                    q0w * q0z * (-s002 + s102) + q0x * q0x * (s012 - s112) + s112 -
                    q0z * q0z * s112 - qperpx * qperpx * s112 -
                    qperpz * qperpz * s112 +
                    q0x * (q0y * (-s002 + s102) + q0w * (s022 - s122)) +
                    q0y * q0z * s122 - qperpw * qperpx * s122 +
                    qperpy * qperpz * s122
            );

            let c21 = DerivativeTerm::new(
                0.,
                qperpx * qperpy * s000 + qperpw * qperpz * s000 + q0z * q0z * s010 -
                    qperpx * qperpx * s010 - qperpz * qperpz * s010 -
                    q0y * q0z * s020 - qperpw * qperpx * s020 +
                    qperpy * qperpz * s020 - qperpx * qperpy * s100 -
                    qperpw * qperpz * s100 + q0x * q0x * (s010 - s110) -
                    q0z * q0z * s110 + qperpx * qperpx * s110 +
                    qperpz * qperpz * s110 + q0y * q0z * s120 +
                    qperpw * qperpx * s120 - qperpy * qperpz * s120 +
                    2. * q0z * qperpw * s000 * theta +
                    2. * q0y * qperpx * s000 * theta -
                    4. * q0z * qperpz * s010 * theta +
                    2. * q0z * qperpy * s020 * theta +
                    2. * q0y * qperpz * s020 * theta +
                    q0x * (q0w * s020 + q0y * (-s000 + s100) - q0w * s120 +
                        2. * qperpy * s000 * theta - 4. * qperpx * s010 * theta -
                        2. * qperpw * s020 * theta) +
                    q0w * (-(q0z * s000) + q0z * s100 + 2. * qperpz * s000 * theta -
                        2. * qperpx * s020 * theta),
                qperpx * qperpy * s001 + qperpw * qperpz * s001 + q0z * q0z * s011 -
                    qperpx * qperpx * s011 - qperpz * qperpz * s011 -
                    q0y * q0z * s021 - qperpw * qperpx * s021 +
                    qperpy * qperpz * s021 - qperpx * qperpy * s101 -
                    qperpw * qperpz * s101 + q0x * q0x * (s011 - s111) -
                    q0z * q0z * s111 + qperpx * qperpx * s111 +
                    qperpz * qperpz * s111 + q0y * q0z * s121 +
                    qperpw * qperpx * s121 - qperpy * qperpz * s121 +
                    2. * q0z * qperpw * s001 * theta +
                    2. * q0y * qperpx * s001 * theta -
                    4. * q0z * qperpz * s011 * theta +
                    2. * q0z * qperpy * s021 * theta +
                    2. * q0y * qperpz * s021 * theta +
                    q0x * (q0w * s021 + q0y * (-s001 + s101) - q0w * s121 +
                        2. * qperpy * s001 * theta - 4. * qperpx * s011 * theta -
                        2. * qperpw * s021 * theta) +
                    q0w * (-(q0z * s001) + q0z * s101 + 2. * qperpz * s001 * theta -
                        2. * qperpx * s021 * theta),
                qperpx * qperpy * s002 + qperpw * qperpz * s002 + q0z * q0z * s012 -
                    qperpx * qperpx * s012 - qperpz * qperpz * s012 -
                    q0y * q0z * s022 - qperpw * qperpx * s022 +
                    qperpy * qperpz * s022 - qperpx * qperpy * s102 -
                    qperpw * qperpz * s102 + q0x * q0x * (s012 - s112) -
                    q0z * q0z * s112 + qperpx * qperpx * s112 +
                    qperpz * qperpz * s112 + q0y * q0z * s122 +
                    qperpw * qperpx * s122 - qperpy * qperpz * s122 +
                    2. * q0z * qperpw * s002 * theta +
                    2. * q0y * qperpx * s002 * theta -
                    4. * q0z * qperpz * s012 * theta +
                    2. * q0z * qperpy * s022 * theta +
                    2. * q0y * qperpz * s022 * theta +
                    q0x * (q0w * s022 + q0y * (-s002 + s102) - q0w * s122 +
                        2. * qperpy * s002 * theta - 4. * qperpx * s012 * theta -
                        2. * qperpw * s022 * theta) +
                    q0w * (-(q0z * s002) + q0z * s102 + 2. * qperpz * s002 * theta -
                        2. * qperpx * s022 * theta)
            );

            let c31 = DerivativeTerm::new(
                0., 2. * (-(q0x * qperpy * s000) - q0w * qperpz * s000 +
                2. * q0x * qperpx * s010 + q0x * qperpw * s020 +
                q0w * qperpx * s020 + q0x * qperpy * s100 +
                q0w * qperpz * s100 - 2. * q0x * qperpx * s110 -
                q0x * qperpw * s120 - q0w * qperpx * s120 +
                q0z * (2. * qperpz * s010 - qperpy * s020 +
                    qperpw * (-s000 + s100) - 2. * qperpz * s110 +
                    qperpy * s120) +
                q0y * (-(qperpx * s000) - qperpz * s020 + qperpx * s100 +
                    qperpz * s120)) *
            theta,
                2. * (-(q0x * qperpy * s001) - q0w * qperpz * s001 +
                2. * q0x * qperpx * s011 + q0x * qperpw * s021 +
                q0w * qperpx * s021 + q0x * qperpy * s101 +
                q0w * qperpz * s101 - 2. * q0x * qperpx * s111 -
                q0x * qperpw * s121 - q0w * qperpx * s121 +
                q0z * (2. * qperpz * s011 - qperpy * s021 +
                qperpw * (-s001 + s101) - 2. * qperpz * s111 +
                qperpy * s121) +
                q0y * (-(qperpx * s001) - qperpz * s021 + qperpx * s101 +
                qperpz * s121)) *
                theta,
                2. * (-(q0x * qperpy * s002) - q0w * qperpz * s002 +
                2. * q0x * qperpx * s012 + q0x * qperpw * s022 +
                q0w * qperpx * s022 + q0x * qperpy * s102 +
                q0w * qperpz * s102 - 2. * q0x * qperpx * s112 -
                q0x * qperpw * s122 - q0w * qperpx * s122 +
                q0z * (2. * qperpz * s012 - qperpy * s022 +
                qperpw * (-s002 + s102) - 2. * qperpz * s112 +
                qperpy * s122) +
                q0y * (-(qperpx * s002) - qperpz * s022 + qperpx * s102 +
                qperpz * s122)) *
                theta
            );

            let c41 = DerivativeTerm::new(
                0.,
                -(q0x * qperpy * s000) - q0w * qperpz * s000 +
                    2. * q0x * qperpx * s010 + q0x * qperpw * s020 +
                    q0w * qperpx * s020 + q0x * qperpy * s100 +
                    q0w * qperpz * s100 - 2. * q0x * qperpx * s110 -
                    q0x * qperpw * s120 - q0w * qperpx * s120 +
                    2. * qperpx * qperpy * s000 * theta +
                    2. * qperpw * qperpz * s000 * theta +
                    2. * q0x * q0x * s010 * theta + 2. * q0z * q0z * s010 * theta -
                    2. * qperpx * qperpx * s010 * theta -
                    2. * qperpz * qperpz * s010 * theta +
                    2. * q0w * q0x * s020 * theta -
                    2. * qperpw * qperpx * s020 * theta +
                    2. * qperpy * qperpz * s020 * theta +
                    q0y * (-(qperpx * s000) - qperpz * s020 + qperpx * s100 +
                        qperpz * s120 - 2. * q0x * s000 * theta) +
                    q0z * (2. * qperpz * s010 - qperpy * s020 +
                        qperpw * (-s000 + s100) - 2. * qperpz * s110 +
                        qperpy * s120 - 2. * q0w * s000 * theta -
                        2. * q0y * s020 * theta),
                -(q0x * qperpy * s001) - q0w * qperpz * s001 +
                    2. * q0x * qperpx * s011 + q0x * qperpw * s021 +
                    q0w * qperpx * s021 + q0x * qperpy * s101 +
                    q0w * qperpz * s101 - 2. * q0x * qperpx * s111 -
                    q0x * qperpw * s121 - q0w * qperpx * s121 +
                    2. * qperpx * qperpy * s001 * theta +
                    2. * qperpw * qperpz * s001 * theta +
                    2. * q0x * q0x * s011 * theta + 2. * q0z * q0z * s011 * theta -
                    2. * qperpx * qperpx * s011 * theta -
                    2. * qperpz * qperpz * s011 * theta +
                    2. * q0w * q0x * s021 * theta -
                    2. * qperpw * qperpx * s021 * theta +
                    2. * qperpy * qperpz * s021 * theta +
                    q0y * (-(qperpx * s001) - qperpz * s021 + qperpx * s101 +
                        qperpz * s121 - 2. * q0x * s001 * theta) +
                    q0z * (2. * qperpz * s011 - qperpy * s021 +
                        qperpw * (-s001 + s101) - 2. * qperpz * s111 +
                        qperpy * s121 - 2. * q0w * s001 * theta -
                        2. * q0y * s021 * theta),
                -(q0x * qperpy * s002) - q0w * qperpz * s002 +
                    2. * q0x * qperpx * s012 + q0x * qperpw * s022 +
                    q0w * qperpx * s022 + q0x * qperpy * s102 +
                    q0w * qperpz * s102 - 2. * q0x * qperpx * s112 -
                    q0x * qperpw * s122 - q0w * qperpx * s122 +
                    2. * qperpx * qperpy * s002 * theta +
                    2. * qperpw * qperpz * s002 * theta +
                    2. * q0x * q0x * s012 * theta + 2. * q0z * q0z * s012 * theta -
                    2. * qperpx * qperpx * s012 * theta -
                    2. * qperpz * qperpz * s012 * theta +
                    2. * q0w * q0x * s022 * theta -
                    2. * qperpw * qperpx * s022 * theta +
                    2. * qperpy * qperpz * s022 * theta +
                    q0y * (-(qperpx * s002) - qperpz * s022 + qperpx * s102 +
                        qperpz * s122 - 2. * q0x * s002 * theta) +
                    q0z * (2. * qperpz * s012 - qperpy * s022 +
                        qperpw * (-s002 + s102) - 2. * qperpz * s112 +
                        qperpy * s122 - 2. * q0w * s002 * theta -
                        2. * q0y * s022 * theta)
            );

            let c51 = DerivativeTerm::new(
                0., -2. * (qperpx * qperpy * s000 + qperpw * qperpz * s000 +
                    q0z * q0z * s010 - qperpx * qperpx * s010 -
                    qperpz * qperpz * s010 - q0y * q0z * s020 -
                    qperpw * qperpx * s020 + qperpy * qperpz * s020 -
                    qperpx * qperpy * s100 - qperpw * qperpz * s100 +
                    q0w * q0z * (-s000 + s100) + q0x * q0x * (s010 - s110) -
                    q0z * q0z * s110 + qperpx * qperpx * s110 +
                    qperpz * qperpz * s110 +
                    q0x * (q0y * (-s000 + s100) + q0w * (s020 - s120)) +
                    q0y * q0z * s120 + qperpw * qperpx * s120 -
                    qperpy * qperpz * s120) *
                theta,
                -2. * (qperpx * qperpy * s001 + qperpw * qperpz * s001 +
                q0z * q0z * s011 - qperpx * qperpx * s011 -
                qperpz * qperpz * s011 - q0y * q0z * s021 -
                qperpw * qperpx * s021 + qperpy * qperpz * s021 -
                qperpx * qperpy * s101 - qperpw * qperpz * s101 +
                q0w * q0z * (-s001 + s101) + q0x * q0x * (s011 - s111) -
                q0z * q0z * s111 + qperpx * qperpx * s111 +
                qperpz * qperpz * s111 +
                q0x * (q0y * (-s001 + s101) + q0w * (s021 - s121)) +
                q0y * q0z * s121 + qperpw * qperpx * s121 -
                qperpy * qperpz * s121) *
                theta,
                -2. * (qperpx * qperpy * s002 + qperpw * qperpz * s002 +
                q0z * q0z * s012 - qperpx * qperpx * s012 -
                qperpz * qperpz * s012 - q0y * q0z * s022 -
                qperpw * qperpx * s022 + qperpy * qperpz * s022 -
                qperpx * qperpy * s102 - qperpw * qperpz * s102 +
                q0w * q0z * (-s002 + s102) + q0x * q0x * (s012 - s112) -
                q0z * q0z * s112 + qperpx * qperpx * s112 +
                qperpz * qperpz * s112 +
                q0x * (q0y * (-s002 + s102) + q0w * (s022 - s122)) +
                q0y * q0z * s122 + qperpw * qperpx * s122 -
                qperpy * qperpz * s122) *
                theta
            );

            let c12 = DerivativeTerm::new(
                -t0z + t1z, qperpw * qperpy * s000 - qperpx * qperpz * s000 -
                q0y * q0z * s010 - qperpw * qperpx * s010 -
                qperpy * qperpz * s010 - s020 + q0y * q0y * s020 +
                qperpx * qperpx * s020 + qperpy * qperpy * s020 -
                qperpw * qperpy * s100 + qperpx * qperpz * s100 +
                q0x * q0z * (-s000 + s100) + q0y * q0z * s110 +
                qperpw * qperpx * s110 + qperpy * qperpz * s110 +
                q0w * (q0y * (s000 - s100) + q0x * (-s010 + s110)) +
                q0x * q0x * (s020 - s120) + s120 - q0y * q0y * s120 -
                qperpx * qperpx * s120 - qperpy * qperpy * s120,
                qperpw * qperpy * s001 - qperpx * qperpz * s001 -
                q0y * q0z * s011 - qperpw * qperpx * s011 -
                qperpy * qperpz * s011 - s021 + q0y * q0y * s021 +
                qperpx * qperpx * s021 + qperpy * qperpy * s021 -
                qperpw * qperpy * s101 + qperpx * qperpz * s101 +
                q0x * q0z * (-s001 + s101) + q0y * q0z * s111 +
                qperpw * qperpx * s111 + qperpy * qperpz * s111 +
                q0w * (q0y * (s001 - s101) + q0x * (-s011 + s111)) +
                q0x * q0x * (s021 - s121) + s121 - q0y * q0y * s121 -
                qperpx * qperpx * s121 - qperpy * qperpy * s121,
                qperpw * qperpy * s002 - qperpx * qperpz * s002 -
                q0y * q0z * s012 - qperpw * qperpx * s012 -
                qperpy * qperpz * s012 - s022 + q0y * q0y * s022 +
                qperpx * qperpx * s022 + qperpy * qperpy * s022 -
                qperpw * qperpy * s102 + qperpx * qperpz * s102 +
                q0x * q0z * (-s002 + s102) + q0y * q0z * s112 +
                qperpw * qperpx * s112 + qperpy * qperpz * s112 +
                q0w * (q0y * (s002 - s102) + q0x * (-s012 + s112)) +
                q0x * q0x * (s022 - s122) + s122 - q0y * q0y * s122 -
                qperpx * qperpx * s122 - qperpy * qperpy * s122
            );

            let c22 = DerivativeTerm::new(
                0.,
                q0w * q0y * s000 - q0x * q0z * s000 - qperpw * qperpy * s000 +
                qperpx * qperpz * s000 - q0w * q0x * s010 - q0y * q0z * s010 +
                qperpw * qperpx * s010 + qperpy * qperpz * s010 +
                q0x * q0x * s020 + q0y * q0y * s020 - qperpx * qperpx * s020 -
                qperpy * qperpy * s020 - q0w * q0y * s100 + q0x * q0z * s100 +
                qperpw * qperpy * s100 - qperpx * qperpz * s100 +
                q0w * q0x * s110 + q0y * q0z * s110 - qperpw * qperpx * s110 -
                qperpy * qperpz * s110 - q0x * q0x * s120 - q0y * q0y * s120 +
                qperpx * qperpx * s120 + qperpy * qperpy * s120 -
                2. * q0y * qperpw * s000 * theta + 2. * q0z * qperpx * s000 * theta -
                2. * q0w * qperpy * s000 * theta + 2. * q0x * qperpz * s000 * theta +
                2. * q0x * qperpw * s010 * theta + 2. * q0w * qperpx * s010 * theta +
                2. * q0z * qperpy * s010 * theta + 2. * q0y * qperpz * s010 * theta -
                4. * q0x * qperpx * s020 * theta - 4. * q0y * qperpy * s020 * theta,
                q0w * q0y * s001 - q0x * q0z * s001 - qperpw * qperpy * s001 +
                qperpx * qperpz * s001 - q0w * q0x * s011 - q0y * q0z * s011 +
                qperpw * qperpx * s011 + qperpy * qperpz * s011 +
                q0x * q0x * s021 + q0y * q0y * s021 - qperpx * qperpx * s021 -
                qperpy * qperpy * s021 - q0w * q0y * s101 + q0x * q0z * s101 +
                qperpw * qperpy * s101 - qperpx * qperpz * s101 +
                q0w * q0x * s111 + q0y * q0z * s111 - qperpw * qperpx * s111 -
                qperpy * qperpz * s111 - q0x * q0x * s121 - q0y * q0y * s121 +
                qperpx * qperpx * s121 + qperpy * qperpy * s121 -
                2. * q0y * qperpw * s001 * theta + 2. * q0z * qperpx * s001 * theta -
                2. * q0w * qperpy * s001 * theta + 2. * q0x * qperpz * s001 * theta +
                2. * q0x * qperpw * s011 * theta + 2. * q0w * qperpx * s011 * theta +
                2. * q0z * qperpy * s011 * theta + 2. * q0y * qperpz * s011 * theta -
                4. * q0x * qperpx * s021 * theta - 4. * q0y * qperpy * s021 * theta,
                q0w * q0y * s002 - q0x * q0z * s002 - qperpw * qperpy * s002 +
                qperpx * qperpz * s002 - q0w * q0x * s012 - q0y * q0z * s012 +
                qperpw * qperpx * s012 + qperpy * qperpz * s012 +
                q0x * q0x * s022 + q0y * q0y * s022 - qperpx * qperpx * s022 -
                qperpy * qperpy * s022 - q0w * q0y * s102 + q0x * q0z * s102 +
                qperpw * qperpy * s102 - qperpx * qperpz * s102 +
                q0w * q0x * s112 + q0y * q0z * s112 - qperpw * qperpx * s112 -
                qperpy * qperpz * s112 - q0x * q0x * s122 - q0y * q0y * s122 +
                qperpx * qperpx * s122 + qperpy * qperpy * s122 -
                2. * q0y * qperpw * s002 * theta + 2. * q0z * qperpx * s002 * theta -
                2. * q0w * qperpy * s002 * theta + 2. * q0x * qperpz * s002 * theta +
                2. * q0x * qperpw * s012 * theta + 2. * q0w * qperpx * s012 * theta +
                2. * q0z * qperpy * s012 * theta + 2. * q0y * qperpz * s012 * theta -
                4. * q0x * qperpx * s022 * theta -
                4. * q0y * qperpy * s022 * theta
            );

            let c32 = DerivativeTerm::new(
                0., -2. * (-(q0w * qperpy * s000) + q0x * qperpz * s000 +
                        q0x * qperpw * s010 + q0w * qperpx * s010 -
                        2. * q0x * qperpx * s020 + q0w * qperpy * s100 -
                        q0x * qperpz * s100 - q0x * qperpw * s110 -
                        q0w * qperpx * s110 +
                        q0z * (qperpx * s000 + qperpy * s010 - qperpx * s100 -
                                qperpy * s110) +
                        2. * q0x * qperpx * s120 +
                        q0y * (qperpz * s010 - 2. * qperpy * s020 +
                                qperpw * (-s000 + s100) - qperpz * s110 +
                                2. * qperpy * s120)) *
                        theta,
                -2. * (-(q0w * qperpy * s001) + q0x * qperpz * s001 +
                    q0x * qperpw * s011 + q0w * qperpx * s011 -
                    2. * q0x * qperpx * s021 + q0w * qperpy * s101 -
                    q0x * qperpz * s101 - q0x * qperpw * s111 -
                    q0w * qperpx * s111 +
                    q0z * (qperpx * s001 + qperpy * s011 - qperpx * s101 -
                            qperpy * s111) +
                    2. * q0x * qperpx * s121 +
                    q0y * (qperpz * s011 - 2. * qperpy * s021 +
                            qperpw * (-s001 + s101) - qperpz * s111 +
                            2. * qperpy * s121)) *
                    theta,
                -2. * (-(q0w * qperpy * s002) + q0x * qperpz * s002 +
                    q0x * qperpw * s012 + q0w * qperpx * s012 -
                    2. * q0x * qperpx * s022 + q0w * qperpy * s102 -
                    q0x * qperpz * s102 - q0x * qperpw * s112 -
                    q0w * qperpx * s112 +
                    q0z * (qperpx * s002 + qperpy * s012 - qperpx * s102 -
                            qperpy * s112) +
                    2. * q0x * qperpx * s122 +
                    q0y * (qperpz * s012 - 2. * qperpy * s022 +
                            qperpw * (-s002 + s102) - qperpz * s112 +
                            2. * qperpy * s122)) *
                    theta
            );

            let c42 = DerivativeTerm::new(
                0.,
                q0w * qperpy * s000 - q0x * qperpz * s000 - q0x * qperpw * s010 -
                    q0w * qperpx * s010 + 2. * q0x * qperpx * s020 -
                    q0w * qperpy * s100 + q0x * qperpz * s100 +
                    q0x * qperpw * s110 + q0w * qperpx * s110 -
                    2. * q0x * qperpx * s120 - 2. * qperpw * qperpy * s000 * theta +
                    2. * qperpx * qperpz * s000 * theta -
                    2. * q0w * q0x * s010 * theta +
                    2. * qperpw * qperpx * s010 * theta +
                    2. * qperpy * qperpz * s010 * theta +
                    2. * q0x * q0x * s020 * theta + 2. * q0y * q0y * s020 * theta -
                    2. * qperpx * qperpx * s020 * theta -
                    2. * qperpy * qperpy * s020 * theta +
                    q0z * (-(qperpx * s000) - qperpy * s010 + qperpx * s100 +
                        qperpy * s110 - 2. * q0x * s000 * theta) +
                    q0y * (-(qperpz * s010) + 2. * qperpy * s020 +
                        qperpw * (s000 - s100) + qperpz * s110 -
                        2. * qperpy * s120 + 2. * q0w * s000 * theta -
                        2. * q0z * s010 * theta),
                q0w * qperpy * s001 - q0x * qperpz * s001 - q0x * qperpw * s011 -
                    q0w * qperpx * s011 + 2. * q0x * qperpx * s021 -
                    q0w * qperpy * s101 + q0x * qperpz * s101 +
                    q0x * qperpw * s111 + q0w * qperpx * s111 -
                    2. * q0x * qperpx * s121 - 2. * qperpw * qperpy * s001 * theta +
                    2. * qperpx * qperpz * s001 * theta -
                    2. * q0w * q0x * s011 * theta +
                    2. * qperpw * qperpx * s011 * theta +
                    2. * qperpy * qperpz * s011 * theta +
                    2. * q0x * q0x * s021 * theta + 2. * q0y * q0y * s021 * theta -
                    2. * qperpx * qperpx * s021 * theta -
                    2. * qperpy * qperpy * s021 * theta +
                    q0z * (-(qperpx * s001) - qperpy * s011 + qperpx * s101 +
                        qperpy * s111 - 2. * q0x * s001 * theta) +
                    q0y * (-(qperpz * s011) + 2. * qperpy * s021 +
                        qperpw * (s001 - s101) + qperpz * s111 -
                        2. * qperpy * s121 + 2. * q0w * s001 * theta -
                        2. * q0z * s011 * theta),
                q0w * qperpy * s002 - q0x * qperpz * s002 - q0x * qperpw * s012 -
                    q0w * qperpx * s012 + 2. * q0x * qperpx * s022 -
                    q0w * qperpy * s102 + q0x * qperpz * s102 +
                    q0x * qperpw * s112 + q0w * qperpx * s112 -
                    2. * q0x * qperpx * s122 - 2. * qperpw * qperpy * s002 * theta +
                    2. * qperpx * qperpz * s002 * theta -
                    2. * q0w * q0x * s012 * theta +
                    2. * qperpw * qperpx * s012 * theta +
                    2. * qperpy * qperpz * s012 * theta +
                    2. * q0x * q0x * s022 * theta + 2. * q0y * q0y * s022 * theta -
                    2. * qperpx * qperpx * s022 * theta -
                    2. * qperpy * qperpy * s022 * theta +
                    q0z * (-(qperpx * s002) - qperpy * s012 + qperpx * s102 +
                        qperpy * s112 - 2. * q0x * s002 * theta) +
                    q0y * (-(qperpz * s012) + 2. * qperpy * s022 +
                        qperpw * (s002 - s102) + qperpz * s112 -
                        2. * qperpy * s122 + 2. * q0w * s002 * theta -
                        2. * q0z * s012 * theta)
            );

            let c52 = DerivativeTerm::new(
                0., 2. * (qperpw * qperpy * s000 - qperpx * qperpz * s000 +
                    q0y * q0z * s010 - qperpw * qperpx * s010 -
                    qperpy * qperpz * s010 - q0y * q0y * s020 +
                    qperpx * qperpx * s020 + qperpy * qperpy * s020 +
                    q0x * q0z * (s000 - s100) - qperpw * qperpy * s100 +
                    qperpx * qperpz * s100 +
                    q0w * (q0y * (-s000 + s100) + q0x * (s010 - s110)) -
                    q0y * q0z * s110 + qperpw * qperpx * s110 +
                    qperpy * qperpz * s110 + q0y * q0y * s120 -
                    qperpx * qperpx * s120 - qperpy * qperpy * s120 +
                    q0x * q0x * (-s020 + s120)) *
                theta,
            2. * (qperpw * qperpy * s001 - qperpx * qperpz * s001 +
                q0y * q0z * s011 - qperpw * qperpx * s011 -
                qperpy * qperpz * s011 - q0y * q0y * s021 +
                qperpx * qperpx * s021 + qperpy * qperpy * s021 +
                q0x * q0z * (s001 - s101) - qperpw * qperpy * s101 +
                qperpx * qperpz * s101 +
                q0w * (q0y * (-s001 + s101) + q0x * (s011 - s111)) -
                q0y * q0z * s111 + qperpw * qperpx * s111 +
                qperpy * qperpz * s111 + q0y * q0y * s121 -
                qperpx * qperpx * s121 - qperpy * qperpy * s121 +
                q0x * q0x * (-s021 + s121)) *
            theta,
            2. * (qperpw * qperpy * s002 - qperpx * qperpz * s002 +
                q0y * q0z * s012 - qperpw * qperpx * s012 -
                qperpy * qperpz * s012 - q0y * q0y * s022 +
                qperpx * qperpx * s022 + qperpy * qperpy * s022 +
                q0x * q0z * (s002 - s102) - qperpw * qperpy * s102 +
                qperpx * qperpz * s102 +
                q0w * (q0y * (-s002 + s102) + q0x * (s012 - s112)) -
                q0y * q0z * s112 + qperpw * qperpx * s112 +
                qperpy * qperpz * s112 + q0y * q0y * s122 -
                qperpx * qperpx * s122 - qperpy * qperpy * s122 +
                q0x * q0x * (-s022 + s122)) *
            theta
            ); //TODO: fix indentation

            return AnimatedTransform { // we include the entire derivative here
                animated: true,
                rotating: true,
                start_trans,
                end_trans,
                start_time,
                end_time,
                t: [t0, t1],
                r: [r0, r1],
                s: [s0, s1],
                c1: [c10, c11, c12],
                c2: [c20, c21, c22],
                c3: [c30, c31, c32],
                c4: [c40, c41, c42],
                c5: [c50, c51, c52]
            };

        }

        AnimatedTransform { // not rotating so we don't need the derivative
            animated: true,
            rotating: false,
            start_trans,
            end_trans,
            start_time,
            end_time,
            t: [t0, t1],
            r: [r0, r1],
            s: [s0, s1],
            c1: [DerivativeTerm::default(), DerivativeTerm::default(), DerivativeTerm::default()],
            c2: [DerivativeTerm::default(), DerivativeTerm::default(), DerivativeTerm::default()],
            c3: [DerivativeTerm::default(), DerivativeTerm::default(), DerivativeTerm::default()],
            c4: [DerivativeTerm::default(), DerivativeTerm::default(), DerivativeTerm::default()],
            c5: [DerivativeTerm::default(), DerivativeTerm::default(), DerivativeTerm::default()]
        }
        
    }

    fn decompose(m: &Matrix4x4) -> (Vector3, Quaternion, Matrix4x4) {

        let t = Vector3{ // translation
            x: m.m[0][3],
            y: m.m[1][3],
            z: m.m[2][3]
        };

        let mut m = m.clone(); //rebind m without translation
        for i in 0..3 {
            m.m[i][3] = 0.;
            m.m[3][i] = 0.;
        }
        m.m[3][3] = 1.;

        let mut norm = 0.;
        let mut count = 0;
        let mut r = m.clone(); // get an extra copy for messing up
        loop {
            // we repeatedly calculate the next matrix in the series
            // by taking the inverse of the transpose.
            let mut r_next = IDENTITY;
            let mut r_it = r.transpose().inverse();
            for i in 0..4 {
                for j in 0..4 {
                    r_next.m[i][j] = 0.5 * (r.m[i][j] + r_it.m[i][j]);
                }
            }
            norm = 0.;
            for i in 0..3 {
                let n = f64::abs(r.m[i][0] - r_next.m[i][0]) + 
                    f64::abs(r.m[i][1] - r_next.m[i][1]) +
                    f64::abs(r.m[i][2] - r_next.m[i][2]);
                norm = math::max(norm, n);
            }
            r = r_next;

            count += 1;
            if count >= 100 || norm < 0.0001 {
                break;
            }
        }
        let r_q = Quaternion::from_transform(&Transform::from_mat(r.m)); //TODO make this less gross
        
        // compute scale
        let s = r.inverse() * m;

        (t, r_q, s)
    }

    /// Computes the transformation matrix for some time t, returning the boundaries of the time range if a time value outside the range is provided
    fn interpolate(&self, t: f64) -> Transform { //TODO is this mem inefficient?
        if !self.animated || t <= self.start_time {
            return self.start_trans;
        } else if t >= self.end_time {
            return self.end_trans;
        }
        let dt = (t - self.start_time) / (self.end_time - self.start_time);
        let trans = self.t[0].mult(1. - dt) + self.t[1].mult(dt);
        let rot = Quaternion::slerp(dt, &self.r[0], &self.r[1]);
        let mut scale = IDENTITY;
        for i in 0..3 {
            for j in 0..3 {
                scale.m[i][j] = math::lerp(dt, self.s[0].m[i][j], self.s[1].m[i][j]);
            }
        }
        Transform::from_mat((Transform::translate(&trans).m * rot.to_transform().m * Transform::from_mat(scale.m).m).m)
    }

    pub fn trans_ray(&self, r: &Ray) -> Ray {
        if !self.animated || r.time <= self.start_time {
            self.start_trans.trans_ray(r)
        } else if r.time >= self.end_time {
            self.end_trans.trans_ray(r)
        } else {
            self.interpolate(r.time).trans_ray(r)
        }
    }

    pub fn trans_point(&self, time: f64, p: &Point3) -> Point3{
        if !self.animated || time <= self.start_time {
            self.start_trans.trans_point(p)
        } else if time >= self.end_time {
            self.end_trans.trans_point(p)
        } else {
            self.interpolate(time).trans_point(p)
        }
    }

    pub fn trans_vec(&self, time: f64, v: &Vector3) -> Vector3 {
        if !self.animated || time <= self.start_time {
            self.start_trans.trans_vec(v)
        } else if time >= self.end_time {
            self.end_trans.trans_vec(v)
        } else {
            self.interpolate(time).trans_vec(v)
        }
    }

    /// Creates a boundary for the motion of a boundary.
    pub fn motion_bounds(&self, b: &Bounds3) -> Bounds3 {
        if !self.animated {
            return self.start_trans.trans_bounds(b);
        }
        if !self.rotating {
            return self.start_trans.trans_bounds(b).union(&self.end_trans.trans_bounds(b));
        }
        let mut bounds = Bounds3::default();
        for c in 0..8 {
            bounds = bounds.union(&self.bound_pt_motion(&b.corner(c)));
        }
        bounds
    }

    /// Creates a boundary for the motion of a single point, given this animated transformation.
    #[allow(clippy::needless_range_loop)] //clippy is wrong here
    pub fn bound_pt_motion(&self, p: &Point3) -> Bounds3 {
        let mut bounds = Bounds3::new(&self.start_trans.trans_point(p), &self.end_trans.trans_point(p));
        let cos = self.r[0].dot(&self.r[1]);
        let theta = f64::acos(math::clamp(cos, -1., 1.));
        for i in 0..3 {
            let mut zeros = [0., 0., 0., 0.];
            let n_zeros: usize = 0; //this will be mutated below
            Interval::interval_find_zeros((self.c1[i].eval(p), self.c2[i].eval(p), self.c3[i].eval(p), self.c4[i].eval(p),
                    self.c5[i].eval(p)), theta, Interval::new(0., 1.), &mut zeros, n_zeros as *mut usize, 8);
            
            for i in 0..n_zeros {
                bounds = bounds.union_pt(&self.trans_point(math::lerp(zeros[i], self.start_time, self.end_time), p));
            }
        }
        bounds
    }

}

#[derive(Copy, Clone, PartialEq)]
struct Interval {
    low: f64,
    high: f64
}

impl Interval {
    fn new(low: f64, high: f64) -> Interval {
        Interval {low, high}
    } 

    fn from_pt(pt: f64) -> Interval {
        Interval {
            low: pt,
            high: pt
        }
    }

    fn sin(&self) -> Interval {
        let mut sin_low = f64::sin(self.low);
        let mut sin_high = f64::sin(self.high);
        if sin_low > sin_high {
            std::mem::swap(&mut sin_low, &mut sin_high);
        }
        if self.low < f64::consts::PI / 2. && self.high > f64::consts::PI / 2. {
            sin_high = 1.;
        }
        if self.low < (3. / 2.) * f64::consts::PI && self.high > (3. / 2.) * f64::consts::PI {
            sin_low = -1.;
        }
        Interval::new(sin_low, sin_high)
    }

    fn cos(&self) -> Interval {
        let mut cos_low = f64::cos(self.low);
        let mut cos_high = f64::cos(self.high);
        if cos_low > cos_high {
            std::mem::swap(&mut cos_low, &mut cos_high);
        }
        if self.low < f64::consts::PI && self.high > f64::consts::PI {
            cos_low = -1.;
        }
        Interval::new(cos_low, cos_high)
    }

    /// Calculates the zeroes of an interval recursively. Uses Newton's method to find the zeroes once it gets close enough.
    /// The recursion should start at 8. (TODO make a helper function for this?)
    /// It's important to note that it owns the interval given to it for performance reasons, so be sure to copy before giving
    /// the interval if you need to keep the interval.
    #[allow(clippy::float_cmp)] //TODO figure out how to deal with floating point error in a meaningful way
    fn interval_find_zeros(c: (f64, f64, f64, f64, f64), theta: f64, 
            t_int: Interval, zeros: &mut [f64; 4], 
            zero_count: *mut usize, depth: u32) {

        // First, we'll make sure that the function has zeroes. We make an interval of the entire range of values, given 5 derivative terms
        // on the function and an angle.
        let range = Interval::from_pt(c.0) + 
            (Interval::from_pt(c.1) + Interval::from_pt(c.2) * t_int) * 
            (Interval::from_pt(2. * theta) * t_int).cos() +
            (Interval::from_pt(c.3) + Interval::from_pt(c.4) * t_int) * 
            (Interval::from_pt(2. * theta) * t_int).sin();
        if range.low > 0. || range.high < 0. || range.low == range.high { // checks if there are zeroes
            return;
        }
        let mut mid = (t_int.low + t_int.high) * 0.5;
        // we haven't recurred enough times and need to go a layer deeper (don't want to blow up the stack, though)
        // every layer deeper, we split the interval in half
        if depth > 0 {
            Interval::interval_find_zeros(c, theta, Interval::new(t_int.low, mid), zeros, zero_count, depth - 1);
            Interval::interval_find_zeros(c, theta, Interval::new(mid, t_int.high), zeros, zero_count, depth - 1);
        } else { // we use Newton's method - taking the second derivative of the motion function  (nasty).
            for _ in 0..4 {
                let f_newton = c.0 + (c.1 + c.2 * mid) * f64::cos(2. * theta * mid);
                let f_prime_newton = 
                    (c.2 + 2. * (c.3 + c.4 * mid) * theta) *
                    f64::cos(2. * mid * theta) +
                    (c.4 - 2. * (c.1 + c.2 * mid) * theta) *
                    f64::sin(2. * mid * theta);
                if math::eq_f64(f_newton, 0.) || math::eq_f64(f_prime_newton, 0.) { //TODO does this even work?
                    break;
                }
                mid -= -f_newton / f_prime_newton;
            }
            unsafe {
                zeros[*zero_count] = mid;
                *zero_count += 1; //TODO does this need to be unsafe? Probably not.
            }
        }
    }
}

impl Add for Interval {
    type Output = Self;
    fn add(self, other: Interval) -> Interval {
        Interval {
            high: self.high + other.high,
            low: self.low + other.low
        }
    }
}

impl Sub for Interval {
    type Output = Self;
    fn sub(self, other: Interval) -> Interval {
        Interval {
            high: self.high - other.high,
            low: self.low - other.low
        }
    }
}

impl Mul for Interval {
    type Output = Self;
    fn mul(self, other: Interval) -> Interval {
        Interval {
            high: math::max(math::max(self.low * other.low, self.high * other.low),
                            math::max(self.low * other.high, self.high * other.high)),
            low: math::min(math::min(self.low * other.low, self.high * other.low),
                            math::min(self.low * other.high, self.high * other.high))
        }
    }
}

struct DerivativeTerm {
    c: f64,
    x: f64,
    y: f64,
    z: f64
}

impl DerivativeTerm {
    fn new(c: f64, x: f64, y: f64, z: f64) -> DerivativeTerm {
        DerivativeTerm {
            c,
            x,
            y,
            z
        }
    }

    fn default() -> DerivativeTerm {
        DerivativeTerm {
            c: 0.,
            x: 0.,
            y: 0.,
            z: 0.
        }
    }

    fn eval(&self, p: &Point3) -> f64 {
        self.c + self.x * p.x + self.y * p.y + self.z * p.z
    }
}