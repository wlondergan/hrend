use core::f32;
use std::random::Random;
use crate::math::{next_down, ONE_MINUS_EPSILON};

use super::math::square;



pub trait Integrator {

    fn new() -> Self;

}

pub trait RayIntegrator: Integrator {

}

struct RandomWalkIntegrator {

}

impl Integrator for RandomWalkIntegrator {

    fn new() -> Self {
        RandomWalkIntegrator {}
    }
}

impl RandomWalkIntegrator {

}

/**
 * Applies the balance heuristic to the sampling technique `f` with `g` being the other sampling technique.
 * Used in Monte Carlo integration.
 * See: pbrt, 2.14
 */
fn balance_heuristic(nf: i32, f_pdf: f32, ng: i32, g_pdf: f32) -> f32 {
    (nf as f32 * f_pdf) / (nf as f32 * f_pdf + ng as f32 * g_pdf)
}

/**
 * Applies the power heuristic to the sampling technique `f` with `g` being the other sampling technique.
 * Uses `beta=2`, since apparently that works fine for this weighting technique.
 * See: pbrt, 2.15
 */
fn power_heuristic(nf: i32, f_pdf: f32, ng: i32, g_pdf: f32) -> f32 {
    let f = nf as f32 * f_pdf;
    let g = ng as f32 * g_pdf;
    square(f) / (square(f) + square(g))
}

/**
  Applies the inversion sampling technique to the discrete PMF given by `weights`.

 `u` is the uniform sample used for the sample, `ret_pmf` provides the value of the PMF in the sample in return value 2, and `ret_u_remap`
 provides a new uniform random sample derived from `u` in return value 3.
 */
fn sample_discrete(weights: &[f32], u: f32, ret_pmf: bool, ret_u_remap: bool) -> (i32, f32, f32) {
    if weights.is_empty() {
        return (-1, 0.0, 0.0);
    }
    let sum_weights: f32 = weights.iter().sum();
    let mut uprime = u * sum_weights;

    if uprime == sum_weights {
        uprime = next_down(uprime);
    }

    let mut offset = 0;
    let mut sum = 0.0;
    while sum + weights[offset] <= uprime {
        sum += weights[offset];
        offset += 1;
    }

    let mut pmf = 0.0;
    let mut u_remap = 0.0;

    if ret_pmf {
        pmf = weights[offset] / sum_weights;
    }
    if ret_u_remap {
        u_remap = f32::min((uprime - sum) / weights[offset], ONE_MINUS_EPSILON);
    }

    (offset, pmf, u_remap)

}

