//! Nonlinear least squares on a few unknowns: Levenberg–Marquardt steps
//! on the normal equations, accumulated on the stack.
//!
//! The model is a state that moves by a step of `N` numbers. The step is
//! local: a state may hold a unit direction moved along the two directions
//! across it, and rebuild that basis as it moves, so no parametrization of
//! the sphere of directions has a pole on the way.

use crate::math::numeric::solve::solve_square;

/// The most accepted steps; a fit from a fair start takes a handful.
const MAX_STEPS: usize = 100;

/// The most damping increases tried for one step before the state counts
/// as settled: each multiplies the damping by ten.
const MAX_DAMPING_TRIALS: usize = 16;

/// A step whose every entry is below this, taken while the damping is
/// light, ends the iteration: the model's unknowns are near unit size, so
/// the next step would move them below rounding.
const SETTLED_STEP: f64 = 1e-12;

/// The state that minimizes the sum of the squared residuals, from
/// `start`, by Levenberg–Marquardt steps.
///
/// `residuals(state, index)` gives the `M` residual rows of item `index`,
/// for `index` in `0..count`: each a residual and its derivatives along
/// the `N` entries of a step. `advance(state, step)` is the state moved by
/// `step`. Never fails: it returns the last state that lowered the sum,
/// once a step settles, no damping lowers the sum, or the steps run out.
#[allow(clippy::needless_range_loop)] // index arithmetic on a symmetric matrix
pub(crate) fn levenberg_marquardt<S: Copy, const N: usize, const M: usize>(
    start: S,
    count: usize,
    residuals: impl Fn(&S, usize) -> [(f64, [f64; N]); M],
    advance: impl Fn(&S, [f64; N]) -> S,
) -> S {
    // The sum of squares, `JᵀJ` and `Jᵀr` at a state.
    let normal_equations = |state: &S| {
        let mut sum = 0.0;
        let mut product = [[0.0; N]; N];
        let mut gradient = [0.0; N];
        for index in 0..count {
            for (residual, derivatives) in residuals(state, index) {
                sum += residual * residual;
                for i in 0..N {
                    gradient[i] += derivatives[i] * residual;
                    for j in i..N {
                        product[i][j] += derivatives[i] * derivatives[j];
                    }
                }
            }
        }
        for i in 0..N {
            for j in 0..i {
                product[i][j] = product[j][i];
            }
        }
        (sum, product, gradient)
    };
    let mut state = start;
    let (mut sum, mut product, mut gradient) = normal_equations(&state);
    let mut damping = 1e-3;
    for _ in 0..MAX_STEPS {
        let largest_diagonal = (0..N).map(|i| product[i][i]).fold(0.0, f64::max);
        let mut settled = true;
        for _ in 0..MAX_DAMPING_TRIALS {
            let mut damped = product;
            for (i, row) in damped.iter_mut().enumerate() {
                row[i] += damping * (product[i][i] + f64::EPSILON * largest_diagonal);
            }
            let Some(step) = solve_square(damped, gradient.map(|g| -g), 0.0) else {
                damping *= 10.0;
                continue;
            };
            let trial = advance(&state, step);
            let (trial_sum, trial_product, trial_gradient) = normal_equations(&trial);
            // Near the minimum the sum is flat to rounding long before the
            // state settles; there a step that shrinks the gradient is
            // progress too. A NaN sum compares false, so it is never
            // accepted.
            let flat = trial_sum <= sum * (1.0 + 8.0 * f64::EPSILON)
                && squared_norm(trial_gradient) < squared_norm(gradient);
            if trial_sum < sum || flat {
                state = trial;
                (sum, product, gradient) = (trial_sum, trial_product, trial_gradient);
                let light = damping <= 1.0;
                damping = (damping * 0.1).max(1e-15);
                settled = light && step.iter().all(|s| s.abs() <= SETTLED_STEP);
                break;
            }
            damping *= 10.0;
        }
        if settled {
            break;
        }
    }
    state
}

/// The sum of the squares of `vector`'s entries.
fn squared_norm<const N: usize>(vector: [f64; N]) -> f64 {
    vector.iter().map(|v| v * v).sum()
}

#[cfg(test)]
mod tests;
