/**
 * Strategic Monte Carlo Optimization (SMCO)
 * C++ Implementation (Eigen Optimized)
 * * Based on the R implementation:
 * https://github.com/wayne-y-gao/SMCO/blob/main/SMCO.R
 * * LOGIC SUMMARY:
 * 1. Optimization via Strategic Law of Large Numbers.
 * 2. Strategic Draw 'Z': Deterministic selection of extended bounds based on gradient sign.
 * - If Grad > 0: Z = Upper_Bound + Buffer
 * - If Grad < 0: Z = Lower_Bound - Buffer
 * 3. Hierarchy: Single -> Refine (Two-stage) -> Boost (iterative).
 * * Dependencies: Eigen 3
 */

#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <limits>
#include <Eigen/Dense>

// ============================================================================
// Configuration Structures
// ============================================================================

struct SMCOParams {
    // Multi-start parameters
    int n_starts = 100;

    // Iteration control
    int iter_max = 100;
    int iter_nstart = 1;        // n_0 in the paper
    int iter_boost = 0;         // k in the paper (virtual iteration offset)

    // Bounds & Buffers
    double bounds_buffer = 0.05; // 5% expansion of bounds for Z
    bool buffer_rand = false;    // If true, multiplies buffer by random U(-1,1)

    // Convergence & Logic
    double tol_conv = 1e-6;
    bool refine_search = true;
    double refine_ratio = 0.5;   // Split iter_max between initial and refine stages
    bool use_runmax = true;      // Keep track of best point seen during gradient steps
    bool minimize = false;       // Internal logic is Maximization. Set true to minimize f(x).

    int random_seed = 123;
};

struct SingleResult {
    Eigen::VectorXd x_optimal; // The "Strategic Mean" X_n
    double f_optimal;

    Eigen::VectorXd x_runmax;  // The greedy best point seen
    double f_runmax;

    int iterations;
};

struct SMCOResult {
    Eigen::VectorXd optimal_x;
    double optimal_value;
    int iterations_performed;
    bool converged;
};

// ============================================================================
// SMCO Optimizer Class
// ============================================================================

class SMCOptimizer {
private:
    std::mt19937 rng;

    // Helper: Check bounds and clamp if necessary
    // Returns clamped vector
    Eigen::VectorXd check_bounds(const Eigen::VectorXd& x,
                                 const Eigen::VectorXd& lb,
                                 const Eigen::VectorXd& ub) {
        return x.cwiseMax(lb).cwiseMin(ub);
    }

public:
    SMCOptimizer() {
        std::random_device rd;
        rng.seed(rd());
    }

    // ------------------------------------------------------------------------
    // Core Routine: SMCO_single
    // Corresponds to R function: SMCO_single
    // ------------------------------------------------------------------------
    SingleResult solve_single(
            const std::function<double(const Eigen::VectorXd&)>& f,
            const Eigen::VectorXd& lb,
            const Eigen::VectorXd& ub,
            const Eigen::VectorXd& start_point,
            double current_bounds_buffer,
            const SMCOParams& params,
            int current_iter_max,
            int current_iter_boost
    ) {
        int d = lb.size();
        Eigen::VectorXd bounds_diff = ub - lb;

        // Initialize State
        Eigen::VectorXd x_current = start_point;
        double f_current = f(x_current);

        // Running Max (Greedy tracking)
        Eigen::VectorXd x_runmax = x_current;
        double f_runmax = f_current;

        // Counters
        int n_start = current_iter_boost + params.iter_nstart;
        int n_end = n_start + current_iter_max;

        // Initial Sum S_n = x_start * n
        // Note: R code does S = start_point * n_boost_1
        Eigen::VectorXd S = start_point * (double)n_start;

        // Main Loop
        int n = n_start;
        for (; n < n_end; ++n) {

            // 1. Compute Adaptive Step Size h
            // R: h_step <- bounds_diff / (n + 1)
            Eigen::VectorXd h_step = bounds_diff / (double)(n + 1);

            // 2. Compute Partial Gradient Signs (Finite Difference)
            // R: compute_partial_signs(..., partial_option="center")
            // We implement "center" (two-sided) strategy here as it's robust.

            Eigen::VectorXi signs(d); // 1 for Positive, 0 for Negative

            Eigen::VectorXd x_plus = x_current;
            Eigen::VectorXd x_minus = x_current;

            for (int j = 0; j < d; ++j) {
                // Constrained perturbation
                double h = h_step[j];
                double val_plus = std::min(x_current[j] + h, ub[j]);
                double val_minus = std::max(x_current[j] - h, lb[j]);

                double original_val = x_current[j];

                // Eval f(x+)
                x_current[j] = val_plus;
                double f_plus = f(x_current);

                // Eval f(x-)
                x_current[j] = val_minus;
                double f_minus = f(x_current);

                // Reset
                x_current[j] = original_val;

                // Determine Sign
                bool sign_positive = (f_plus > f_minus);
                signs[j] = sign_positive ? 1 : 0;

                // Update Greedy Runmax
                if (params.use_runmax) {
                    if (sign_positive && f_plus > f_runmax) {
                        f_runmax = f_plus;
                        x_runmax = x_current; // Base vector
                        x_runmax[j] = val_plus;
                    } else if (!sign_positive && f_minus > f_runmax) {
                        f_runmax = f_minus;
                        x_runmax = x_current; // Base vector
                        x_runmax[j] = val_minus;
                    }
                }
            }

            // 3. Calculate Strategic Draw 'Z'
            // R: Z <- signs * bounds_upper_out + (1-signs) * bounds_lower_out
            // Logic: If ascent (sign=1), pull towards Upper Bound. Else Lower Bound.

            // Calculate Expanded Bounds
            Eigen::VectorXd bound_out_amt = bounds_diff * current_bounds_buffer;

            if (params.buffer_rand) {
                // If random buffer enabled (R: runif(d, -1, 1))
                std::uniform_real_distribution<double> dist(-1.0, 1.0);
                for(int k=0; k<d; ++k) bound_out_amt[k] *= dist(rng);
            }

            Eigen::VectorXd ub_out = ub + bound_out_amt;
            Eigen::VectorXd lb_out = lb - bound_out_amt;

            Eigen::VectorXd Z(d);
            for (int j = 0; j < d; ++j) {
                Z[j] = (signs[j] == 1) ? ub_out[j] : lb_out[j];
            }

            // 4. Update Sum and Average
            // S_{n+1} = S_n + Z
            // X_{n+1} = S_{n+1} / (n+1)
            S += Z;
            Eigen::VectorXd x_next = S / (double)(n + 1);
            double f_next = f(x_next);

            // Update Runmax with the new average
            if (params.use_runmax && f_next > f_runmax) {
                f_runmax = f_next;
                x_runmax = x_next;
            }

            // 5. Convergence Check
            // R: checks convergence after 50% of iterations
            double diff = std::abs(f_next - f_current);
            if (n >= (n_start + current_iter_max/2) && diff < params.tol_conv) {
                x_current = x_next;
                f_current = f_next;
                break;
            }

            x_current = x_next;
            f_current = f_next;
        }

        return {x_current, f_current, x_runmax, f_runmax, n - current_iter_boost};
    }

    // ------------------------------------------------------------------------
    // Refinement Wrapper: SMCO_single_refine
    // ------------------------------------------------------------------------
    SingleResult solve_refine(
            const std::function<double(const Eigen::VectorXd&)>& f,
            const Eigen::VectorXd& lb,
            const Eigen::VectorXd& ub,
            const Eigen::VectorXd& start_point,
            double bounds_buffer,
            const SMCOParams& params,
            int iter_max,
            int iter_boost
    ) {
        // 1. Initial Search
        // If refining, split iterations. Else use full.
        int iter_initial = params.refine_search ?
                           std::round(iter_max * (1.0 - params.refine_ratio)) :
                           iter_max;

        SingleResult res = solve_single(f, lb, ub, start_point, bounds_buffer,
                                        params, iter_initial, iter_boost);

        // Clamp results to valid bounds
        res.x_optimal = check_bounds(res.x_optimal, lb, ub);
        res.f_optimal = f(res.x_optimal);
        if (params.use_runmax) {
            res.x_runmax = check_bounds(res.x_runmax, lb, ub);
            res.f_runmax = f(res.x_runmax);
        }

        // 2. Refined Search (Optional)
        if (params.refine_search) {
            int iter_refine = std::round(iter_max * params.refine_ratio);

            // Pick best point to start refinement
            Eigen::VectorXd start_refine = res.x_optimal;
            if (params.use_runmax && res.f_runmax > res.f_optimal) {
                start_refine = res.x_runmax;
            }

            // R: iter_boost_refine = iter_boost + 1000
            // R: bounds_buffer = 0 (Tighten bounds for refinement)
            int boost_refine = iter_boost + 1000;

            SingleResult res_refine = solve_single(f, lb, ub, start_refine,
                                                   0.0, // Zero buffer
                                                   params, iter_refine, boost_refine);

            // Merge Refined Results
            if (params.use_runmax) {
                // Check if refined runmax is better
                Eigen::VectorXd r_rm_in = check_bounds(res_refine.x_runmax, lb, ub);
                double f_r_rm_in = f(r_rm_in);

                if (f_r_rm_in > res_refine.f_optimal) {
                    res_refine.f_optimal = f_r_rm_in;
                    res_refine.x_optimal = r_rm_in;
                }
            }
            return res_refine;
        }

        // If no refinement, just return cleaned result
        // Logic to swap optimal with runmax if runmax is better
        if (params.use_runmax && res.f_runmax > res.f_optimal) {
            res.f_optimal = res.f_runmax;
            res.x_optimal = res.x_runmax;
        }

        return res;
    }

    // ------------------------------------------------------------------------
    // Boost Wrapper: SMCO_single_boost
    // ------------------------------------------------------------------------
    SingleResult solve_boost(
            const std::function<double(const Eigen::VectorXd&)>& f,
            const Eigen::VectorXd& lb,
            const Eigen::VectorXd& ub,
            const Eigen::VectorXd& start_point,
            const SMCOParams& params
    ) {
        // 1. Run regular search (Boost = 0)
        SingleResult res = solve_refine(f, lb, ub, start_point, params.bounds_buffer,
                                        params, params.iter_max, 0);

        // 2. Run additional boosted search if configured
        if (params.iter_boost > 0) {
            SingleResult res_boost = solve_refine(f, lb, ub, start_point, params.bounds_buffer,
                                                  params, params.iter_max, params.iter_boost);

            if (res_boost.f_optimal > res.f_optimal) {
                res = res_boost;
            }
        }
        return res;
    }

    // ------------------------------------------------------------------------
    // Public Entry Point: solve (Multi-Start)
    // ------------------------------------------------------------------------
    SMCOResult solve(
            std::function<double(const Eigen::VectorXd&)> user_func,
            const Eigen::VectorXd& lb,
            const Eigen::VectorXd& ub,
            SMCOParams params = SMCOParams()
    ) {
        // Handle Minimization Wrapper
        auto objective = [&](const Eigen::VectorXd& x) {
            return params.minimize ? -user_func(x) : user_func(x);
        };

        rng.seed(params.random_seed);
        int d = lb.size();

        SingleResult global_best;
        global_best.f_optimal = -std::numeric_limits<double>::infinity();

        int total_iter = 0;

        // Multi-start Loop
        // R uses Sobol. We use Uniform Random here for simplicity in a single file,
        // but iterating this matches the "apply" logic in R.
        for (int i = 0; i < params.n_starts; ++i) {

            // Generate Start Point (Sobol equivalent placeholder)
            Eigen::VectorXd start_point(d);
            for (int j = 0; j < d; ++j) {
                std::uniform_real_distribution<double> dist(lb[j], ub[j]);
                start_point[j] = dist(rng);
            }

            // Execute SMCO Boosted for this point
            SingleResult res = solve_boost(objective, lb, ub, start_point, params);

            total_iter += res.iterations;

            if (res.f_optimal > global_best.f_optimal) {
                global_best = res;
            }
        }

        // Final cleanup
        SMCOResult result;
        result.optimal_x = global_best.x_optimal;
        // Revert sign if minimizing
        result.optimal_value = params.minimize ? -global_best.f_optimal : global_best.f_optimal;
        result.iterations_performed = total_iter;
        result.converged = true; // Simplified flag

        return result;
    }
};

// ============================================================================
// Example Usage
// ============================================================================

double rastrigin_eigen(const Eigen::VectorXd& x) {
    double sum = (x.array().square() - 10 * (2 * M_PI * x.array()).cos()).sum();
    return 10 * x.size() + sum;
}

int main() {
    SMCOptimizer optimizer;
    SMCOParams params;

    // Setup Rastrigin Problem (50 Dimensions)
    // R Default: n_starts = 100, iter_max = 200
    int dim = 2;
    std::cout << "=== SMCO (Deterministic Strategic Arms) on " << dim << "D Rastrigin ===" << std::endl;

    Eigen::VectorXd lb = Eigen::VectorXd::Constant(dim, -5.12);
    Eigen::VectorXd ub = Eigen::VectorXd::Constant(dim, 5.12);

    params.n_starts = 1;      // Reduced for quick demo
    params.iter_max = 500;
    params.minimize = true;    // Rastrigin is minimization
    params.random_seed = 24;

    // R defaults:
    params.bounds_buffer = 0.05;
    params.refine_search = true;
    params.use_runmax = true;

    auto start = std::chrono::high_resolution_clock::now();
    SMCOResult res = optimizer.solve(rastrigin_eigen, lb, ub, params);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Best Value: " << res.optimal_value << " (Target: 0.0)" << std::endl;
    std::cout << "Total Iterations (Across all starts): " << res.iterations_performed << std::endl;
    std::cout << "Time: " << elapsed.count() << "s" << std::endl;
   // std::cout << "First 5 X: " << res.optimal_x.head(5).transpose() << std::endl;

    return 0;
}
