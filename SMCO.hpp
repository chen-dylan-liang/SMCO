/**
 * Strategic Monte Carlo Optimization (SMCO) - Eigen Optimized Version
 * * Description:
 * A high-performance implementation of SMCO using the Eigen linear algebra library.
 * This version leverages SIMD vectorization for faster vector updates, making it
 * suitable for high-dimensional optimization problems.
 *
 * Dependencies: 
 * - Eigen 3 (http://eigen.tuxfamily.org/)
 * * Compile command example:
 * g++ -I /usr/local/include/eigen3 smco_eigen.cpp -o smco_opt -O3
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <random>
#include <iomanip>
#include <Eigen/Dense> // Requires Eigen library

// ============================================================================
// Configuration
// ============================================================================

struct SMCOParams {
    int max_iterations = 10000;
    double tolerance = 1e-6;
    double step_size_grad = 1e-5;
    bool minimize = false;
    int random_seed = 42;
};

struct SMCOResult {
    Eigen::VectorXd optimal_x;
    double optimal_value;
    int iterations_performed;
    bool converged;
};

// ============================================================================
// SMCO Optimizer Class (Eigen Version)
// ============================================================================

class SMCOptimizer {
private:
    std::mt19937 rng;

public:
    SMCOptimizer() {
        std::random_device rd;
        rng.seed(rd());
    }

    SMCOResult solve(
        std::function<double(const Eigen::VectorXd&)> func,
        const Eigen::VectorXd& lower_bounds,
        const Eigen::VectorXd& upper_bounds,
        SMCOParams params = SMCOParams()
    ) {
        // 1. Validation
        long dim = lower_bounds.size();
        if (upper_bounds.size() != dim) {
            throw std::invalid_argument("Bounds dimension mismatch.");
        }

        rng.seed(params.random_seed);

        // Wrapper to handle Minimization vs Maximization internally
        auto objective = [&](const Eigen::VectorXd& x) {
            return params.minimize ? -func(x) : func(x);
        };

        // 2. Initialization
        // Initialize x randomly within bounds
        Eigen::VectorXd current_x(dim);
        for (int i = 0; i < dim; ++i) {
            std::uniform_real_distribution<double> dist(lower_bounds[i], upper_bounds[i]);
            current_x[i] = dist(rng);
        }

        // S_0 = x_0
        Eigen::VectorXd S = current_x;
        double current_val = objective(current_x);
        bool converged = false;
        int n = 0;

        // Pre-allocate vector Z to avoid re-allocation inside loop
        Eigen::VectorXd Z(dim);

        // 3. Main Loop
        for (n = 0; n < params.max_iterations; ++n) {
            
            // A. Compute Gradient (Finite Difference)
            // Note: We still need a loop here because we perturb one dimension at a time.
            // However, the vector addition inside is now optimized by Eigen.
            Eigen::VectorXd grad(dim);
            for (int j = 0; j < dim; ++j) {
                double original_xj = current_x[j];
                
                // Perturb x_j
                current_x[j] += params.step_size_grad;
                double f_plus = objective(current_x);
                
                // Restore x_j
                current_x[j] = original_xj;

                grad[j] = (f_plus - current_val) / params.step_size_grad;
            }

            // B. Strategic Sampling (The "Two-Armed" Bandit)
            // This part is inherently conditional per dimension, so we loop.
            // Ideally, we could use Eigen::select, but random generation usually requires sequential calls.
            for (int j = 0; j < dim; ++j) {
                double L = lower_bounds[j];
                double U = upper_bounds[j];
                double Mid = (L + U) * 0.5;

                if (grad[j] >= 0) {
                    // Positive gradient -> Sample upper half
                    std::uniform_real_distribution<double> dist(Mid, U);
                    Z[j] = dist(rng);
                } else {
                    // Negative gradient -> Sample lower half
                    std::uniform_real_distribution<double> dist(L, Mid);
                    Z[j] = dist(rng);
                }
            }

            // C. Update Step (Eigen Optimization shines here)
            // Vectorized Addition
            S += Z;

            // Vectorized Scalar Division
            // x_{n+1} = S_{n+1} / (n + 2)
            Eigen::VectorXd next_x = S / (n + 2.0);

            // D. Check Convergence
            double next_val = objective(next_x);
            
            if (std::abs(next_val - current_val) <= params.tolerance) {
                current_x = next_x;
                current_val = next_val;
                converged = true;
                break;
            }

            current_x = next_x;
            current_val = next_val;
        }

        return {current_x, (params.minimize ? -current_val : current_val), n, converged};
    }
};

// ============================================================================
// Examples
// ============================================================================

double rastrigin_eigen(const Eigen::VectorXd& x) {
    // Array-based reduction (fast in Eigen)
    // formula: 10n + sum(x^2 - 10cos(2pi*x))
    double sum = (x.array().square() - 10 * (2 * M_PI * x.array()).cos()).sum();
    return 10 * x.size() + sum;
}

int main() {
    SMCOptimizer optimizer;
    SMCOParams params;

    // Setup Rastrigin Problem (50 Dimensions - where Eigen starts to help)
    int dim = 50;
    std::cout << "=== Running SMCO with Eigen on " << dim << "D Rastrigin Function ===" << std::endl;

    // Vectorized bounds setup
    Eigen::VectorXd lb = Eigen::VectorXd::Constant(dim, -5.12);
    Eigen::VectorXd ub = Eigen::VectorXd::Constant(dim, 5.12);

    params.max_iterations = 5000;
    params.minimize = true;
    params.random_seed = 999;

    // Solve
    auto start = std::chrono::high_resolution_clock::now();
    SMCOResult res = optimizer.solve(rastrigin_eigen, lb, ub, params);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;

    // Output
    std::cout << "Converged: " << (res.converged ? "Yes" : "No") << std::endl;
    std::cout << "Iterations: " << res.iterations_performed << std::endl;
    std::cout << "Optimal Value: " << res.optimal_value << " (Target: 0.0)" << std::endl;
    std::cout << "Time Elapsed: " << elapsed.count() << "s" << std::endl;

    // Show first 5 coordinates
    std::cout << "First 5 coords: " << res.optimal_x.head(5).transpose() << std::endl;

    return 0;
}