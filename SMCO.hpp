/**
 * Strategic Monte Carlo Optimization (SMCO)
 * C++ Implementation
 *
 * Based on the paper: "Optimization via Strategic Law of Large Numbers" (arXiv:2412.05604)
 * and the GitHub repository: wayne-y-gao/SMCO
 *
 * Author: Gemini
 * Date: 2025
 *
 * DESCRIPTION:
 * This class implements the SMCO algorithm, which optimizes a function by treating
 * the optimization variable as a running average of "strategic" random draws.
 * At each step, it estimates the gradient sign and chooses to draw from a "High Mean"
 * or "Low Mean" distribution (Arms X and Y) to steer the average towards the optimum.
 *
 * FEATURES:
 * - Gradient-free interface (uses internal finite difference if gradient not provided).
 * - Support for arbitrary dimensions.
 * - Customizable bounds per dimension.
 * - "Uniform Split" default strategy for arms.
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <random>
#include <numeric>
#include <iomanip>
#include <string>

// ============================================================================
// SMCO Configuration & Helper Structures
// ============================================================================

struct SMCOParams {
    int max_iterations = 10000;     // Maximum number of iterations (n in paper)
    double tolerance = 1e-6;        // Convergence tolerance (epsilon)
    double step_size_grad = 1e-5;   // h for finite difference gradient estimation
    bool minimize = false;          // If true, minimizes f(x). If false, maximizes.
    int random_seed = 42;           // Seed for reproducibility
};

struct SMCOResult {
    std::vector<double> optimal_x;  // The best x found
    double optimal_value;           // f(optimal_x)
    int iterations_performed;       // Number of steps taken
    bool converged;                 // True if |f(x_new) - f(x_old)| <= tolerance
};

// ============================================================================
// SMCO Optimizer Class
// ============================================================================

class SMCOptimizer {
private:
    std::mt19937 rng;

public:
    SMCOptimizer() {
        std::random_device rd;
        rng.seed(rd());
    }

    void set_seed(int seed) {
        rng.seed(seed);
    }

    /**
     * Main Optimization Function
     *
     * @param func The objective function f(x) -> double.
     * @param lower_bounds Vector of lower bounds for each dimension.
     * @param upper_bounds Vector of upper bounds for each dimension.
     * @param params Configuration parameters (optional).
     * @return SMCOResult containing the best point and value.
     */
    SMCOResult solve(
        std::function<double(const std::vector<double>&)> func,
        const std::vector<double>& lower_bounds,
        const std::vector<double>& upper_bounds,
        SMCOParams params = SMCOParams()
    ) {
        // 1. Input Validation
        size_t dim = lower_bounds.size();
        if (upper_bounds.size() != dim) {
            throw std::invalid_argument("Lower and upper bounds must have the same dimension.");
        }

        // Reset RNG if seed provided in params
        rng.seed(params.random_seed);

        // 2. Internal Objective Wrapper
        // SMCO is naturally a Maximization algorithm (moves towards positive gradient).
        // If user wants to Minimize, we optimize -f(x).
        auto objective = [&](const std::vector<double>& x) {
            return params.minimize ? -func(x) : func(x);
        };

        // 3. Initialization
        // In SMCO, S_0 is the initial sum. \hat{x}_0 is the initial estimate.
        // We initialize randomly within bounds.
        std::vector<double> current_x(dim);
        std::vector<double> S(dim); // Running Sum

        for (size_t j = 0; j < dim; ++j) {
            std::uniform_real_distribution<double> dist(lower_bounds[j], upper_bounds[j]);
            current_x[j] = dist(rng);
            S[j] = current_x[j]; // S_0 = x_0
        }

        double current_val = objective(current_x);
        bool converged = false;
        int n = 0;

        // 4. Main Loop (The Strategic Law of Large Numbers process)
        for (n = 0; n < params.max_iterations; ++n) {
            std::vector<double> Z(dim); // The "Strategic Draw" for this step

            // Loop over dimensions
            for (size_t j = 0; j < dim; ++j) {
                // A. Estimate Partial Derivative (Gradient)
                // We use finite difference: (f(x + h*e_j) - f(x)) / h
                double grad_j = compute_partial_derivative(objective, current_x, current_val, j, params.step_size_grad);

                // B. Strategic Choice (The "Two-Armed" Decision)
                // If grad >= 0 (Function increases this way), we pick the "High Mean" arm.
                // If grad < 0 (Function decreases this way), we pick the "Low Mean" arm.
                
                double L = lower_bounds[j];
                double U = upper_bounds[j];
                double Mid = (L + U) / 2.0;

                if (grad_j >= 0) {
                    // Arm X (High): Sample from upper half [Mid, U]
                    // This pulls the average UP.
                    std::uniform_real_distribution<double> dist_high(Mid, U);
                    Z[j] = dist_high(rng);
                } else {
                    // Arm Y (Low): Sample from lower half [L, Mid]
                    // This pulls the average DOWN.
                    std::uniform_real_distribution<double> dist_low(L, Mid);
                    Z[j] = dist_low(rng);
                }
            }

            // C. Update Running Sum and Estimate
            // Formula: S_{n+1} = S_n + Z_{n+1}
            //          x_{n+1} = S_{n+1} / (n + 2)  <-- Note: n starts at 0, so we divide by count (n+1 previous + 1 new)
            // Note on indexing: If S_0 is the 1st sample, then at n=0 we add Z (2nd sample). Total count is n+2.
            
            std::vector<double> next_x(dim);
            for (size_t j = 0; j < dim; ++j) {
                S[j] += Z[j];
                next_x[j] = S[j] / (n + 2.0); 
            }

            // D. Check Convergence
            // |f(next_x) - f(current_x)| <= tolerance
            double next_val = objective(next_x);
            
            if (std::abs(next_val - current_val) <= params.tolerance) {
                converged = true;
                current_x = next_x;
                current_val = next_val;
                break; 
            }

            // Update for next iteration
            current_x = next_x;
            current_val = next_val;
        }

        // 5. Construct Result
        SMCOResult result;
        result.optimal_x = current_x;
        // If we were minimizing, we negated the function. Negate back for the true value.
        result.optimal_value = params.minimize ? -current_val : current_val;
        result.iterations_performed = n;
        result.converged = converged;

        return result;
    }

private:
    /**
     * Helper: Compute Partial Derivative via Finite Difference
     */
    double compute_partial_derivative(
        const std::function<double(const std::vector<double>&)>& func,
        const std::vector<double>& x,
        double current_f_val,
        size_t dim_index,
        double h
    ) {
        std::vector<double> x_plus = x;
        x_plus[dim_index] += h;
        
        double f_plus = func(x_plus);
        
        return (f_plus - current_f_val) / h;
    }
};

// ============================================================================
// Example Usage / Demo
// ============================================================================

// Test Function: Rastrigin Function (Standard Optimization Benchmark)
// Usually defined on [-5.12, 5.12]. Global min at (0,0,...,0) with value 0.
// It has many local minima, making it hard for standard Gradient Descent.
double rastrigin(const std::vector<double>& x) {
    double sum = 0.0;
    double A = 10.0;
    for (double val : x) {
        sum += (val * val - A * std::cos(2 * M_PI * val));
    }
    return A * x.size() + sum;
}

// Test Function: Sphere Function (Simple convex)
double sphere(const std::vector<double>& x) {
    double sum = 0.0;
    for (double val : x) sum += val * val;
    return sum;
}

void print_vector(const std::string& label, const std::vector<double>& v) {
    std::cout << label << " [";
    for (size_t i = 0; i < v.size(); ++i) {
        std::cout << std::fixed << std::setprecision(4) << v[i] << (i < v.size() - 1 ? ", " : "");
    }
    std::cout << "]" << std::endl;
}

int main() {
    SMCOptimizer optimizer;
    SMCOParams params;

    // ---------------------------------------------------------
    // Experiment 1: Sphere Function (Easy)
    // Goal: Minimize Sphere function in 5 Dimensions
    // ---------------------------------------------------------
    std::cout << "=== Experiment 1: Sphere Function (5D) ===" << std::endl;
    
    int dim = 5;
    std::vector<double> lb(dim, -5.0);
    std::vector<double> ub(dim, 5.0);

    params.max_iterations = 2000;
    params.minimize = true; // Sphere has a minimum at 0
    params.random_seed = 123;

    SMCOResult res1 = optimizer.solve(sphere, lb, ub, params);

    print_vector("Optimal X found:", res1.optimal_x);
    std::cout << "Optimal Value: " << res1.optimal_value << std::endl;
    std::cout << "Iterations: " << res1.iterations_performed << std::endl;
    std::cout << "Converged: " << (res1.converged ? "Yes" : "No") << std::endl;
    std::cout << "------------------------------------------\n" << std::endl;


    // ---------------------------------------------------------
    // Experiment 2: Rastrigin Function (Hard, Multi-modal)
    // Goal: Minimize Rastrigin in 2 Dimensions
    // Global Minimum is at [0, 0] -> Value 0
    // ---------------------------------------------------------
    std::cout << "=== Experiment 2: Rastrigin Function (2D) ===" << std::endl;
    
    dim = 2;
    lb = {-5.12, -5.12};
    ub = {5.12, 5.12};

    params.max_iterations = 5000; // Needs more iterations to average out local traps
    params.minimize = true;
    params.random_seed = 999;

    SMCOResult res2 = optimizer.solve(rastrigin, lb, ub, params);

    print_vector("Optimal X found:", res2.optimal_x);
    std::cout << "Optimal Value: " << res2.optimal_value << std::endl;
    std::cout << "Iterations: " << res2.iterations_performed << std::endl;
    
    // Check distance to true global optimum (0,0)
    double dist = std::sqrt(std::pow(res2.optimal_x[0], 2) + std::pow(res2.optimal_x[1], 2));
    std::cout << "Distance to Global Min (0,0): " << dist << std::endl; 

    return 0;
}