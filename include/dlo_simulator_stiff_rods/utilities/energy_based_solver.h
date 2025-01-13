/*
 * Author: Burak Aksoy
 * Implements an alternative solver for the DLO simulator. 
 * This solver is based on the minimization of the potential energy of the system.
 */

#ifndef ENERGY_BASED_SOLVER_H
#define ENERGY_BASED_SOLVER_H

#include <ros/ros.h> // For ROS logging
#include "dlo_simulator_stiff_rods/utilities/dlo.h"

#include <memory> // needed for std::shared_ptr
#include <chrono>  // Include this header for timing functions

#include <ceres/ceres.h>


namespace utilities
{
    class EnergyBasedSolver
    {
    public:
        // ----------------------------------------------
        struct DloEnergyModelParams
        {

            Real density = 1.0e+3; // dlo mass per meter cube (kg/m^3)
            Real radius = 0.5e-2; // radius of dlo assuming it's cylndrical (m)

            Real stretch_stiffness = 1.0e+14; 
            Real young_modulus = 3.0e+10; // E (Pa)
            Real torsion_modulus = 1.0e+10; // G (Pa)

            Eigen::Matrix<Real,3,1> gravity = Eigen::Matrix<Real,3,1>(0.0, 0.0, -9.804);
            
            // Eigen::Vector2d log_twist_stiffness_range = Eigen::Vector2d(-4, 0); // lower_bound, upper_bound
            // Eigen::Vector2d log_density_range = Eigen::Vector2d(-4, 1);
        };

        // ----------------------------------------------
        struct OptimizeOptions
        {
            std::string solver = "Ceres"; // Initialize the solver here

            std::string minimizer_type = "trust_region"; // default: trust_region, other options: line_search

            std::string line_search_direction_type = "LBFGS"; // default: LBFGS, other options: BFGS, NCG, STEEPEST_DESCENT

            std::string trust_region_strategy_type = "LEVENBERG_MARQUARDT"; // default: LEVENBERG_MARQUARDT, other options: DOGLEG

            std::string linear_solver_type = "SPARSE_NORMAL_CHOLESKY"; // default: SPARSE_NORMAL_CHOLESKY, other options: DENSE_QR

            double max_num_iterations = 10000; 
            double max_num_line_search_step_size_iterations = 100; // default: 20

            double function_tolerance = 1e-9; // default: 1e-6
            double gradient_tolerance = 1e-10; // default: 1e-10
            double parameter_tolerance = 1e-8; // default: 1e-8
            
            int num_threads = 1; // default: 1
            double initial_trust_region_radius = 1e4; // default: 1e4
            double max_trust_region_radius = 1e16; // default: 1e16

    
            bool print_log = true; // Prints the optimization report to the console with time information
            bool minimizer_progress_to_stdout = true; // default: false, prints the iterations of the minimizer to the console
        };

        // ----------------------------------------------
        // Functions

        // Default constructor
        EnergyBasedSolver(); 

        ~EnergyBasedSolver();

        // Functions
        // Optimize the shape of the DLO given the DLO object
        pbd_object::MeshDLO optimizeShape(pbd_object::Dlo &dlo,
                                            OptimizeOptions &optimize_options,
                                            std::vector<int> &static_particles);

        // Optimize the shape of the DLO given the DLO state and the energy model parameters
        pbd_object::MeshDLO optimizeShape(pbd_object::MeshDLO &dlo_state,
                                        DloEnergyModelParams &energy_model_params,
                                        OptimizeOptions &optimize_options,
                                        std::vector<int> &static_particles);

        // void renderState(pbd_object::MeshDLO &dlo_state, ...);

        void saveDLOStateAsNpy(pbd_object::MeshDLO &dlo_state, 
                            std::string file_name, 
                            bool append_current_time);
        
        void SetCeresOptions(const OptimizeOptions& optimize_options, 
                            ceres::Solver::Options* ceres_options);

    private:
        // Variables
        // pbd_object::Dlo &dlo_;
        

        // Functions
        Eigen::Matrix<Real,3,1> gravityRefPoint(const std::vector<Eigen::Matrix<Real, 3, 1>> &vertices,
                                            const Real &dlo_length);
        
    };

} // namespace utilities

#endif /* !ENERGY_BASED_SOLVER_H */