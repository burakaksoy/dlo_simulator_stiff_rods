/*
 * Author: Burak Aksoy
 * Implements an alternative solver for the DLO simulator. 
 * This solver is based on the minimization of the potential energy of the system.
 */

#include "dlo_simulator_stiff_rods/utilities/energy_based_solver.h"
#include "dlo_simulator_stiff_rods/utilities/energy_based_model_ceres.h"

using namespace utilities;


// Default constructor
EnergyBasedSolver::EnergyBasedSolver(){
	ROS_INFO("EnergyBasedSolver constructor called.");
}

// Destructor
EnergyBasedSolver::~EnergyBasedSolver(){
}

pbd_object::MeshDLO EnergyBasedSolver::optimizeShape(pbd_object::Dlo &dlo,						
													OptimizeOptions &optimize_options,
													std::vector<int> &static_particles)
{
	// Creates dlo state and energy model parameters from the passed dlo object
	// Then calls the other optimizeShape function

	ROS_INFO("EnergyBasedSolver::optimizeShape() with a passed DLO object is called.");

	// Extract the mesh from the dlo object
	pbd_object::MeshDLO dlo_state;
	dlo_state.name = "current_dlo_state"; // this is a generic name, it can be changed later if needed 
	dlo_state.vertices = *(dlo.getPosPtr()); // std::vector<Eigen::Matrix<Real,3,1>>, dereference the pointer, std::vector creates a full copy 
	dlo_state.quaternions = *(dlo.getOriPtr()); // std::vector<Eigen::Quaternion<Real>>
	dlo_state.segment_lengths = *(dlo.getSegmentLengthsPtr()); // std::vector<Real>

	// Create the energy model parameters from the dlo object as well
	DloEnergyModelParams energy_model_params;
	energy_model_params.density = dlo.getDensity();
	energy_model_params.radius = dlo.getRadius();
	energy_model_params.stretch_stiffness = dlo.getZeroStretchStiffness();
	energy_model_params.young_modulus = dlo.getYoungModulus();
	energy_model_params.torsion_modulus = dlo.getTorsionModulus();
	energy_model_params.gravity = dlo.getGravity();

	return optimizeShape(dlo_state, 
						 energy_model_params,
						 optimize_options,
						 static_particles);
}


// pbd_object::MeshDLO EnergyBasedSolver::optimizeShape(std::vector<Eigen::Matrix<Real,3,1>> )
// {
// 	// Assume that the orientations are not provided except the static particles
// 	// as would be the case, for example, when capturing the DLO shape from the real world using a rgbd camera
// 	// Since the orientations are not provided, the orientations are assumed to be the identity quaternion to initialize the optimization
// 	// TODO


// 	pbd_object::MeshDLO dlo_state;

// 	return optimizeShape(dlo_state, 
// 						 energy_model_params,
// 						 optimize_options,
// 						 static_particles);
// }

pbd_object::MeshDLO EnergyBasedSolver::optimizeShape(pbd_object::MeshDLO &dlo_state,
										DloEnergyModelParams &energy_model_params,
										OptimizeOptions &optimize_options,
										std::vector<int> &static_particles){
	ROS_INFO("EnergyBasedSolver::optimizeShape() called.");

	// Measure time for optimizeShapeDerm function
    auto start_optimize = std::chrono::high_resolution_clock::now();

	// --------------------------------------------------------------------------
	// Create the energy minimization problem
    ceres::Problem problem;

	ceres::LossFunction* loss_function = nullptr;
    ceres::Manifold* quaternion_manifold = new ceres::EigenQuaternionManifold;
	
	int num_particles = dlo_state.vertices.size();
	int num_quaternions = dlo_state.quaternions.size(); 
	
	// -----------------------------------------
	// Create bending and twisting constraints

	// Set the zero-stretch bending twisting constraints
	for (int i = 0; i < (num_quaternions - 1); i++)
	{
		// Segment id pair of each constraint
		const int id0 = i;
		const int id1 = i + 1;
		
		// compute average segment length
		const Real average_segment_length = 0.5 * (dlo_state.segment_lengths[id0] + dlo_state.segment_lengths[id1]);

		// compute rest Darboux vector 
		// if initially the dlo is not straight, use this:
		// Eigen::Quaternion<Real> rest_darboux_vect = dlo_state.quaternions[id0].conjugate() * dlo_state.quaternions[id1]; 
		// if initially the dlo is straight, use this
		Eigen::Quaternion<Real> rest_darboux_vect = Eigen::Quaternion<Real>::Identity(); 
		const Eigen::Matrix<Real,3,1> rest_darboux_vect_vec = (2.0 / average_segment_length) * rest_darboux_vect.vec();

		// locally each constraint position is along the z axis at the half length of the segment
		const Eigen::Matrix<Real, 3, 1> connector0(0, 0, 0.5 * dlo_state.segment_lengths[id0]);
		const Eigen::Matrix<Real, 3, 1> connector1(0, 0, -0.5 * dlo_state.segment_lengths[id1]);

		// Set sqrt stiffness matrix
		Eigen::Matrix<double, 6, 6> sqrt_stiffness_mat = Eigen::Matrix<double, 6, 6>::Identity();

		Real E = energy_model_params.young_modulus;
		Real G = energy_model_params.torsion_modulus;
		Real r = energy_model_params.radius;
		// Real area(static_cast<Real>(M_PI) * std::pow(r, static_cast<Real>(2.0)));
		Real secondMomentOfArea(static_cast<Real>(M_PI_4) * std::pow(r, static_cast<Real>(4.0)));
		
		Real stretchStiffness = energy_model_params.stretch_stiffness;
		Real bendingStiffness(E * secondMomentOfArea);
		Real torsionStiffness(static_cast<Real>(2.0) * G * secondMomentOfArea);

		// Real sqrt_stretch_stiffness = std::sqrt(stretchStiffness * (static_cast<Real>(1.0)/average_segment_length));
		Real sqrt_stretch_stiffness = std::sqrt(stretchStiffness);
		Real sqrt_bending_stiffness = std::sqrt(bendingStiffness * (static_cast<Real>(1.0)/average_segment_length));
		Real sqrt_torsion_stiffness = std::sqrt(torsionStiffness * (static_cast<Real>(1.0)/average_segment_length));
		// Above 1/average_segment_length multiplication is to remove the dependency of the stiffness to the segment length.

		// Set the diagonal elements of the stiffness matrix
		sqrt_stiffness_mat(0, 0) = sqrt_stretch_stiffness;
		sqrt_stiffness_mat(1, 1) = sqrt_stretch_stiffness;
		sqrt_stiffness_mat(2, 2) = sqrt_stretch_stiffness;
		sqrt_stiffness_mat(3, 3) = sqrt_bending_stiffness;
		sqrt_stiffness_mat(4, 4) = sqrt_bending_stiffness;
		sqrt_stiffness_mat(5, 5) = sqrt_torsion_stiffness;

		// ROS info for the stiffness matrix with squared values
		if (i == 0)
		{
			ROS_INFO("sqrt_stiffness_mat: ");
			for (int k = 0; k < 6; k++)
			{
				ROS_INFO("%f ", sqrt_stiffness_mat(k, k)*sqrt_stiffness_mat(k, k));
			}
		}

		// Create the cost function
		ceres::CostFunction* stretch_bend_twist_cost = StretchBendTwistCost::Create(average_segment_length,
																					connector0,
																					connector1,
																					rest_darboux_vect_vec,
																					sqrt_stiffness_mat);

		// Add the cost function to the problem
		problem.AddResidualBlock(stretch_bend_twist_cost, 
								loss_function, 
								dlo_state.vertices[id0].data(),
								dlo_state.quaternions[id0].coeffs().data(),
								dlo_state.vertices[id1].data(),
								dlo_state.quaternions[id1].coeffs().data());

		// Add the quaternion manifold to the problem
		problem.SetManifold(dlo_state.quaternions[id0].coeffs().data(), quaternion_manifold);
		problem.SetManifold(dlo_state.quaternions[id1].coeffs().data(), quaternion_manifold);

	}

	// -----------------------------------------

	// -----------------------------------------
	// Create gravity energy

	Real dlo_length = 0.0;

	// Set masses and calculate the total length of the DLO
	std::vector<Real> masses(num_particles);

	for (int i = 0; i < num_particles; i++)
	{
		const Real &l = dlo_state.segment_lengths[i]; // segment length
        Real V = M_PI*(energy_model_params.radius*energy_model_params.radius)*l;  // segment volume
        Real mass = V * energy_model_params.density; // segment mass

		masses[i] = mass;
		dlo_length += l;
	}
	
	// // Calculate the reference height for gravity energy
	// Real gravity_ref_height = gravityRefHeight(dlo_state.vertices, dlo_length);
	// ROS_INFO("gravity_ref_height: %f", gravity_ref_height);

	// Calculate the reference height for gravity energy (vector version)
	const Eigen::Matrix<Real,3,1> gravity_ref_point = gravityRefPoint(dlo_state.vertices, dlo_length);
	
	// // Define the gravity ref point to a minimum number of the system
	// // std::numeric_limits<Real>::min();
	// const Eigen::Matrix<Real,3,1> gravity_ref_point(-1000000.0,
	// 												-1000000.0,
	// 												-1000000.0);
	
	ROS_INFO("Gravity ref point: %f, %f, %f", gravity_ref_point(0), gravity_ref_point(1), gravity_ref_point(2));


	// gravity energy cost
	for (int k = 0; k < num_particles; k++)
	{
		ceres::CostFunction* gravity_energy_cost = GravityEnergyCost::Create(masses[k], 
																			energy_model_params.gravity*1.43e+2,
																			gravity_ref_point);

		problem.AddResidualBlock(gravity_energy_cost,
								loss_function,
								dlo_state.vertices[k].data());
	}
	

	/* TODO
		// Set the zero-stretch bending twisting constraints
	for (int i = 1; i < (num_quaternions - 1); i++)
	{
		// Segment id pair of each constraint
		const int id0 = i - 1;
		const int id1 = i;
		const int id2 = i + 1;
		
		// // compute average segment length
		const Real average_segment_length0 = 0.5 * (dlo_state.segment_lengths[id0] + dlo_state.segment_lengths[id1]);
		const Real average_segment_length1 = 0.5 * (dlo_state.segment_lengths[id1] + dlo_state.segment_lengths[id2]);

		// // compute rest Darboux vector 
		// // if initially the dlo is not straight, use this:
		// // Eigen::Quaternion<Real> rest_darboux_vect = dlo_state.quaternions[id0].conjugate() * dlo_state.quaternions[id1]; 
		// // if initially the dlo is straight, use this
		// Eigen::Quaternion<Real> rest_darboux_vect = Eigen::Quaternion<Real>::Identity(); 
		// const Eigen::Matrix<Real,3,1> rest_darboux_vect_vec = (2.0 / average_segment_length) * rest_darboux_vect.vec();

		// locally each constraint position is along the z axis at the half length of the segment
		const Eigen::Matrix<Real, 3, 1> connector00(0, 0, 0.5 * dlo_state.segment_lengths[id0]);
		const Eigen::Matrix<Real, 3, 1> connector01(0, 0, -0.5 * dlo_state.segment_lengths[id1]);
		const Eigen::Matrix<Real, 3, 1> connector10(0, 0, 0.5 * dlo_state.segment_lengths[id1]);
		const Eigen::Matrix<Real, 3, 1> connector11(0, 0, -0.5 * dlo_state.segment_lengths[id2]);

		// Set sqrt stiffness matrix
		Real stiffness = 1.0e+4;
		Eigen::Matrix<double, 3, 3> sqrt_stiffness_mat = Eigen::Matrix<double, 3, 3>::Identity();
		
		// Set the diagonal elements of the stiffness matrix
		sqrt_stiffness_mat(0, 0) = stiffness;
		sqrt_stiffness_mat(1, 1) = stiffness;
		sqrt_stiffness_mat(2, 2) = stiffness;


		// // Create the cost function
		// ceres::CostFunction* gravity_force_cost = GravityForceCost::Create(average_segment_length,
		// 																			connector0,
		// 																			connector1,
		// 																			rest_darboux_vect_vec,
		// 																			sqrt_stiffness_mat);

		// // Add the cost function to the problem
		// problem.AddResidualBlock(gravity_force_cost, 
		// 						loss_function, 
		// 						dlo_state.vertices[id0].data(),
		// 						dlo_state.quaternions[id0].coeffs().data(),
		// 						dlo_state.vertices[id1].data(),
		// 						dlo_state.quaternions[id1].coeffs().data());

		// // Add the quaternion manifold to the problem
		// problem.SetManifold(dlo_state.quaternions[id0].coeffs().data(), quaternion_manifold);
		// problem.SetManifold(dlo_state.quaternions[id1].coeffs().data(), quaternion_manifold);

	}
	*/



	// -----------------------------------------


	// Constraint the problem for static particles
	for (int i = 0; i < static_particles.size(); i++)
	{
		int id = static_particles[i];
		problem.SetParameterBlockConstant(dlo_state.vertices[id].data());
		problem.SetParameterBlockConstant(dlo_state.quaternions[id].coeffs().data());
	}


	// -----------------------------------------
	// Run the solver!
	ceres::Solver::Options ceres_options;

	// Set up ceres_options using the function
	SetCeresOptions(optimize_options, &ceres_options);

	ceres::Solver::Summary summary;

	// solving
	std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
	
	ceres::Solve(ceres_options, &problem, &summary);
	
	std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
	std::chrono::duration<double> time_used_solving = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

	if (optimize_options.print_log)
	{
		std::cout << summary.FullReport() << ", Time: " << time_used_solving.count() << "s" << std::endl;
	}
	// -----------------------------------------

	
	// --------------------------------------------------------------------------
	auto end_optimize = std::chrono::high_resolution_clock::now();
    auto duration_optimize = std::chrono::duration_cast<std::chrono::microseconds>(end_optimize - start_optimize).count();
	ROS_INFO("Time taken by optimizeShape: %f ms", duration_optimize/1000.0);

	return dlo_state;
}

// Function to set up ceres::Solver::Options based on OptimizeOptions
void EnergyBasedSolver::SetCeresOptions(const OptimizeOptions& optimize_options, ceres::Solver::Options* ceres_options)
{
    // Set the minimizer type
    if (optimize_options.minimizer_type == "line_search")
        ceres_options->minimizer_type = ceres::LINE_SEARCH;
    else if (optimize_options.minimizer_type == "trust_region")
        ceres_options->minimizer_type = ceres::TRUST_REGION;
    else
        std::cerr << "Ceres doesn't support minimizer_type: " << optimize_options.minimizer_type << std::endl;

    // Set the line search direction type
    if (optimize_options.line_search_direction_type == "BFGS")
        ceres_options->line_search_direction_type = ceres::BFGS;
    else if (optimize_options.line_search_direction_type == "LBFGS")
        ceres_options->line_search_direction_type = ceres::LBFGS;
    else if (optimize_options.line_search_direction_type == "NCG")
        ceres_options->line_search_direction_type = ceres::NONLINEAR_CONJUGATE_GRADIENT;
    else if (optimize_options.line_search_direction_type == "STEEPEST_DESCENT")
        ceres_options->line_search_direction_type = ceres::STEEPEST_DESCENT;
    else
        std::cerr << "Ceres doesn't support line_search_direction_type: " << optimize_options.line_search_direction_type << std::endl;

    // Set the trust region strategy type
    if (optimize_options.trust_region_strategy_type == "DOGLEG")
        ceres_options->trust_region_strategy_type = ceres::DOGLEG;
    else if (optimize_options.trust_region_strategy_type == "LEVENBERG_MARQUARDT")
        ceres_options->trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    else
        std::cerr << "Ceres doesn't support trust_region_strategy_type: " << optimize_options.trust_region_strategy_type << std::endl;

    // Set the linear solver type
    if (optimize_options.linear_solver_type == "SPARSE_NORMAL_CHOLESKY")
        ceres_options->linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    else if (optimize_options.linear_solver_type == "DENSE_QR")
        ceres_options->linear_solver_type = ceres::DENSE_QR;
    else
        std::cerr << "Ceres doesn't support linear_solver_type: " << optimize_options.linear_solver_type << std::endl;

    // Set other options
    ceres_options->max_num_iterations = static_cast<int>(optimize_options.max_num_iterations);
    ceres_options->max_num_line_search_step_size_iterations = static_cast<int>(optimize_options.max_num_line_search_step_size_iterations);

    ceres_options->function_tolerance = optimize_options.function_tolerance;
    ceres_options->gradient_tolerance = optimize_options.gradient_tolerance;
    ceres_options->parameter_tolerance = optimize_options.parameter_tolerance;

    ceres_options->num_threads = optimize_options.num_threads;
    ceres_options->initial_trust_region_radius = optimize_options.initial_trust_region_radius;
    ceres_options->max_trust_region_radius = optimize_options.max_trust_region_radius;

    ceres_options->minimizer_progress_to_stdout = optimize_options.minimizer_progress_to_stdout;
}


void EnergyBasedSolver::saveDLOStateAsNpy(pbd_object::MeshDLO &dlo_state, 
										std::string file_name, 
										bool append_current_time=false)
{
	// Save the DLO state as a .npy file
	// Similar to below: 
	// TODO: Implement this function
    // std::vector<Real> data;
    // int num_fps = dlo_state.fps_pos_.size() / 3;
    // size_t one_state_dim = 3 * num_fps + 4 + 4 + 1;

    // Utils::pushBackEigenVecToStdVec(data, dlo_state.fps_pos_);
    // Utils::pushBackEigenVecToStdVec(data, dlo_state.end_quat_0_);
    // Utils::pushBackEigenVecToStdVec(data, dlo_state.end_quat_1_);
    // data.push_back(dlo_state.angle_b2m_(dlo_state.angle_b2m_.size() - 1));

    // Utils::createDirectory(file_name);
    // cnpy::npy_save(file_name, &data[0], {one_state_dim}, "w");
}


/**
 * @brief Compute the reference point for gravity energy in 3D space.
 */
Eigen::Matrix<Real, 3, 1> EnergyBasedSolver::gravityRefPoint(
    const std::vector<Eigen::Matrix<Real, 3, 1>> &vertices,
    const Real &dlo_length)
{
    // Initialize min_x, min_y, min_z with very large values
    Real min_x = std::numeric_limits<Real>::max();
    Real min_y = std::numeric_limits<Real>::max();
    Real min_z = std::numeric_limits<Real>::max();

    // Iterate through each vertex to find the minimum x, y, z coordinates
    for (const auto& vertex : vertices)
    {
        Real x = vertex(0);
        Real y = vertex(1);
        Real z = vertex(2);

        if (x < min_x) { min_x = x; }
        if (y < min_y) { min_y = y; }
        if (z < min_z) { min_z = z; }
    }

    // Compute the reference point by subtracting dlo_length from each min coordinate
    Eigen::Matrix<Real, 3, 1> gravity_ref_point;
    gravity_ref_point << (min_x - dlo_length), (min_y - dlo_length), (min_z - dlo_length);

    return gravity_ref_point;
}
