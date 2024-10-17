/*
 * Author: Burak Aksoy
 */

#include "dlo_simulator_stiff_rods/dlo_simulator.h"

// using namespace std;
using namespace dlo_simulator;

DloSimulator::DloSimulator(ros::NodeHandle &nh, ros::NodeHandle &nh_local, boost::recursive_mutex &mtx): 
    nh_(nh), 
    nh_local_(nh_local),
    mtx_(mtx)
{
    p_active_ = false;

    time_frames_ = 0;
    time_sum_ = 0.0;

    marker_id_ = 3; // for bookkeeping of the visualization markers used in the code.

    is_auto_sim_rate_set_ = false; 

    is_rob_01_attached_ = false;
    is_rob_02_attached_ = false;
    is_rob_03_attached_ = false;
    is_rob_04_attached_ = false;

    rob_01_attached_id_ = -1; //note: ids start from 0. -1 would be a null id.
    rob_02_attached_id_ = -1;
    rob_03_attached_id_ = -1;
    rob_04_attached_id_ = -1;
    /*
    rob_01_attached_force_.setZero();
    rob_02_attached_force_.setZero();
    rob_03_attached_force_.setZero();
    rob_04_attached_force_.setZero();
    */

    // Initialize Timers with deafault period (note: last 2 false mean: oneshot=false, autostart=false)
    timer_render_ = nh_.createTimer(ros::Duration(1.0), &DloSimulator::render, this,false, false); 
    timer_simulate_ = nh_.createTimer(ros::Duration(1.0), &DloSimulator::simulate, this,false, false); 
    /*
    timer_wrench_pub_ = nh_.createTimer(ros::Duration(1.0), &DloSimulator::publishWrenches, this,false, false); 
    */
    timer_render_rb_ = nh_.createTimer(ros::Duration(1.0), &DloSimulator::renderRigidBodies, this,false, false); 
    timer_min_dist_to_rb_pub_ = nh_.createTimer(ros::Duration(1.0), &DloSimulator::publishMinDistancesToRigidBodies, this,false, false); 
    timer_dlo_state_pub_ = nh_.createTimer(ros::Duration(1.0), &DloSimulator::publishDloState, this,false, false); 

    // Initilize parameters
    params_srv_ = nh_local_.advertiseService("params", &DloSimulator::updateParams, this);
    initialize(); // Calls updateParams() 
}

DloSimulator::~DloSimulator() {
    nh_local_.deleteParam("gravity_x");
    nh_local_.deleteParam("gravity_y");
    nh_local_.deleteParam("gravity_z");

    nh_local_.deleteParam("dt");
    nh_local_.deleteParam("set_sim_rate_auto");
    nh_local_.deleteParam("num_substeps");
    nh_local_.deleteParam("num_steps");

    nh_local_.deleteParam("is_collision_handling_enabled");
    nh_local_.deleteParam("visualize_min_distances");

    nh_local_.deleteParam("dlo_visualization_mode");
    nh_local_.deleteParam("rb_visualization_mode");

    nh_local_.deleteParam("dlo_l");
    nh_local_.deleteParam("dlo_r");
    nh_local_.deleteParam("dlo_density");
    nh_local_.deleteParam("dlo_num_segments");

    nh_local_.deleteParam("dlo_zero_stretch_stiffness");
    nh_local_.deleteParam("dlo_young_modulus");
    nh_local_.deleteParam("dlo_torsion_modulus");

    nh_local_.deleteParam("use_zero_stretch_stiffness");

    nh_local_.deleteParam("contact_tolerance");
    nh_local_.deleteParam("contact_sdf_domain_offset");

    nh_local_.deleteParam("global_damp_coeff_v");
    nh_local_.deleteParam("global_damp_coeff_w");

    nh_local_.deleteParam("initial_height");
    nh_local_.deleteParam("dlo_translation");
    nh_local_.deleteParam("dlo_rotationAxis");
    nh_local_.deleteParam("dlo_rotationAngle");
    nh_local_.deleteParam("dlo_scale");

    nh_local_.deleteParam("num_hang_corners");
    nh_local_.deleteParam("custom_static_particles");

    nh_local_.deleteParam("custom_static_particles_odom_topic_prefix");
    nh_local_.deleteParam("custom_static_particles_cmd_vel_topic_prefix");

    nh_local_.deleteParam("simulation_rate");
    nh_local_.deleteParam("rendering_rate");
    /*
    nh_local_.deleteParam("wrench_pub_rate");
    */
    nh_local_.deleteParam("rendering_rb_rate");
    nh_local_.deleteParam("min_dist_to_rb_pub_rate");
    nh_local_.deleteParam("dlo_state_pub_rate");

    nh_local_.deleteParam("rb_scene_config_path");

    nh_local_.deleteParam("dlo_state_topic_name");
    nh_local_.deleteParam("dlo_markers_topic_name");
    nh_local_.deleteParam("dlo_frame_id");

    nh_local_.deleteParam("point_marker_scale");

    nh_local_.deleteParam("point_marker_color_rgba");

    nh_local_.deleteParam("line_marker_scale_multiplier");
    
    nh_local_.deleteParam("line_marker_color_rgba");

    nh_local_.deleteParam("frame_marker_scale");
    nh_local_.deleteParam("frame_marker_axis_length");

    nh_local_.deleteParam("rb_markers_topic_name");

    nh_local_.deleteParam("min_dist_to_rb_topic_name");
    nh_local_.deleteParam("min_dist_markers_topic_name");

    nh_local_.deleteParam("change_particle_dynamicity_topic_name");
    nh_local_.deleteParam("set_particle_dynamicity_service_name");
 
    nh_local_.deleteParam("rb_line_marker_scale_multiplier");

    nh_local_.deleteParam("rb_line_marker_color_rgba");

    nh_local_.deleteParam("rb_mesh_marker_color_rgba");

    nh_local_.deleteParam("min_dist_line_marker_scale_multiplier");

    nh_local_.deleteParam("min_dist_line_marker_color_rgba");

    /*
    nh_local_.deleteParam("wrench_01_topic_name");
    nh_local_.deleteParam("wrench_02_topic_name");
    nh_local_.deleteParam("wrench_03_topic_name");
    nh_local_.deleteParam("wrench_04_topic_name");

    nh_local_.deleteParam("wrench_01_frame_id");
    nh_local_.deleteParam("wrench_02_frame_id");
    nh_local_.deleteParam("wrench_03_frame_id");
    nh_local_.deleteParam("wrench_04_frame_id");
    */

}

bool DloSimulator::updateParams(std_srvs::Empty::Request& req, std_srvs::Empty::Response& res)
{
    bool prev_active = p_active_;

    // Get parameters from the parameter server
    nh_local_.param<bool>("active", p_active_, true);
    nh_local_.param<bool>("reset", p_reset_, false);

    nh_local_.param<Real>("gravity_x", gravity_x_, 0.0);
    nh_local_.param<Real>("gravity_y", gravity_y_, 0.0);
    nh_local_.param<Real>("gravity_z", gravity_z_, -9.81);
    
    nh_local_.param<Real>("dt", dt_, 1.0 / 100.0); //200
    nh_local_.param<bool>("set_sim_rate_auto", set_sim_rate_auto_, false); // to set the simulation rate and dt automatically

    nh_local_.param<int>("num_substeps", num_substeps_, 30); 
    nh_local_.param<int>("num_steps", num_steps_, 1);

    nh_local_.param<bool>("is_collision_handling_enabled", is_collision_handling_enabled_, true);
    nh_local_.param<bool>("visualize_min_distances", visualize_min_distances_, true);

    nh_local_.param<int>("dlo_visualization_mode", dlo_visualization_mode_, 4); // 0: Points Only, 1: Line Segments Only, 2: Orientation Frames Only (TODO), 3: Points and Line Segments, 4: All
    nh_local_.param<int>("rb_visualization_mode", rb_visualization_mode_, 2); // 0: Mesh Only, 1: Wireframe Only, 2: Both
    
    nh_local_.param<Real>("dlo_l", dlo_l_, 3.); 
    nh_local_.param<Real>("dlo_r", dlo_r_, 0.005); 
    nh_local_.param<Real>("dlo_density", dlo_density_, 1000.0);
    nh_local_.param<int>("dlo_num_segments", dlo_num_segments_, 10);

    nh_local_.param<Real>("dlo_zero_stretch_stiffness", dlo_zero_stretch_stiffness_, 0.01);
    nh_local_.param<Real>("dlo_young_modulus", dlo_young_modulus_, 0.01);
    nh_local_.param<Real>("dlo_torsion_modulus", dlo_torsion_modulus_, 0.01);

    nh_local_.param<bool>("use_zero_stretch_stiffness", use_zero_stretch_stiffness_, true);

    nh_local_.param<Real>("contact_tolerance", contact_tolerance_, static_cast<Real>(0.1));
    nh_local_.param<Real>("contact_sdf_domain_offset", contact_sdf_domain_offset_, static_cast<Real>(2.0));

    nh_local_.param<Real>("global_damp_coeff_v", global_damp_coeff_v_, 0.0);
    nh_local_.param<Real>("global_damp_coeff_w", global_damp_coeff_w_, 0.0);

    nh_local_.param<Real>("initial_height", initial_height_, 0.0);
    nh_local_.param("dlo_translation", dlo_translation_, std::vector<Real>({0.0, 0.0, 0.0}));
    nh_local_.param("dlo_rotationAxis", dlo_rotationAxis_, std::vector<Real>({0.0, 0.0, 1.0}));
    nh_local_.param<Real>("dlo_rotationAngle", dlo_rotationAngle_, 0.0);
    nh_local_.param("dlo_scale", dlo_scale_, std::vector<Real>({1.0, 1.0, 1.0}));

    nh_local_.param<int>("num_hang_corners", num_hang_corners_, 0);
    nh_local_.param("custom_static_particles", custom_static_particles_, std::vector<int>());

    nh_local_.param<std::string>("custom_static_particles_odom_topic_prefix", custom_static_particles_odom_topic_prefix_, std::string("custom_static_particles_odom_"));
    nh_local_.param<std::string>("custom_static_particles_cmd_vel_topic_prefix", custom_static_particles_cmd_vel_topic_prefix_, std::string("custom_static_particles_cmd_vel_"));

    nh_local_.param<Real>("simulation_rate", simulation_rate_, 90.0); 
    nh_local_.param<Real>("rendering_rate", rendering_rate_, 30.0);
    /*
    nh_local_.param<Real>("wrench_pub_rate", wrench_pub_rate_, 60.0); //60
    */
    nh_local_.param<Real>("rendering_rb_rate", rendering_rb_rate_, 1.0); //1
    nh_local_.param<Real>("min_dist_to_rb_pub_rate", min_dist_to_rb_pub_rate_, 60.0); //60.0
    nh_local_.param<Real>("dlo_state_pub_rate", dlo_state_pub_rate_, 60.0); //60.0

    nh_local_.param<std::string>("rb_scene_config_path", rb_scene_config_path_, std::string(""));

    nh_local_.param<std::string>("dlo_state_topic_name", dlo_state_topic_name_, std::string("dlo_state"));
    nh_local_.param<std::string>("dlo_markers_topic_name", dlo_markers_topic_name_, std::string("dlo_markers"));
    nh_local_.param<std::string>("dlo_frame_id", dlo_frame_id_, std::string("map"));
    
    /*
    nh_local_.param<std::string>("wrench_01_topic_name", wrench_01_topic_name_, std::string("d1/dlo_wrench_stamped"));
    nh_local_.param<std::string>("wrench_02_topic_name", wrench_02_topic_name_, std::string("d2/dlo_wrench_stamped"));
    nh_local_.param<std::string>("wrench_03_topic_name", wrench_03_topic_name_, std::string("d3/dlo_wrench_stamped"));
    nh_local_.param<std::string>("wrench_04_topic_name", wrench_04_topic_name_, std::string("d4/dlo_wrench_stamped"));

    nh_local_.param<std::string>("wrench_01_frame_id", wrench_01_frame_id_, std::string("d1_tf_dlo_mount_link"));
    nh_local_.param<std::string>("wrench_02_frame_id", wrench_02_frame_id_, std::string("d2_tf_dlo_mount_link"));
    nh_local_.param<std::string>("wrench_03_frame_id", wrench_03_frame_id_, std::string("d3_tf_dlo_mount_link"));
    nh_local_.param<std::string>("wrench_04_frame_id", wrench_04_frame_id_, std::string("d4_tf_dlo_mount_link"));
    */

    nh_local_.param<Real>("point_marker_scale",   point_marker_scale_,   0.015);

    nh_local_.param("point_marker_color_rgba", point_marker_color_rgba_, std::vector<Real>({1.0, 0.0, 1.0, 1.0}));

    nh_local_.param<Real>("line_marker_scale_multiplier", line_marker_scale_multiplier_, 1.0);

    nh_local_.param("line_marker_color_rgba", line_marker_color_rgba_, std::vector<Real>({0.0, 1.0, 1.0, 1.0}));

    nh_local_.param<Real>("frame_marker_scale",   frame_marker_scale_, 0.01);
    nh_local_.param<Real>("frame_marker_axis_length",   frame_marker_axis_length_, 0.05);

    nh_local_.param<std::string>("rb_markers_topic_name", rb_markers_topic_name_, std::string("rigid_body_markers"));
    nh_local_.param<std::string>("min_dist_to_rb_topic_name", min_dist_to_rb_topic_name_, std::string("min_dist_to_rigid_bodies"));
    nh_local_.param<std::string>("min_dist_markers_topic_name", min_dist_markers_topic_name_, std::string("min_dist_markers"));

    nh_local_.param<std::string>("change_particle_dynamicity_topic_name", change_particle_dynamicity_topic_name_, std::string("change_particle_dynamicity"));
    nh_local_.param<std::string>("set_particle_dynamicity_service_name", set_particle_dynamicity_service_name_, std::string("set_particle_dynamicity"));

    nh_local_.param<std::string>("set_dlo_young_modulus_service_name", set_dlo_young_modulus_service_name_, std::string("set_dlo_young_modulus"));
    nh_local_.param<std::string>("set_dlo_torsion_modulus_service_name", set_dlo_torsion_modulus_service_name_, std::string("set_dlo_torsion_modulus"));
    nh_local_.param<std::string>("get_dlo_young_modulus_service_name", get_dlo_young_modulus_service_name_, std::string("get_dlo_young_modulus"));
    nh_local_.param<std::string>("get_dlo_torsion_modulus_service_name", get_dlo_torsion_modulus_service_name_, std::string("get_dlo_torsion_modulus"));
    nh_local_.param<std::string>("change_dlo_young_modulus_topic_name", change_dlo_young_modulus_topic_name_, std::string("change_dlo_young_modulus"));
    nh_local_.param<std::string>("change_dlo_torsion_modulus_topic_name", change_dlo_torsion_modulus_topic_name_, std::string("change_dlo_torsion_modulus"));
    
    nh_local_.param<std::string>("enable_collision_handling_service_name", enable_collision_handling_service_name_, std::string("enable_collision_handling"));
    
    nh_local_.param<Real>("rb_line_marker_scale_multiplier", rb_line_marker_scale_multiplier_, 1.0);
    
    nh_local_.param("rb_line_marker_color_rgba", rb_line_marker_color_rgba_, std::vector<Real>({0.0, 0.0, 1.0, 1.0}));
    
    nh_local_.param("rb_mesh_marker_color_rgba", rb_mesh_marker_color_rgba_, std::vector<Real>({0.5, 0.5, 0.5, 0.7}));

    nh_local_.param<Real>("min_dist_line_marker_scale_multiplier", min_dist_line_marker_scale_multiplier_, 1.0);

    nh_local_.param("min_dist_line_marker_color_rgba", min_dist_line_marker_color_rgba_, std::vector<Real>({0.0, 0.0, 0.0, 0.5}));

    // Set timer periods based on the parameters
    timer_render_.setPeriod(ros::Duration(1.0/rendering_rate_));
    timer_simulate_.setPeriod(ros::Duration(1.0/simulation_rate_));
    /*
    timer_wrench_pub_.setPeriod(ros::Duration(1.0/wrench_pub_rate_));
    */
    timer_render_rb_.setPeriod(ros::Duration(1.0/rendering_rb_rate_));
    timer_min_dist_to_rb_pub_.setPeriod(ros::Duration(1.0/min_dist_to_rb_pub_rate_));
    timer_dlo_state_pub_.setPeriod(ros::Duration(1.0/dlo_state_pub_rate_));

    // Initilize gravity vector
    gravity_ << gravity_x_, gravity_y_, gravity_z_;    

    //Create mesh
    std::string dlo_name = "dlo";
    pbd_object::MeshDLO dlo_mesh = DloSimulator::createMeshDLO(dlo_name, dlo_l_, initial_height_, dlo_num_segments_);

    // std::cout << "dlo_mesh.name: " << dlo_mesh.name << std::endl;

    // Transform the created mesh
    dlo_mesh = DloSimulator::transformMeshDLO(dlo_mesh,
                                            dlo_translation_,
                                            dlo_rotationAxis_,
                                            dlo_rotationAngle_,
                                            dlo_scale_); // scale is not used

    // Create dlo
    dlo_ = pbd_object::Dlo(dlo_mesh, 
                           dlo_zero_stretch_stiffness_, 
                           dlo_young_modulus_, 
                           dlo_torsion_modulus_, 
                           dlo_density_,dlo_r_,
                           use_zero_stretch_stiffness_,
                           global_damp_coeff_v_,
                           global_damp_coeff_w_);

    // Set static particles
    dlo_.setStaticParticles(custom_static_particles_);

    // Rigid body scene setting 
    // (sets rigid_bodies_ vector which holds the data related to all the rigid bodies in the scene)
    if(!rb_scene_config_path_.empty()){
        // read the rigid body scene json file with rigid_body_scene_loader utility class
        utilities::RigidBodySceneLoader *rb_scene_loader = new utilities::RigidBodySceneLoader();
        rb_scene_loader->readScene(rb_scene_config_path_, rigid_bodies_);

        // Print all the read data from rigid_bodies_ vector for debugging
        /*
        std::cout << "---------------------------------------" << std::endl;
        for (const utilities::RigidBodySceneLoader::RigidBodyData& rbd : rigid_bodies_) {
            std::cout << "ID: " << rbd.m_id << std::endl;
            std::cout << "Model File: " << rbd.m_modelFile << std::endl;
            std::cout << "Is Dynamic: " << std::boolalpha << rbd.m_isDynamic << std::endl;
            std::cout << "Density: " << rbd.m_density << std::endl;
            std::cout << "Translation: " << rbd.m_x << std::endl;
            std::cout << "Rotation Quaternion: " << rbd.m_q.coeffs() << std::endl; 
            // Note: Eigen quaternions store coefficients as (x, y, z, w)
            std::cout << "Scale: " << rbd.m_scale << std::endl;
            std::cout << "Linear Velocity: " << rbd.m_v << std::endl;
            std::cout << "Angular Velocity: " << rbd.m_omega << std::endl;
            std::cout << "Restitution Coefficient: " << rbd.m_restitutionCoeff << std::endl;
            std::cout << "friction Static Coefficient: " << rbd.m_frictionCoeffStatic << std::endl;
            std::cout << "friction Dynamic Coefficient: " << rbd.m_frictionCoeffDynamic << std::endl;
            std::cout << "Collision Object File Name: " << rbd.m_collisionObjectFileName << std::endl;
            std::cout << "Collision Object Scale: " << rbd.m_collisionObjectScale << std::endl;
            std::cout << "Resolution SDF: " << rbd.m_resolutionSDF << std::endl;
            std::cout << "Invert SDF: " << std::boolalpha << rbd.m_invertSDF << std::endl;
            std::cout << "---------------------------------------" << std::endl;
        }
        */

        // For each rigid body entry in the json file
        // load its mesh (read from .obj file using initial height = 0
        // transform its mesh (translate, rotate, scale)  (overload the transform Mesh function with accepting a quaternion)
        // update the the created rigid body data's pbd_object::Mesh object with the loaded mesh
        for (utilities::RigidBodySceneLoader::RigidBodyData& rbd : rigid_bodies_) {
            // Load the mesh
            pbd_object::Mesh mesh = loadMesh("RigidBodyMesh_" + std::to_string(rbd.m_id), rbd.m_modelFile, 0);

            // Transform the mesh
            mesh = transformMesh(mesh, rbd.m_x, rbd.m_q.normalized(), rbd.m_scale);

            // Update the RigidBodyData's mesh
            rbd.m_mesh = mesh;
        }

        delete rb_scene_loader; // Don't forget to delete the loader to free memory
    }

    // For each rigid body entry in the json file
    // Set the Signed Distance Field function for the rigid bodies
    // TODO: Make this for loop a function
    for (utilities::RigidBodySceneLoader::RigidBodyData& rbd : rigid_bodies_) {
        // use Discregrid Library to set discrete signed distance fields of the rigid bodies

        // This is for the parallel discretization of (preferably smooth) functions on regular grids. 
	    // This is especially well-suited for the discretization of signed distance functions. 

        // get the sdf file name entry in the json file
        std::string sdf_file_name = rbd.m_collisionObjectFileName; 

        // check if collision sdf file is speficied in the json file
        std::ifstream sdf_file(sdf_file_name); 

        // if the sdf is specifed and it actually exists
		if ((sdf_file_name != "") && (sdf_file.good()))
		{
            // load the sdf file
            // append the loaded sdf file to the rigidBodyData vector(rigid_bodies_)

            std::cout << "Load SDF file: " << sdf_file_name << std::endl;
            rbd.m_discregrid_ptr = std::make_shared<utilities::CollisionHandler::Grid>(sdf_file_name);

		}
        else // sdf file  is not specified in the json file or not readable
        {
            // Therefore,
            // With the Discregrid library generate a (cubic) polynomial discretization given: 
            // a box-shaped domain (A), 
            // a grid resolution (B), 
            // and a function that maps a 3D position in space to a real scalar value.(C)
            // It can also serialize and save the discretization to a file (D)


            // First we Create Discregrid's TriangleMesh object from 
            // the object's vertexData, mesh faces data, and number of faces.

            // std::vector<double> doubleVec;
            // doubleVec.resize(3 * rbd.m_mesh.vertices.rows());
            // for (unsigned int i = 0; i < rbd.m_mesh.vertices.rows(); i++)
            //     for (unsigned int j = 0; j < 3; j++)
            //         doubleVec[3 * i + j] = rbd.m_mesh.vertices.row(i)[j];

            // std::vector<unsigned int> facesVec;
            // facesVec.resize(3 * rbd.m_mesh.face_tri_ids.rows());
            // for (unsigned int i = 0; i < rbd.m_mesh.face_tri_ids.rows(); i++)
            //     for (unsigned int j = 0; j < 3; j++)
            //         facesVec[3 * i + j] = rbd.m_mesh.face_tri_ids.row(i)[j];
            
            // Discregrid::TriangleMesh sdfMesh(&doubleVec[0], facesVec.data(), rbd.m_mesh.vertices.rows(), rbd.m_mesh.face_tri_ids.rows());
            // Note that if you do this the SDF object pose and scale are defaulted to the transformed object, 
            // This is because the mesh is already transformed. So you may need to update your json file accordingly.

            // OR another option is to Create Discregrid's TriangleMesh object directly 
            // reading the obj file using the discregrid library's file reader function:
            Discregrid::TriangleMesh sdfMesh(rbd.m_modelFile);

            // --------------- (A) -------------------
            // Now create the box shaped domain of the discretization from the size of 
            // the AABB(Axis Aligned Bounding Box) of the Discregrid's TriangleMesh object.
            Eigen::AlignedBox3d domain;
            for (auto const& x : sdfMesh.vertices())
            {
                domain.extend(x);
            }

            // This offset defines the minimum distance to a rigid body that will start reporting minimum distance readings
            // Eigen::Vector3d offsetVector = contact_sdf_domain_offset_ * Eigen::Vector3d::Ones();
            Eigen::Vector3d offsetVector = (contact_sdf_domain_offset_ * Eigen::Vector3d::Ones()).array() / rbd.m_scale.transpose().array();
            std::cout << "SDF domain offset :" << offsetVector << std::endl;
            domain.max() += offsetVector;
            domain.min() -= offsetVector;

            // --------------- (B) -------------------
            // Specify the grid resolution from the scene (.json) file
            std::cout << "Set SDF resolution: " 
                        << rbd.m_resolutionSDF[0] << ", " 
                        << rbd.m_resolutionSDF[1] << ", " 
                        << rbd.m_resolutionSDF[2] << std::endl;
            
            std::array<unsigned int, 3> resolution({rbd.m_resolutionSDF[0], 
                                                    rbd.m_resolutionSDF[1], 
                                                    rbd.m_resolutionSDF[2] });

            // Now generate the SDF Grid object by creating the 
            // Discregrid::CubicLagrangeDiscreteGrid object,
            // with the domain and the specified resolution:
            rbd.m_discregrid_ptr = std::make_shared<utilities::CollisionHandler::Grid>(domain, resolution);

            // --------------- (C) -------------------
            // Now specify the function that maps a 3D position in space to a real scalar value.
            
            // We need the signed distances hence we
            // create object of a TriangleMeshDistance class data structure that directly provides 
            // the capability to compute and discretize signed distance fields to triangle meshes.
            Discregrid::TriangleMeshDistance md(sdfMesh); 

            // create a function object "func" using the default constructor of Discregrid::DiscreteGrid::ContinuousFunction
            auto func = Discregrid::DiscreteGrid::ContinuousFunction{}; 
            // Reassign func to a lambda function that takes an Eigen::Vector3d input and returns the distance from the input point to the mesh.
            func = [&md](Eigen::Vector3d const& xi) {return md.signed_distance(xi).distance; };
            
            // add the function "func" to the distance field grid corresponding to sdfFileName.
            // This generated the discretization on the initiated grid.
            rbd.m_discregrid_ptr->addFunction(func, true); 

            // --------------- (D) -------------------
            // Create the sdf file for later use.

            // We can serialize this discretization to an output file:
            
            // Log the message that an SDF is being generated for the model file.
            std::cout << "Generate SDF for model file: " << rbd.m_modelFile << std::endl;

            // create a file path and name for the collision object sdf file if it is not specified in the json file
            const std::string resStr = std::to_string(rbd.m_resolutionSDF[0]) + "_" + 
                                        std::to_string(rbd.m_resolutionSDF[1]) + "_" + 
                                        std::to_string(rbd.m_resolutionSDF[2]); // sdf resolution string to append to end of obj file name
            sdf_file_name =  rbd.m_modelFile + "_" + resStr + ".csdf"; // sdf file name with full path next to the obj file.

            // Log the message that the SDF is being saved.
            std::cout << "Save SDF: " << sdf_file_name << std::endl;

            // save the sdf file to the specified path (sdf_file_name)
            rbd.m_discregrid_ptr->save(sdf_file_name);
        }
    }


    collision_handler_ = new utilities::CollisionHandler(dlo_, rigid_bodies_);
 
    collision_handler_->setContactTolerance(contact_tolerance_);

    collision_handler_->setContactCallback(contactCallbackFunction, this);

    // Add cloth object to the collision handler
    collision_handler_->addCollisionObjectWithoutGeometry(0, // unsigned int bodyIndex
                                                        utilities::CollisionHandler::CollisionObject::TriangleModelCollisionObjectType,
                                                        &(*dlo_.getPosPtr())[0], // Address of first element in the vector
                                                        dlo_.getPosPtr()->size(), // unsigned int numVertices
                                                        true);// bool testMesh)


    // Add rigid bodies to the collision handler
    // TODO: Make this for loop a function
    int i = 0;
    for (utilities::RigidBodySceneLoader::RigidBodyData& rbd : rigid_bodies_) {

        collision_handler_->addCubicSDFCollisionObject(i, // unsigned int bodyIndex
                                                    utilities::CollisionHandler::CollisionObject::RigidBodyCollisionObjectType, 
                                                    &rbd.m_mesh.vertices, // Eigen::Matrix<Real,Eigen::Dynamic,3> *
                                                    rbd.m_mesh.vertices.rows(), // int
                                                    rbd.m_discregrid_ptr, //std::shared_ptr<Discregrid::CubicLagrangeDiscreteGrid>
                                                    rbd.m_collisionObjectScale, //const Eigen::Matrix<Real, 3, 1> 
                                                    true, //bool test mesh
                                                    rbd.m_invertSDF); // bool invert sdf

        i++;
    }


    if (p_active_ != prev_active) {
        if (p_active_) {
            // Create visualization marker publisher
            pub_dlo_state_ = nh_.advertise<dlo_simulator_stiff_rods::SegmentStateArray>(dlo_state_topic_name_, 1);
            pub_dlo_marker_array_ = nh_.advertise<visualization_msgs::MarkerArray>(dlo_markers_topic_name_, 1);

            // Create a subscriber for each custom static particle
            for (const int& particle_id : custom_static_particles_) {
                // Check if a subscriber for this particle ID already exists
                if (custom_static_particles_odom_subscribers_.find(particle_id) == custom_static_particles_odom_subscribers_.end()) {
                    std::string topic = custom_static_particles_odom_topic_prefix_ + std::to_string(particle_id);
                    ros::Subscriber sub = nh_.subscribe<nav_msgs::Odometry>(topic, 1,
                                                                            [this, particle_id](const nav_msgs::Odometry::ConstPtr& odom_msg) { 
                                                                                this->odometryCb_custom_static_particles(odom_msg, particle_id); }
                                                                            );    
                    // Add the new subscriber to the map
                    custom_static_particles_odom_subscribers_[particle_id] = sub;
                }
                // Else, a subscriber for this particle ID already exists

                if (custom_static_particles_cmd_vel_subscribers_.find(particle_id) == custom_static_particles_cmd_vel_subscribers_.end()) {
                    std::string topic = custom_static_particles_cmd_vel_topic_prefix_ + std::to_string(particle_id);
                    ros::Subscriber sub = nh_.subscribe<geometry_msgs::Twist>(topic, 10,
                                                                            [this, particle_id](const geometry_msgs::Twist::ConstPtr& twist_msg) { 
                                                                                this->cmdVelCb_custom_static_particles(twist_msg, particle_id); }
                                                                            );    
                    // Add the new subscriber to the map
                    custom_static_particles_cmd_vel_subscribers_[particle_id] = sub;
                }
                // Else, a subscriber for this particle ID already exists
            }

            //Create subscribers
            sub_change_particle_dynamicity_ = nh_.subscribe(change_particle_dynamicity_topic_name_, 
                                                            1, 
                                                            &DloSimulator::changeParticleDynamicityCb, 
                                                            this);

            sub_change_young_modulus_ = nh_.subscribe(change_dlo_young_modulus_topic_name_, 
                                                        1, 
                                                        &DloSimulator::changeYoungModulusCb, 
                                                        this);

            sub_change_torsion_modulus_ = nh_.subscribe(change_dlo_torsion_modulus_topic_name_,
                                                        1,
                                                        &DloSimulator::changeTorsionModulusCb,
                                                        this);                                                        

            /*
            // Create publishers
            pub_wrench_stamped_01_ = nh_.advertise<geometry_msgs::WrenchStamped>(wrench_01_topic_name_, 1);
            pub_wrench_stamped_02_ = nh_.advertise<geometry_msgs::WrenchStamped>(wrench_02_topic_name_, 1);
            pub_wrench_stamped_03_ = nh_.advertise<geometry_msgs::WrenchStamped>(wrench_03_topic_name_, 1);
            pub_wrench_stamped_04_ = nh_.advertise<geometry_msgs::WrenchStamped>(wrench_04_topic_name_, 1);
            */

            pub_rb_marker_array_ = nh_.advertise<visualization_msgs::MarkerArray>(rb_markers_topic_name_, 1);

            pub_min_dist_to_rb_ = nh_.advertise<dlo_simulator_stiff_rods::MinDistanceDataArray>(min_dist_to_rb_topic_name_, 1);

            pub_min_dist_marker_array_ = nh_.advertise<visualization_msgs::MarkerArray>(min_dist_markers_topic_name_, 1);

            // Create Services
            set_particle_dynamicity_srv_ = nh_local_.advertiseService(set_particle_dynamicity_service_name_, 
                                                                    &DloSimulator::setParticleDynamicityCallback, this);

            set_young_modulus_srv_ = nh_local_.advertiseService(set_dlo_young_modulus_service_name_,
                                                                    &DloSimulator::setDloYoungModulusCallback, this);

            set_torsion_modulus_srv_ = nh_local_.advertiseService(set_dlo_torsion_modulus_service_name_,
                                                                    &DloSimulator::setDloTorsionModulusCallback, this);

            get_young_modulus_srv_ = nh_local_.advertiseService(get_dlo_young_modulus_service_name_,
                                                                    &DloSimulator::getDloYoungModulusCallback, this);

            get_torsion_modulus_srv_ = nh_local_.advertiseService(get_dlo_torsion_modulus_service_name_,
                                                                    &DloSimulator::getDloTorsionModulusCallback, this);

            enable_collision_handling_srv_ = nh_local_.advertiseService(enable_collision_handling_service_name_,
                                                                    &DloSimulator::enableCollisionHandlingCallback, this);

            // Start timers
            timer_simulate_.start();
            timer_render_.start();
            /*
            timer_wrench_pub_.start();
            */
            timer_render_rb_.start();
            timer_min_dist_to_rb_pub_.start();
            timer_dlo_state_pub_.start();
        }
        else {
            // Send empty message?/*

            // Stop publishers
            pub_dlo_state_.shutdown();
            pub_dlo_marker_array_.shutdown();

            // Stop subscribers
            sub_change_particle_dynamicity_.shutdown();

            sub_change_young_modulus_.shutdown();
            sub_change_torsion_modulus_.shutdown();

            // Iterate through the map and shut down each subscriber
            for (auto& kv : custom_static_particles_odom_subscribers_) {
                ros::Subscriber& subscriber = kv.second; // Get the subscriber (which is the value in the key-value pair)
                subscriber.shutdown();
            }

            for (auto& kv : custom_static_particles_cmd_vel_subscribers_) {
                ros::Subscriber& subscriber = kv.second; // Get the subscriber (which is the value in the key-value pair)
                subscriber.shutdown();
            }

            /*
            // Stop publishers
            pub_wrench_stamped_01_.shutdown();
            pub_wrench_stamped_02_.shutdown();
            pub_wrench_stamped_03_.shutdown();
            pub_wrench_stamped_04_.shutdown();
            */

            pub_rb_marker_array_.shutdown();

            pub_min_dist_to_rb_.shutdown();

            // Stop timers
            timer_render_.stop();
            timer_simulate_.stop();
            /*
            timer_wrench_pub_.stop();
            */
            timer_render_rb_.stop();
            timer_min_dist_to_rb_pub_.stop();
            timer_dlo_state_pub_.stop();
            
        }
    }

    if (p_reset_)
        reset();

    return true;
}

void DloSimulator::reset(){
    time_frames_ = 0;
    time_sum_ = 0.0;

    marker_id_ = 3;

    is_auto_sim_rate_set_ = false; 

    is_rob_01_attached_ = false;
    is_rob_02_attached_ = false;
    is_rob_03_attached_ = false;
    is_rob_04_attached_ = false;

    rob_01_attached_id_ = -1; //note: ids start from 0. -1 would be a null id.
    rob_02_attached_id_ = -1;
    rob_03_attached_id_ = -1;
    rob_04_attached_id_ = -1;

    /*
    rob_01_attached_force_.setZero();
    rob_02_attached_force_.setZero();
    rob_03_attached_force_.setZero();
    rob_04_attached_force_.setZero();
    */

    p_reset_ = false;
    nh_local_.setParam("reset",false);
}

bool DloSimulator::setParticleDynamicityCallback(dlo_simulator_stiff_rods::SetParticleDynamicity::Request &req,
                                                 dlo_simulator_stiff_rods::SetParticleDynamicity::Response &res)
{
    int particle_id = req.particle_id;
    bool is_dynamic = req.is_dynamic;
    
    try {
        dlo_.changeParticleDynamicity(particle_id, is_dynamic);
    } catch (const std::out_of_range& e) {
        ROS_ERROR("Error in setParticleDynamicity: %s", e.what());
        res.success = false; // Indicate failure
        return true; // Still return true to indicate the service call was processed
    }

    // Check if a subscriber for this particle ID exists
    auto sub_iter = custom_static_particles_odom_subscribers_.find(particle_id);
    auto sub_iter_cmd_vel = custom_static_particles_cmd_vel_subscribers_.find(particle_id);

    if (is_dynamic){  // particle is being tried to set dynamic
        if (sub_iter != custom_static_particles_odom_subscribers_.end()) {
            // Subscriber exists, shut it down
            sub_iter->second.shutdown();
            custom_static_particles_odom_subscribers_.erase(sub_iter); // Remove the subscriber from the map
            res.success = true;
        } else {
            // No subscriber was found for the given ID
            res.success = false; 
        }

        if (sub_iter_cmd_vel != custom_static_particles_cmd_vel_subscribers_.end()) {
            // Subscriber exists, shut it down
            sub_iter_cmd_vel->second.shutdown();
            custom_static_particles_cmd_vel_subscribers_.erase(sub_iter_cmd_vel); // Remove the subscriber from the map
            res.success = true;
        } else {
            // No subscriber was found for the given ID
            res.success = false; 
        }
    }
    else{ // particle is being tried to set static
        if (sub_iter != custom_static_particles_odom_subscribers_.end()) {
            // Subscriber already exists, so just return with success
            res.success = true;
        } else {
            // Subscriber does not exist, create a new one
            std::string topic = custom_static_particles_odom_topic_prefix_ + std::to_string(particle_id);
            ros::Subscriber sub = nh_.subscribe<nav_msgs::Odometry>(topic, 1,
                                                                    [this, particle_id](const nav_msgs::Odometry::ConstPtr& odom_msg) { 
                                                                        this->odometryCb_custom_static_particles(odom_msg, particle_id); }
                                                                    );
            custom_static_particles_odom_subscribers_[particle_id] = sub; // Add the new subscriber to the map

            res.success = true; // Indicate success
        }

        if (sub_iter_cmd_vel != custom_static_particles_cmd_vel_subscribers_.end()) {
            // Subscriber already exists, so just return with success
            res.success = true;
        } else {
            // Subscriber does not exist, create a new one
            std::string topic = custom_static_particles_cmd_vel_topic_prefix_ + std::to_string(particle_id);
            ros::Subscriber sub = nh_.subscribe<geometry_msgs::Twist>(topic, 10,
                                                                    [this, particle_id](const geometry_msgs::Twist::ConstPtr& twist_msg) { 
                                                                        this->cmdVelCb_custom_static_particles(twist_msg, particle_id); }
                                                                    );
            custom_static_particles_cmd_vel_subscribers_[particle_id] = sub; // Add the new subscriber to the map

            res.success = true; // Indicate success
        }
    }

    return true;
}

void DloSimulator::changeParticleDynamicityCb(const dlo_simulator_stiff_rods::ChangeParticleDynamicity::ConstPtr msg){
    int particle_id = msg->particle_id;
    bool is_dynamic = msg->is_dynamic;
    
    try {
        dlo_.changeParticleDynamicity(particle_id, is_dynamic);
    } catch (const std::out_of_range& e) {
        ROS_ERROR("Error in changeParticleDynamicityCb: %s", e.what()); // Indicate failure
        return;
    }

    // Check if a subscriber for this particle ID exists
    auto sub_iter = custom_static_particles_odom_subscribers_.find(particle_id);
    auto sub_iter_cmd_vel = custom_static_particles_cmd_vel_subscribers_.find(particle_id);

    if (is_dynamic){  // particle is being tried to set dynamic
        if (sub_iter != custom_static_particles_odom_subscribers_.end()) {
            // Subscriber exists, shut it down
            sub_iter->second.shutdown();
            custom_static_particles_odom_subscribers_.erase(sub_iter); // Remove the subscriber from the map
        } else {
            // No subscriber was found for the given ID
        }

        if (sub_iter_cmd_vel != custom_static_particles_cmd_vel_subscribers_.end()) {
            // Subscriber exists, shut it down
            sub_iter_cmd_vel->second.shutdown();
            custom_static_particles_cmd_vel_subscribers_.erase(sub_iter_cmd_vel); // Remove the subscriber from the map
        } else {
            // No subscriber was found for the given ID
        }
    }
    else{ // particle is being tried to set static
        if (sub_iter != custom_static_particles_odom_subscribers_.end()) {
            // Subscriber already exists, so just return with success
        } else {
            // Subscriber does not exist, create a new one
            std::string topic = custom_static_particles_odom_topic_prefix_ + std::to_string(particle_id);
            ros::Subscriber sub = nh_.subscribe<nav_msgs::Odometry>(topic, 1,
                                                                    [this, particle_id](const nav_msgs::Odometry::ConstPtr& odom_msg) { 
                                                                        this->odometryCb_custom_static_particles(odom_msg, particle_id); }
                                                                    );
            custom_static_particles_odom_subscribers_[particle_id] = sub; // Add the new subscriber to the map
        }

        if (sub_iter_cmd_vel != custom_static_particles_cmd_vel_subscribers_.end()) {
            // Subscriber already exists, so just return with success
        } else {
            // Subscriber does not exist, create a new one
            std::string topic = custom_static_particles_cmd_vel_topic_prefix_ + std::to_string(particle_id);
            ros::Subscriber sub = nh_.subscribe<geometry_msgs::Twist>(topic, 10,
                                                                    [this, particle_id](const geometry_msgs::Twist::ConstPtr& twist_msg) { 
                                                                        this->cmdVelCb_custom_static_particles(twist_msg, particle_id); }
                                                                    );
            custom_static_particles_cmd_vel_subscribers_[particle_id] = sub; // Add the new subscriber to the map
        }
    }
}

void DloSimulator::changeYoungModulusCb(const std_msgs::Float32::ConstPtr msg){
    Real young_modulus = msg->data;
    dlo_.setYoungModulus(young_modulus);
}

void DloSimulator::changeTorsionModulusCb(const std_msgs::Float32::ConstPtr msg){
    Real torsion_modulus = msg->data;
    dlo_.setTorsionModulus(torsion_modulus);
}

// Create setYoungModulus service callback
bool DloSimulator::setDloYoungModulusCallback(dlo_simulator_stiff_rods::SetDloYoungModulus::Request &req, 
                                                dlo_simulator_stiff_rods::SetDloYoungModulus::Response &res) {
    dlo_.setYoungModulus(req.young_modulus);
    res.success = true;
    return true;
}

// Create setTorsionModulus service callback
bool DloSimulator::setDloTorsionModulusCallback(dlo_simulator_stiff_rods::SetDloTorsionModulus::Request &req, 
                                                dlo_simulator_stiff_rods::SetDloTorsionModulus::Response &res) {
    dlo_.setTorsionModulus(req.torsion_modulus);
    res.success = true;
    return true;
}

// Create getYoungModulus service callback
bool DloSimulator::getDloYoungModulusCallback(dlo_simulator_stiff_rods::GetDloYoungModulus::Request &req, 
                                                dlo_simulator_stiff_rods::GetDloYoungModulus::Response &res) {
    res.young_modulus = dlo_.getYoungModulus();
    return true;
}

// Create getTorsionModulus service callback
bool DloSimulator::getDloTorsionModulusCallback(dlo_simulator_stiff_rods::GetDloTorsionModulus::Request &req, 
                                                dlo_simulator_stiff_rods::GetDloTorsionModulus::Response &res) {
    res.torsion_modulus = dlo_.getTorsionModulus();
    return true;
}

// Create enableCollisionHandling service callback
bool DloSimulator::enableCollisionHandlingCallback(dlo_simulator_stiff_rods::EnableCollisionHandling::Request &req,
                                                    dlo_simulator_stiff_rods::EnableCollisionHandling::Response &res) {
    boost::recursive_mutex::scoped_lock lock(mtx_);

    is_collision_handling_enabled_ = req.is_enable;

    // log the message
    if (is_collision_handling_enabled_) {
        ROS_INFO("Collision handling is enabled.");
    } else {
        collision_handler_->resetContacts();
        ROS_INFO("Collision handling is disabled.");
    }

    res.success = true;
    return true;
}

pbd_object::MeshDLO DloSimulator::createMeshDLO(const std::string &name, 
                                                const Real &dlo_l, 
                                                const Real &dlo_z, 
                                                const int &dlo_num_segments){
    // Assume: a dlo centered at the origin, and parallel with the y-axis, 
    // Assume: particle ids start from 0, and each id is a segment
    // Each segment has an orientation specified with a quaternion

    int num_segments = dlo_num_segments;
    int num_quaternions = num_segments;
    Real segment_l = dlo_l/static_cast<Real>(num_segments);  // segment length

    // create a linear spaced coordinate vectors for the coordinates

    // Create vector of 3D Eigen vectors to hold the "list of vertices" 
    std::vector<Eigen::Matrix<Real,3,1>> vertices(num_segments);
    // Create vector of queternions to hold the "list of orientations" 
    std::vector<Eigen::Quaternion<Real>> quaternions(num_quaternions);
    // Create vector of segment lengths
    std::vector<Real> segment_lengths(num_segments);

    
    // Generate the y_coords
    Eigen::Matrix<Real,1,Eigen::Dynamic> y_coords(num_segments);
    Real start = (dlo_l - segment_l)/2.0 ;
    Real end = -(dlo_l - segment_l)/2.0;
    Real step = (end - start) / (static_cast<Real>(num_segments) - 1.0);

    for (int i = 0; i < num_segments; i++) {
        y_coords(i) = start + i * step;
    }

    // Create vertices with x,y,z coordinates (x=0 based on the assumptions above)
    for (int j = 0; j < y_coords.size(); j++) {
        Eigen::Matrix<Real,3,1> v(0.0, y_coords(j), dlo_z);
        vertices[j] = v;
    }

    // Create quaterninons (orientations) 
    Eigen::Matrix<Real, 1, 3> from(0.0, 0.0, 1.0); // dlo is expanded following z-axes
    Eigen::Matrix<Real, 1, 3> to;
    Eigen::Quaternion<Real> dq;

    for (int j = 0; j < num_quaternions; j++) {
        // std::cout << "j: " << j << std::endl;
        if (j < num_quaternions - 1) { // ie. until the last vertex position is used
            to = (vertices[j+1] - vertices[j]).normalized();
            // std::cout << "to: " << to << std::endl;

            dq = Eigen::Quaternion<Real>::FromTwoVectors(from,to);
            // std::cout << "dq: " << dq.coeffs() << std::endl;

            if (j == 0){
                quaternions[j] = dq;
            } else {
                quaternions[j] = (dq * quaternions[j-1]); // notice the reverse order!
            }
            from = to;
        } else {
            // Assume the last one has the same orientation with the previous one.
            quaternions[j] = quaternions[j-1]; 
        }
        // std::cout << "quaternions[j]: " << quaternions[j].coeffs() << std::endl;
    }

    for (int i = 0; i < num_segments; i++) {
        segment_lengths[i] = segment_l; // they have the same length 
    }

    pbd_object::MeshDLO mesh;
    mesh.name = name;
    mesh.vertices = vertices;
    mesh.quaternions = quaternions;
    mesh.segment_lengths = segment_lengths; 

    return mesh;
}

pbd_object::MeshDLO DloSimulator::transformMeshDLO(const pbd_object::MeshDLO &mesh,
                                                const std::vector<Real> &translation,
                                                const std::vector<Real> &rotationAxis,
                                                const Real &rotationAngle,
                                                const std::vector<Real> &scale){
    // Translates, rotates, and scales a given DLO mesh with Rotation given with Axis-Angle
    // Note: scale is not used. 

    // Check size of vectors
    assert(translation.size() == 3);
    assert(rotationAxis.size() == 3);
    assert(scale.size() == 3);

    // Prepare translation, rotation and scale matrices
    Eigen::Matrix<Real,3,1> translationVector(translation[0], translation[1], translation[2]);
    Eigen::AngleAxis<Real> rotationMatrix(rotationAngle, Eigen::Matrix<Real,3,1>(rotationAxis[0], rotationAxis[1], rotationAxis[2]));
    Eigen::Matrix<Real,3,1> scaleVector(scale[0], scale[1], scale[2]);

    // Create a copy of the input mesh
    pbd_object::MeshDLO transformed_mesh = mesh;

    // Convert angle-axis to quaternion
    Eigen::Quaternion<Real> rotationQuaternion(rotationMatrix);

    // Apply the transformation matrix to each vertex and quaternion in the mesh
    for(int i = 0; i < transformed_mesh.vertices.size(); i++){
        Eigen::Matrix<Real,3,1> &vertex = transformed_mesh.vertices[i];
        vertex = rotationQuaternion * vertex + translationVector;

        Eigen::Quaternion<Real> &quaternion = transformed_mesh.quaternions[i];
        // Update the quaternion here with the new orientation
        quaternion = rotationQuaternion * quaternion; // Apply rotation
        // std::cout << "quaternion[i]: " << quaternion.coeffs() << std::endl;
    }
    return transformed_mesh;
}

pbd_object::Mesh DloSimulator::loadMesh(const std::string &name, 
                                        const std::string &path, 
                                        const Real &z) {
    // Function to load an .obj file 
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file " + path);
    }

    std::vector<Eigen::Matrix<Real, 1, 3>> vertices;
    std::vector<Eigen::Matrix<int, 1, 3>> face_tri_ids;
    std::vector<Eigen::Matrix<Real, 1, 2>> tex_coords;
    std::vector<Eigen::Matrix<Real, 1, 3>> normals;

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string prefix;
        ss >> prefix;

        if (prefix == "v") {
            Eigen::Matrix<Real, 1, 3> vertex;
            ss >> vertex(0, 0) >> vertex(0, 1) >> vertex(0, 2);
            vertex(0, 2) += z; // Adjust z-coordinate
            vertices.push_back(vertex);
        } else if (prefix == "vt") {
            Eigen::Matrix<Real, 1, 2> tex_coord;
            ss >> tex_coord(0, 0) >> tex_coord(0, 1);
            tex_coords.push_back(tex_coord);
        } else if (prefix == "vn") {
            Eigen::Matrix<Real, 1, 3> normal;
            ss >> normal(0, 0) >> normal(0, 1) >> normal(0, 2);
            normals.push_back(normal);
        } else if (prefix == "f") {
            Eigen::Matrix<int, 1, 3> face;
            for (int i = 0; i < 3; ++i) {
                std::string vertex_spec;
                ss >> vertex_spec;
                std::stringstream vertex_ss(vertex_spec);
                std::string vertex_index;
                std::getline(vertex_ss, vertex_index, '/');
                face(0, i) = std::stoi(vertex_index) - 1; // Convert to 0-based index
            }
            face_tri_ids.push_back(face);
        }
    }
    file.close();

    // Convert vectors to Eigen matrices
    Eigen::Matrix<Real, Eigen::Dynamic, 3> vertices_mat(vertices.size(), 3);
    for (size_t i = 0; i < vertices.size(); ++i) {
        vertices_mat.row(i) = vertices[i];
    }

    Eigen::Matrix<int, Eigen::Dynamic, 3> face_tri_ids_mat(face_tri_ids.size(), 3);
    for (size_t i = 0; i < face_tri_ids.size(); ++i) {
        face_tri_ids_mat.row(i) = face_tri_ids[i];
    }

    Eigen::Matrix<Real, Eigen::Dynamic, 2> tex_coords_mat(tex_coords.size(), 2);
    for (size_t i = 0; i < tex_coords.size(); ++i) {
        tex_coords_mat.row(i) = tex_coords[i];
    }

    Eigen::Matrix<Real, Eigen::Dynamic, 3> normals_mat(normals.size(), 3);
    for (size_t i = 0; i < normals.size(); ++i) {
        normals_mat.row(i) = normals[i];
    }

    pbd_object::Mesh mesh;
    mesh.name = name;
    mesh.vertices = vertices_mat;
    mesh.face_tri_ids = face_tri_ids_mat;
    mesh.tex_coords = tex_coords_mat;
    mesh.normals = normals_mat;

    return mesh;
}


pbd_object::Mesh DloSimulator::transformMesh(const pbd_object::Mesh &mesh, 
                                   const Eigen::Matrix<Real, 1, 3> &translation,
                                   const Eigen::Quaternion<Real> &rotation,
                                   const Eigen::Matrix<Real, 1, 3> &scale){
    // Translates, rotates, and scales a given mesh with Rotation given with Quaternion 

    // Prepare translation and scale matrices
    Eigen::Matrix<Real,3,1> translationVector(translation[0], translation[1], translation[2]);
    Eigen::Matrix<Real,3,1> scaleVector(scale[0], scale[1], scale[2]);

    // Create a copy of the input mesh
    pbd_object::Mesh transformed_mesh = mesh;

    // Apply the transformation matrix to each vertex in the mesh
    for(int i=0; i<transformed_mesh.vertices.rows(); i++)
    {
        Eigen::Matrix<Real,3,1> vertex = transformed_mesh.vertices.row(i).transpose();
        vertex = rotation.normalized() * (vertex.cwiseProduct(scaleVector)) + translationVector;
        transformed_mesh.vertices.row(i) = vertex.transpose();
    }

    return transformed_mesh;
}

void DloSimulator::simulate(const ros::TimerEvent& e){
    // With some kind of self lock to prevent collision with rendering
    boost::recursive_mutex::scoped_lock lock(mtx_);

    // // Attempt to acquire the lock without blocking
    // boost::unique_lock<boost::recursive_mutex> lock(mtx_, boost::try_to_lock);

    // if (!lock.owns_lock()) {
    //     // Could not acquire lock, previous call is still running
    //     ROS_WARN("simulate function skipped: previous call still in progress.");
    //     return;
    // }

    Real sdt = dt_ / num_substeps_;

    ros::Time start_time = ros::Time::now();
    
    // --------------------------------------------------------------
    // Small steps (Substep XPBD) 2019 implementation

    // Real change_in_max_error = std::numeric_limits<Real>::infinity(); // Init change in error to infinity
    // Real prev_max_error = 0.0;
    // Real current_max_error = 0.0;

    for (int i = 0; i< num_steps_; i++){

        dlo_.resetForces();
        dlo_.resetLambdas();

        int j;
        for (j = 0; j < num_substeps_; j++){
            dlo_.preSolve(sdt,gravity_);

            // Collision Handling, detect collisions
            if (is_collision_handling_enabled_){
                collision_handler_->collisionDetection();  
            }
            collision_handler_->solveContactPositionConstraints(sdt);
            dlo_.solve(sdt);

            dlo_.postSolve(sdt);

            collision_handler_->solveContactVelocityConstraints(sdt);

            // current_max_error = dlo_.getMaxError();
            // if (current_max_error < 1.0e-2){
            //     i = num_steps_;
            //     break;
            // }
            // // Calculate change in max error
            // change_in_max_error = std::abs(current_max_error - prev_max_error);

            // if (change_in_max_error < 1.0e-9){
            //     i = num_steps_;
            //     break;
            // // }
            // // Update prev_max_error
            // prev_max_error = current_max_error;

        }
        // std::cout << "itr: " << j << ",Max error: " << current_max_error << std::endl;
        // std::cout << "itr: " << j << ",Max error change: " << change_in_max_error << std::endl;

    }
    // --------------------------------------------------------------


    // // --------------------------------------------------------------
    // // XPBD 2016 implementation
    // dlo_.resetForces();
    // dlo_.resetLambdas();

    // dlo_.preSolve(dt_,gravity_);
    
    // for (int i = 0; i< num_steps_; i++){
    //     int j;
    //     for (j = 0; j < num_substeps_; j++){
    //         // Collision Handling, detect collisions
    //         if (is_collision_handling_enabled_){
    //             collision_handler_->collisionDetection();  
    //         }
    //         collision_handler_->solveContactPositionConstraints(sdt);
    //         dlo_.solve(sdt);
    //     }
    // }

    // dlo_.postSolve(dt_);
    // collision_handler_->solveContactVelocityConstraints(dt_);
    // // --------------------------------------------------------------


    // // To debug force readings from hanged corners (use only when robots are not attached)
    // std::vector<int> *attached_ids_ptr = dlo_.getAttachedIdsPtr();
    // Eigen::Matrix<Real,Eigen::Dynamic,3> *for_ptr = dlo_.getForPtr();
    // int id = (*attached_ids_ptr)[0]; // First attached id
    // // force at that attached id
    // std::cout << "id: " << id << ". Force = " << for_ptr->col(id)/num_substeps_ << " N." << std::endl;

    // +++++++++++++++++++++++++++++++++++++++++++++++++++==
    // readAttachedRobotForces();
    dlo_.normalizeForces(num_substeps_);
    // +++++++++++++++++++++++++++++++++++++++++++++++++++==

    ros::Time finish_time = ros::Time::now();
    ros::Duration elapsed_time = finish_time - start_time;
    time_sum_ += elapsed_time.toSec();
    time_frames_ += 1;
    if (time_frames_ > 2) {
        time_sum_ /= time_frames_;
        // ROS_INFO("[Dlo Simulator]: %-4.2lf ms per simulation iteration", time_sum_*1000);
        ROS_INFO_THROTTLE(5, "[Dlo Simulator]: %-4.2lf ms per simulation iteration", time_sum_*1000);
        // Smart dt and simulation rate selection (debug)
        if (!is_auto_sim_rate_set_ && set_sim_rate_auto_) {
            dt_ = time_sum_;
            timer_simulate_.setPeriod(ros::Duration(time_sum_));
            is_auto_sim_rate_set_ = true; 
            std::cout << "Automatic setup. dt = " << dt_ << " seconds." << std::endl;
        }
        time_frames_ = 0;
        time_sum_ = 0;
    }
}

void DloSimulator::render(const ros::TimerEvent& e){
    if (dlo_visualization_mode_ == -1){
        return;
    }
    // // With some kind of self lock to prevent collision with simulation
    boost::recursive_mutex::scoped_lock lock(mtx_);

    visualization_msgs::MarkerArray markerArray;
    int marker_id = 0;
    // int dlo_visualization_mode_: 
    // 0: Points Only, 
    // 1: Line Segments Only, 
    // 2: Orientation Frames and Line Segments, 
    // 3: Points and Line Segments, 
    // 4: All

    const std::vector<Eigen::Matrix<Real,3,1>> *pos_ptr = dlo_.getPosPtr();
    const std::vector<Eigen::Quaternion<Real>> *ori_ptr = dlo_.getOriPtr();
    const std::vector<Real> *len_ptr = dlo_.getSegmentLengthsPtr();

    visualization_msgs::Marker pointsMarker, linesMarker, framesMarker;
    createRvizPointsMarker(pos_ptr, pointsMarker);
    createRvizLinesMarker(pos_ptr, ori_ptr, len_ptr, linesMarker);
    createRvizFramesMarker(pos_ptr, ori_ptr, framesMarker); 
    
    if (dlo_visualization_mode_ == 0 || 
        dlo_visualization_mode_ == 3 || 
        dlo_visualization_mode_ == 4) {
        // Render segment center points
        pointsMarker.id = marker_id++;
        markerArray.markers.push_back(pointsMarker);
    }
    if (dlo_visualization_mode_ == 1 || 
        dlo_visualization_mode_ == 2 || 
        dlo_visualization_mode_ == 3 || 
        dlo_visualization_mode_ == 4) {
        // Render each line segment seperately based on their orientation and length
        linesMarker.id = marker_id++;
        markerArray.markers.push_back(linesMarker);
    }
    if (dlo_visualization_mode_ == 2 || 
        dlo_visualization_mode_ == 4) {
        // Render each orientation frame
        framesMarker.id = marker_id++;
        markerArray.markers.push_back(framesMarker);
    }

    pub_dlo_marker_array_.publish(markerArray);
}

void DloSimulator::createRvizPointsMarker(const std::vector<Eigen::Matrix<Real,3,1>> *poses, 
                                          visualization_msgs::Marker &marker){
    marker.header.frame_id = dlo_frame_id_;
    marker.header.stamp = ros::Time::now();
    marker.type = visualization_msgs::Marker::POINTS;
    marker.action = visualization_msgs::Marker::ADD;
    marker.pose.orientation.w = 1.0;

    marker.scale.x = point_marker_scale_;
    marker.scale.y = point_marker_scale_;
    marker.scale.z = point_marker_scale_;

    marker.color.r = point_marker_color_rgba_[0];
    marker.color.g = point_marker_color_rgba_[1];
    marker.color.b = point_marker_color_rgba_[2];
    marker.color.a = point_marker_color_rgba_[3];

    for (int i = 0; i < poses->size(); i++) {
        geometry_msgs::Point p;
        p.x = (*poses)[i](0);
        p.y = (*poses)[i](1);
        p.z = (*poses)[i](2);
        marker.points.push_back(p);
    }
}

void DloSimulator::createRvizLinesMarker(const std::vector<Eigen::Matrix<Real,3,1>> *poses,
                                         const std::vector<Eigen::Quaternion<Real>> *orients,
                                         const std::vector<Real> *lengths,
                                         visualization_msgs::Marker &marker){
    marker.header.frame_id = dlo_frame_id_;
    marker.header.stamp = ros::Time::now();
    marker.type = visualization_msgs::Marker::LINE_LIST;
    marker.action = visualization_msgs::Marker::ADD;
    marker.pose.orientation.w = 1.0;

    // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
    marker.scale.x = 0.005*line_marker_scale_multiplier_;

    marker.color.r = line_marker_color_rgba_[0];
    marker.color.g = line_marker_color_rgba_[1];
    marker.color.b = line_marker_color_rgba_[2];
    marker.color.a = line_marker_color_rgba_[3];

    for (int i = 0; i < poses->size(); i++){
        // tip vector of the segment in segment frame
        Eigen::Matrix<Real,3,1>v(0, 0, 0.5*(*lengths)[i]);

        // tip vector of the segment rotated to world frame
        v = (*orients)[i] * v;

        // tips of segment in world frame
        const Eigen::Matrix<Real,3,1> &p = (*poses)[i];
        const Eigen::Matrix<Real,3,1> p1w = p + v;
        const Eigen::Matrix<Real,3,1> p2w = p - v;

        geometry_msgs::Point p1;
        p1.x = p1w(0);
        p1.y = p1w(1);
        p1.z = p1w(2);
        marker.points.push_back(p1);

        geometry_msgs::Point p2;
        p2.x = p2w(0);
        p2.y = p2w(1);
        p2.z = p2w(2);
        marker.points.push_back(p2);
    }
}

void DloSimulator::createRvizFramesMarker(const std::vector<Eigen::Matrix<Real,3,1>> *poses,
                                         const std::vector<Eigen::Quaternion<Real>> *orients,
                                         visualization_msgs::Marker &marker){
    marker.header.frame_id = dlo_frame_id_;
    marker.header.stamp = ros::Time::now();

    //frame_marker_scale_
    // Define the marker properties for frame visualization
    marker.type = visualization_msgs::Marker::LINE_LIST;
    marker.action = visualization_msgs::Marker::ADD;
    marker.pose.orientation.w = 1.0;

    // LINE_LIST markers use only the x component of scale, for the line width
    marker.scale.x = frame_marker_scale_;  // Assuming frame_marker_scale_ is defined elsewhere

    for (size_t i = 0; i < poses->size(); ++i) {
        // Define the base of the frame
        const Eigen::Matrix<Real,3,1> &base = (*poses)[i];

        // Calculate the end points for each axis of the frame
        Eigen::Matrix<Real,3,1> x_axis(1.0, 0.0, 0.0);
        Eigen::Matrix<Real,3,1> y_axis(0.0, 1.0, 0.0);
        Eigen::Matrix<Real,3,1> z_axis(0.0, 0.0, 1.0);

        x_axis = (*orients)[i] * x_axis * frame_marker_axis_length_; 
        y_axis = (*orients)[i] * y_axis * frame_marker_axis_length_;
        z_axis = (*orients)[i] * z_axis * frame_marker_axis_length_;

        // Add lines for X-axis in red
        addAxisLineToMarker(marker, base, base + x_axis, 1.0, 0.0, 0.0, 1.0);
        // Add lines for Y-axis in green
        addAxisLineToMarker(marker, base, base + y_axis, 0.0, 1.0, 0.0, 1.0);
        // Add lines for Z-axis in blue
        addAxisLineToMarker(marker, base, base + z_axis, 0.0, 0.0, 1.0, 1.0);
    }   
}

void DloSimulator::addAxisLineToMarker(visualization_msgs::Marker &marker, 
                                       const Eigen::Matrix<Real,3,1> &start,
                                       const Eigen::Matrix<Real,3,1> &end,
                                       float r, float g, float b, float a) {
    geometry_msgs::Point p1, p2;
    std_msgs::ColorRGBA color;

    p1.x = start(0);
    p1.y = start(1);
    p1.z = start(2);

    p2.x = end(0);
    p2.y = end(1);
    p2.z = end(2);

    color.r = r;
    color.g = g;
    color.b = b;
    color.a = a;

    // Add the start and end points twice, as each line needs start and end colors
    marker.points.push_back(p1);
    marker.points.push_back(p2);
    marker.colors.push_back(color);
    marker.colors.push_back(color);
}

// Publish message to RVIZ to visualize the rigid bodies considered in the simulation
void DloSimulator::renderRigidBodies(const ros::TimerEvent& e){
    if (rb_visualization_mode_ == -1){
        return;
    }

    visualization_msgs::MarkerArray markerArray;
    marker_id_ = 0; 

    // int rb_visualization_mode_ = 0; // 0: Mesh Only, 1: Wireframe Only, 2: Both

    for (const utilities::RigidBodySceneLoader::RigidBodyData& rbd : rigid_bodies_) {
        visualization_msgs::Marker meshMarker, wireframeMarker;
        createMeshAndWireframeMarkers(&rbd.m_mesh.vertices, &rbd.m_mesh.face_tri_ids, meshMarker, wireframeMarker);

        if (rb_visualization_mode_ == 0 || rb_visualization_mode_ == 2) {
            meshMarker.id = marker_id_++;
            markerArray.markers.push_back(meshMarker);
        }
        if (rb_visualization_mode_ == 1 || rb_visualization_mode_ == 2) {
            wireframeMarker.id = marker_id_++;
            markerArray.markers.push_back(wireframeMarker);
        }
    }

    // Publish the entire marker array
    pub_rb_marker_array_.publish(markerArray);

}

void DloSimulator::createMeshAndWireframeMarkers(const Eigen::Matrix<Real,Eigen::Dynamic,3> *vertices, 
                                                    const Eigen::MatrixX3i *face_tri_ids, 
                                                    visualization_msgs::Marker &meshMarker, 
                                                    visualization_msgs::Marker &wireframeMarker){
    // Set up mesh marker (triangles)
    setupMeshMarker(meshMarker);
    // Set up wireframe marker (lines)
    setupWireframeMarker(wireframeMarker);

    for (int i = 0; i < face_tri_ids->rows(); i++) {
        for (int j = 0; j < 3; j++) {
            int idx = (*face_tri_ids)(i, j);
            geometry_msgs::Point p;
            p.x = (*vertices)(idx, 0);
            p.y = (*vertices)(idx, 1);
            p.z = (*vertices)(idx, 2);
            meshMarker.points.push_back(p);

            // Wireframe points
            int next_idx = (*face_tri_ids)(i, (j + 1) % 3);
            geometry_msgs::Point p2;
            p2.x = (*vertices)(next_idx, 0);
            p2.y = (*vertices)(next_idx, 1);
            p2.z = (*vertices)(next_idx, 2);
            wireframeMarker.points.push_back(p);
            wireframeMarker.points.push_back(p2);
        }
    }
}

void DloSimulator::setupMeshMarker(visualization_msgs::Marker &marker){
    marker.header.frame_id = dlo_frame_id_;
    marker.header.stamp = ros::Time::now();

    marker.type = visualization_msgs::Marker::TRIANGLE_LIST;
    // marker.id = 1; ID IS SET IN renderRigidBodies FUNCTION
    marker.action = visualization_msgs::Marker::ADD;
    
    marker.pose.orientation.w = 1.0;
    
    // marker.points = points; POINTS ARE PUSHED IN createMeshAndWireframeMarkers FUNCTION
    
    // Set other properties like color, scale, etc.
    marker.scale.x = 1.0; // As it's a list of triangles, scale shouldn't matter
    marker.scale.y = 1.0;
    marker.scale.z = 1.0;

    marker.color.r = rb_mesh_marker_color_rgba_[0]; 
    marker.color.g = rb_mesh_marker_color_rgba_[1];
    marker.color.b = rb_mesh_marker_color_rgba_[2];
    marker.color.a = rb_mesh_marker_color_rgba_[3];
}

void DloSimulator::setupWireframeMarker(visualization_msgs::Marker &marker){
    marker.header.frame_id = dlo_frame_id_;
    marker.header.stamp = ros::Time::now();

    marker.type = visualization_msgs::Marker::LINE_LIST;
    // marker.id = 1; ID IS SET IN renderRigidBodies FUNCTION
    marker.action = visualization_msgs::Marker::ADD;
    
    marker.pose.orientation.w = 1.0;
    
    // marker.points = points; POINTS ARE PUSHED IN createMeshAndWireframeMarkers FUNCTION

    // Set other properties like color, scale, etc.
    // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
    marker.scale.x = 0.005*rb_line_marker_scale_multiplier_;

    marker.color.r = rb_line_marker_color_rgba_[0];
    marker.color.g = rb_line_marker_color_rgba_[1];
    marker.color.b = rb_line_marker_color_rgba_[2];
    marker.color.a = rb_line_marker_color_rgba_[3];
}

void DloSimulator::odometryCb_custom_static_particles(const nav_msgs::Odometry::ConstPtr& odom_msg, const int& id) {
    // // With some kind of self lock to prevent collision with simulation
    boost::recursive_mutex::scoped_lock lock(mtx_);

    // Update Velocity
    const Eigen::Matrix<Real,3,1> vel(odom_msg->twist.twist.linear.x, 
                                        odom_msg->twist.twist.linear.y, 
                                        odom_msg->twist.twist.linear.z);

    const Eigen::Matrix<Real,3,1> omega(odom_msg->twist.twist.angular.x, 
                                        odom_msg->twist.twist.angular.y, 
                                        odom_msg->twist.twist.angular.z);


    // Update pose
    const Real x = odom_msg->pose.pose.position.x;
    const Real y = odom_msg->pose.pose.position.y;
    const Real z = odom_msg->pose.pose.position.z;

    const Real qw = odom_msg->pose.pose.orientation.w;
    const Real qx = odom_msg->pose.pose.orientation.x;
    const Real qy = odom_msg->pose.pose.orientation.y;
    const Real qz = odom_msg->pose.pose.orientation.z;
    
    Eigen::Matrix<Real,3,1> pos(x, y, z);
    Eigen::Quaternion<Real> ori(qw,qx,qy,qz);

    dlo_.updateAttachedVelocity(id, vel, omega);
    dlo_.updateAttachedPose(id, pos, ori);
}

void DloSimulator::cmdVelCb_custom_static_particles(const geometry_msgs::Twist::ConstPtr& twist_msg, const int& id) {
    // // With some kind of self lock to prevent collision with simulation
    boost::recursive_mutex::scoped_lock lock(mtx_);

    // Update Velocity
    const Eigen::Matrix<Real,3,1> vel(twist_msg->linear.x, 
                                      twist_msg->linear.y, 
                                      twist_msg->linear.z);

    const Eigen::Matrix<Real,3,1> omega(twist_msg->angular.x, 
                                        twist_msg->angular.y, 
                                        twist_msg->angular.z);

    dlo_.updateAttachedVelocity(id, vel, omega);
}


/*
void DloSimulator::readAttachedRobotForces(){
    Eigen::Matrix<Real,Eigen::Dynamic,3> *for_ptr = dlo_.getForPtr();

    // std::cout << "Here" << std::endl;

    if (is_rob_01_attached_){
        rob_01_attached_force_ = for_ptr->col(rob_01_attached_id_)/num_substeps_;
        // std::cout << "id: " << rob_01_attached_id_ << ". Force = " << ": " << for_ptr->col(rob_01_attached_id_)/num_substeps_ << " N." << std::endl;
    }
    if (is_rob_02_attached_){
        rob_02_attached_force_ = for_ptr->col(rob_02_attached_id_)/num_substeps_;
    }
    if (is_rob_03_attached_){
        rob_03_attached_force_ = for_ptr->col(rob_03_attached_id_)/num_substeps_;
    }
    if (is_rob_04_attached_){
        rob_04_attached_force_ = for_ptr->col(rob_04_attached_id_)/num_substeps_;
    }
}

// Publish forces on each robot by the dlo
void DloSimulator::publishWrenches(const ros::TimerEvent& e){
    geometry_msgs::WrenchStamped msg;
    msg.header.stamp = ros::Time::now();

    // std::cout << "Here22" << std::endl;

    if (is_rob_01_attached_){
        // std::cout << "id: " << rob_01_attached_id_ << ". Force = " << rob_01_attached_force_ << " N." << std::endl;

        msg.header.frame_id = wrench_01_frame_id_;
        msg.wrench.force.x = rob_01_attached_force_(0);
        msg.wrench.force.y = rob_01_attached_force_(1);
        msg.wrench.force.z = rob_01_attached_force_(2);
        msg.wrench.torque.x = 0.0;
        msg.wrench.torque.y = 0.0;
        msg.wrench.torque.z = 0.0;
        pub_wrench_stamped_01_.publish(msg);
    }
    if (is_rob_02_attached_){
        // std::cout << "id: " << rob_02_attached_id_ << ". Force = " << rob_02_attached_force_ << " N." << std::endl;

        msg.header.frame_id = wrench_02_frame_id_;
        msg.wrench.force.x = rob_02_attached_force_(0);
        msg.wrench.force.y = rob_02_attached_force_(1);
        msg.wrench.force.z = rob_02_attached_force_(2);
        msg.wrench.torque.x = 0.0;
        msg.wrench.torque.y = 0.0;
        msg.wrench.torque.z = 0.0;
        pub_wrench_stamped_02_.publish(msg);
    }
    if (is_rob_03_attached_){
        // std::cout << "id: " << rob_03_attached_id_ << ". Force = " << rob_03_attached_force_ << " N." << std::endl;

        msg.header.frame_id = wrench_03_frame_id_;
        msg.wrench.force.x = rob_03_attached_force_(0);
        msg.wrench.force.y = rob_03_attached_force_(1);
        msg.wrench.force.z = rob_03_attached_force_(2);
        msg.wrench.torque.x = 0.0;
        msg.wrench.torque.y = 0.0;
        msg.wrench.torque.z = 0.0;
        pub_wrench_stamped_03_.publish(msg);
    }
    if (is_rob_04_attached_){
        // std::cout << "id: " << rob_04_attached_id_ << ". Force = " << rob_04_attached_force_ << " N." << std::endl;

        msg.header.frame_id = wrench_04_frame_id_;
        msg.wrench.force.x = rob_04_attached_force_(0);
        msg.wrench.force.y = rob_04_attached_force_(1);
        msg.wrench.force.z = rob_04_attached_force_(2);
        msg.wrench.torque.x = 0.0;
        msg.wrench.torque.y = 0.0;
        msg.wrench.torque.z = 0.0;
        pub_wrench_stamped_04_.publish(msg);
    }    
}
*/

void DloSimulator::contactCallbackFunction(const unsigned int contactType,
                                              const unsigned int bodyIndex1, 
                                              const unsigned int bodyIndex2,
                                              const Eigen::Matrix<Real, 3, 1> &cp1, 
                                              const Eigen::Matrix<Real, 3, 1> &cp2,
                                              const Eigen::Matrix<Real, 3, 1> &normal, 
                                              const Real dist,
                                              const Real restitutionCoeff, 
                                              const Real frictionCoeffStatic,
                                              const Real frictionCoeffDynamic,
                                              void *userData)
{
    DloSimulator *dloSimulator = (DloSimulator*)userData; 

	if (contactType == utilities::CollisionHandler::RigidBodyContactType)
    {
        dloSimulator->collision_handler_->addRigidBodyContactConstraint(bodyIndex1, 
                                                          bodyIndex2, 
                                                          cp1,
                                                          cp2,
                                                          normal, 
                                                          dist,
                                                          restitutionCoeff,
                                                          frictionCoeffStatic,
                                                          frictionCoeffDynamic);
    }
	else if (contactType == utilities::CollisionHandler::ParticleRigidBodyContactType)
    {
        dloSimulator->collision_handler_->addParticleRigidBodyContactConstraint(bodyIndex1, 
                                                                  bodyIndex2, 
                                                                  cp1, 
                                                                  cp2, 
                                                                  normal, 
                                                                  dist, 
                                                                  restitutionCoeff, 
                                                                  frictionCoeffStatic,
                                                                  frictionCoeffDynamic);
    }
		
}


// Publish the minimum distances to rigid bodies in the simulation message
void DloSimulator::publishMinDistancesToRigidBodies(const ros::TimerEvent& e){
    // With some kind of self lock to prevent collision with rendering
    boost::recursive_mutex::scoped_lock lock(mtx_);
    
    std::vector<std::vector<utilities::CollisionHandler::MinDistanceData>> min_distances_mt;
    collision_handler_->computeMinDistancesToRigidBodies(min_distances_mt);

    // // Print min distance data for debug
    // // Iterate through all threads
    // for (size_t threadIdx = 0; threadIdx < min_distances_mt.size(); ++threadIdx) {
    //     // std::cout << "Thread " << threadIdx << ":" << std::endl;        
    //     // Iterate through all minimum distance data in the current thread
    //     for (size_t i = 0; i < min_distances_mt[threadIdx].size(); ++i) {
    //         const auto& minDistData = min_distances_mt[threadIdx][i];
    //         std::cout << "  Min Distance Data " << i << ":" << std::endl;
    //         // std::cout << "    Type: " << static_cast<int>(minDistData.m_type) << std::endl;
    //         // std::cout << "    Index on Object 1: " << minDistData.m_index1 << std::endl;
    //         // std::cout << "    Index on Object 2: " << minDistData.m_index2 << std::endl;
    //         std::cout << "    Minimum Distance: " << minDistData.m_minDistance << std::endl;
    //         // std::cout << "    Point on Object 1: (" 
    //                 //   << minDistData.m_pointOnObject1.x() << ", "
    //                 //   << minDistData.m_pointOnObject1.y() << ", "
    //                 //   << minDistData.m_pointOnObject1.z() << ")" << std::endl;
    //         // std::cout << "    Point on Object 2: (" 
    //                 //   << minDistData.m_pointOnObject2.x() << ", "
    //                 //   << minDistData.m_pointOnObject2.y() << ", "
    //                 //   << minDistData.m_pointOnObject2.z() << ")" << std::endl;
    //     }
    // }

    //Create a publisher for the min_distances data
    dlo_simulator_stiff_rods::MinDistanceDataArray min_distances_msg;

    for (const auto& thread_data : min_distances_mt) {
        for (const auto& minDistData : thread_data) {
            dlo_simulator_stiff_rods::MinDistanceData msg_data;

            // Fill the message fields
            msg_data.header.stamp = ros::Time::now();

            msg_data.type = minDistData.m_type;
            
            // msg_data.index1 = minDistData.m_index1; // index of the particle in dlo
            msg_data.index1 = 0; // index of the dlo in the scene, when there is single dlo in the simulation, it's always 0.  
            msg_data.index2 = minDistData.m_index2;

            msg_data.pointOnObject1.x = minDistData.m_pointOnObject1[0];
            msg_data.pointOnObject1.y = minDistData.m_pointOnObject1[1];
            msg_data.pointOnObject1.z = minDistData.m_pointOnObject1[2];

            msg_data.pointOnObject2.x = minDistData.m_pointOnObject2[0];
            msg_data.pointOnObject2.y = minDistData.m_pointOnObject2[1];
            msg_data.pointOnObject2.z = minDistData.m_pointOnObject2[2];

            // For the normal vector
            msg_data.normal.x = minDistData.m_normal.x();
            msg_data.normal.y = minDistData.m_normal.y();
            msg_data.normal.z = minDistData.m_normal.z();

            msg_data.minDistance = minDistData.m_minDistance;

            min_distances_msg.data.push_back(msg_data);
        }
    }

    pub_min_dist_to_rb_.publish(min_distances_msg);

    //-------------------------------------------------------------------------
    // Visualize the min distance line segments if desired:
    if (visualize_min_distances_) {
        publishMinDistLineMarkers(min_distances_mt);
    }
}

void DloSimulator::publishMinDistLineMarkers(
    const std::vector<std::vector<utilities::CollisionHandler::MinDistanceData>>& min_distances_mt) {

    visualization_msgs::MarkerArray marker_array;
    visualization_msgs::Marker line_marker;
    setupMinDistLineMarker(line_marker);

    int marker_id = 0;
    for (const auto& thread_data : min_distances_mt) {
        for (const auto& minDistData : thread_data) {
            line_marker.id = marker_id++;

            geometry_msgs::Point p1, p2;
            p1.x = minDistData.m_pointOnObject1[0];
            p1.y = minDistData.m_pointOnObject1[1];
            p1.z = minDistData.m_pointOnObject1[2];

            p2.x = minDistData.m_pointOnObject2[0];
            p2.y = minDistData.m_pointOnObject2[1];
            p2.z = minDistData.m_pointOnObject2[2];

            line_marker.points.push_back(p1);
            line_marker.points.push_back(p2);

            marker_array.markers.push_back(line_marker);
            line_marker.points.clear();
        }
    }

    pub_min_dist_marker_array_.publish(marker_array);
}

void DloSimulator::setupMinDistLineMarker(visualization_msgs::Marker &marker){
    marker.header.frame_id = dlo_frame_id_;
    marker.header.stamp = ros::Time::now();

    marker.type = visualization_msgs::Marker::LINE_LIST;
    marker.action = visualization_msgs::Marker::ADD;

    marker.pose.orientation.w = 1.0;

    marker.scale.x = 0.005 * min_dist_line_marker_scale_multiplier_;

    marker.color.r = min_dist_line_marker_color_rgba_[0];
    marker.color.g = min_dist_line_marker_color_rgba_[1];
    marker.color.b = min_dist_line_marker_color_rgba_[2];
    marker.color.a = min_dist_line_marker_color_rgba_[3];
}

void DloSimulator::publishDloState(const ros::TimerEvent& e){
    // With some kind of self lock to prevent collision with simulation
    boost::recursive_mutex::scoped_lock lock(mtx_);

    dlo_simulator_stiff_rods::SegmentStateArray states_msg;

    // Set the dlo_id 
    states_msg.dlo_id = 0;  // TODO: Replace in future with dlo_id obtaining method

    const std::vector<Eigen::Matrix<Real,3,1>> &pos_ptr = *dlo_.getPosPtr();
    const std::vector<Eigen::Quaternion<Real>> &ori_ptr = *dlo_.getOriPtr();
    const Eigen::Matrix<Real,3,Eigen::Dynamic> &vel_ptr = *dlo_.getVelPtr();
    const Eigen::Matrix<Real,3,Eigen::Dynamic> &ang_vel_ptr = *dlo_.getAngVelPtr();
    const Eigen::Matrix<Real,3,Eigen::Dynamic> &for_ptr = *dlo_.getForPtr();
    const Eigen::Matrix<Real,3,Eigen::Dynamic> &tor_ptr = *dlo_.getTorPtr();

    for (size_t i = 0; i < pos_ptr.size(); ++i) {
        dlo_simulator_stiff_rods::SegmentState segment_state;

        // Set the segment id
        segment_state.id = i;

        // Set the header
        segment_state.header.stamp = ros::Time::now();
        segment_state.header.frame_id = dlo_frame_id_;

        // Set pose position
        segment_state.pose.position.x = pos_ptr[i](0);
        segment_state.pose.position.y = pos_ptr[i](1);
        segment_state.pose.position.z = pos_ptr[i](2);

        // Set pose orientation
        Eigen::Quaternion<Real> quat = ori_ptr[i];
        segment_state.pose.orientation.x = quat.x();
        segment_state.pose.orientation.y = quat.y();
        segment_state.pose.orientation.z = quat.z();
        segment_state.pose.orientation.w = quat.w();

        // Set twist (linear and angular velocity)
        segment_state.twist.linear.x = vel_ptr(0, i);
        segment_state.twist.linear.y = vel_ptr(1, i);
        segment_state.twist.linear.z = vel_ptr(2, i);

        segment_state.twist.angular.x = ang_vel_ptr(0, i);
        segment_state.twist.angular.y = ang_vel_ptr(1, i);
        segment_state.twist.angular.z = ang_vel_ptr(2, i);

        // Set wrench (force and torque)
        segment_state.wrench.force.x = for_ptr(0, i);
        segment_state.wrench.force.y = for_ptr(1, i);
        segment_state.wrench.force.z = for_ptr(2, i);

        segment_state.wrench.torque.x = tor_ptr(0, i);
        segment_state.wrench.torque.y = tor_ptr(1, i);
        segment_state.wrench.torque.z = tor_ptr(2, i);

        // Add the segment state to the message
        states_msg.states.push_back(segment_state);
    }

    // Publish the message
    pub_dlo_state_.publish(states_msg);
}