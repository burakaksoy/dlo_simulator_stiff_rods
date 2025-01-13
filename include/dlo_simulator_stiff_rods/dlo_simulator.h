/*
 * Author: Burak Aksoy
 */

#ifndef DLO_SIMULATOR_H
#define DLO_SIMULATOR_H

#include <ros/ros.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/Vector3.h>
#include <geometry_msgs/WrenchStamped.h>
#include <geometry_msgs/Twist.h>
#include <std_msgs/Float32.h>
#include <std_msgs/Int32MultiArray.h>
#include <std_msgs/Header.h>
#include <nav_msgs/Odometry.h>

#include <dlo_simulator_stiff_rods/MinDistanceDataArray.h> 
#include <dlo_simulator_stiff_rods/ChangeParticleDynamicity.h> 
#include <dlo_simulator_stiff_rods/SegmentStateArray.h> 

#include <std_srvs/Empty.h>
#include <dlo_simulator_stiff_rods/SetParticleDynamicity.h> // service
#include <dlo_simulator_stiff_rods/SetDloYoungModulus.h> // service
#include <dlo_simulator_stiff_rods/GetDloYoungModulus.h> // service
#include <dlo_simulator_stiff_rods/SetDloTorsionModulus.h> // service
#include <dlo_simulator_stiff_rods/GetDloTorsionModulus.h> // service
#include <dlo_simulator_stiff_rods/EnableCollisionHandling.h> // service

#include <time.h>
#include <math.h>
#include <numeric>
#include <vector>
#include <chrono>

#include <Eigen/Dense>
#include <Eigen/Geometry>
// #include <scipy/spatial.h>

// #include <thread>
// #include <mutex>
#include <boost/thread/recursive_mutex.hpp>

#include <string>
#include <iostream>

#include <fstream>
#include <sstream>

#include "dlo_simulator_stiff_rods/utilities/dlo.h"
#include "dlo_simulator_stiff_rods/utilities/rigid_body_scene_loader.h"
#include "dlo_simulator_stiff_rods/utilities/collision_handler.h"

#include <memory> // needed for std::unique_ptr and std::make_unique

namespace dlo_simulator
{

class DloSimulator
{
public:
    DloSimulator(ros::NodeHandle &nh, ros::NodeHandle &nh_local, boost::recursive_mutex &mtx);
    ~DloSimulator();

    pbd_object::Dlo& getDlo() { return dlo_; }
    std::vector<utilities::RigidBodySceneLoader::RigidBodyData>& getRigidBodies() { return rigid_bodies_; }
private:
    // Functions ---------------------------------
    
    pbd_object::MeshDLO createMeshDLO(const std::string &name, 
                                    const Real &dlo_l, 
                                    const Real &dlo_z, 
                                    const int &dlo_num_segments);

    pbd_object::MeshDLO transformMeshDLO(const pbd_object::MeshDLO &mesh,
                                        const std::vector<Real> &translation,
                                        const std::vector<Real> &rotationAxis,
                                        const Real &rotationAngle,
                                        const std::vector<Real> &scale);
                                        
    pbd_object::Mesh loadMesh(const std::string &name, const std::string &path, const Real &z);

    pbd_object::Mesh transformMesh(const pbd_object::Mesh &mesh, 
                                   const Eigen::Matrix<Real, 1, 3> &translation,
                                   const Eigen::Quaternion<Real> &rotation,
                                   const Eigen::Matrix<Real, 1, 3> &scale);

    /*
    void readAttachedRobotForces();
    */

    // Creates the markers to publish
    void createRvizPointsMarker(const std::vector<Eigen::Matrix<Real,3,1>> *poses, 
                                visualization_msgs::Marker &marker);

    void createRvizLinesMarker(const std::vector<Eigen::Matrix<Real,3,1>> *poses,
                               const std::vector<Eigen::Quaternion<Real>> *orients,
                               const std::vector<Real> *lengths,
                               visualization_msgs::Marker &marker);

    void createRvizFramesMarker(const std::vector<Eigen::Matrix<Real,3,1>> *poses,
                               const std::vector<Eigen::Quaternion<Real>> *orients,
                               visualization_msgs::Marker &marker);

    void addAxisLineToMarker(visualization_msgs::Marker &marker, 
                            const Eigen::Matrix<Real,3,1> &start,
                            const Eigen::Matrix<Real,3,1> &end,
                            float r, float g, float b, float a); 

    void setupMinDistLineMarker(visualization_msgs::Marker &marker);

    // Timer callback functions
    void simulate(const ros::TimerEvent& e);
    void render(const ros::TimerEvent& e);
    
    void renderMarkers(const std::vector<Eigen::Matrix<Real,3,1>>* pos_ptr,
                        const std::vector<Eigen::Quaternion<Real>>* ori_ptr,
                        const std::vector<Real>* len_ptr,
                        ros::Publisher& publisher,
                        int marker_id);

    /*
    void publishWrenches(const ros::TimerEvent& e);
    */
    void renderRigidBodies(const ros::TimerEvent& e);

    void createMeshAndWireframeMarkers(const Eigen::Matrix<Real,Eigen::Dynamic,3> *vertices, 
                                                        const Eigen::MatrixX3i *face_tri_ids, 
                                                        visualization_msgs::Marker &meshMarker, 
                                                        visualization_msgs::Marker &wireframeMarker);
    void setupMeshMarker(visualization_msgs::Marker &marker);
    void setupWireframeMarker(visualization_msgs::Marker &marker);

    void publishMinDistancesToRigidBodies(const ros::TimerEvent& e);
    void publishDloState(const ros::TimerEvent& e);

    void publishMinDistLineMarkers(const std::vector<std::vector<utilities::CollisionHandler::MinDistanceData>>& min_distances_mt);

    // Odometry callback functions 
    
    void odometryCb_custom_static_particles(const nav_msgs::Odometry::ConstPtr& odom_msg, const int& id);
    
    void cmdVelCb_custom_static_particles(const geometry_msgs::Twist::ConstPtr& twist_msg, const int& id);
    
    // Change Dynamicity callback function
    void changeParticleDynamicityCb(const dlo_simulator_stiff_rods::ChangeParticleDynamicity::ConstPtr change_particle_dynamicity_msg);

    // Change Young Modulus and Torsion Modulus callback functions
    void changeYoungModulusCb(const std_msgs::Float32::ConstPtr young_modulus_msg);
    void changeTorsionModulusCb(const std_msgs::Float32::ConstPtr torsion_modulus_msg);
    
    // Service functions
    bool updateParams(std_srvs::Empty::Request& req, std_srvs::Empty::Response& res);
    void reset();
    void initialize() { std_srvs::Empty empt; updateParams(empt.request, empt.response); }
                                
    bool setParticleDynamicityCallback(dlo_simulator_stiff_rods::SetParticleDynamicity::Request &req,
                                       dlo_simulator_stiff_rods::SetParticleDynamicity::Response &res);

    // Create setYoungModulus service callback and setTorsionModulus service callback
    bool setDloYoungModulusCallback(dlo_simulator_stiff_rods::SetDloYoungModulus::Request &req,
                                    dlo_simulator_stiff_rods::SetDloYoungModulus::Response &res);
    bool setDloTorsionModulusCallback(dlo_simulator_stiff_rods::SetDloTorsionModulus::Request &req,
                                        dlo_simulator_stiff_rods::SetDloTorsionModulus::Response &res);

    // Create getYoungModulus service callback and getTorsionModulus service callback
    bool getDloYoungModulusCallback(dlo_simulator_stiff_rods::GetDloYoungModulus::Request &req,
                                    dlo_simulator_stiff_rods::GetDloYoungModulus::Response &res);
    bool getDloTorsionModulusCallback(dlo_simulator_stiff_rods::GetDloTorsionModulus::Request &req,
                                        dlo_simulator_stiff_rods::GetDloTorsionModulus::Response &res);

    bool enableCollisionHandlingCallback(dlo_simulator_stiff_rods::EnableCollisionHandling::Request &req,
                                        dlo_simulator_stiff_rods::EnableCollisionHandling::Response &res);

    // ROS variables---------------------------------
    ros::NodeHandle nh_;
    ros::NodeHandle nh_local_;
    boost::recursive_mutex &mtx_;

    ros::Publisher pub_dlo_state_;
    ros::Publisher pub_dlo_marker_array_;
    /*
    ros::Publisher pub_wrench_stamped_01_;
    ros::Publisher pub_wrench_stamped_02_;
    ros::Publisher pub_wrench_stamped_03_;
    ros::Publisher pub_wrench_stamped_04_;
    */

    ros::Publisher pub_rb_marker_array_;

    ros::Publisher pub_min_dist_to_rb_;
    ros::Publisher pub_min_dist_marker_array_;

    ros::ServiceServer params_srv_;

    ros::Subscriber sub_change_particle_dynamicity_;
    ros::ServiceServer set_particle_dynamicity_srv_;

    ros::Subscriber sub_change_young_modulus_;
    ros::Subscriber sub_change_torsion_modulus_;
    ros::ServiceServer set_young_modulus_srv_;
    ros::ServiceServer set_torsion_modulus_srv_;
    ros::ServiceServer get_young_modulus_srv_;
    ros::ServiceServer get_torsion_modulus_srv_;

    ros::ServiceServer enable_collision_handling_srv_;

    // Map to hold particle ID and its corresponding subscriber
    std::map<int, ros::Subscriber> custom_static_particles_odom_subscribers_;
    std::map<int, ros::Subscriber> custom_static_particles_cmd_vel_subscribers_;
    
    ros::Timer timer_render_;
    ros::Timer timer_simulate_;
    /*
    ros::Timer timer_wrench_pub_;
    */
    ros::Timer timer_render_rb_;
    ros::Timer timer_min_dist_to_rb_pub_; // renderMinDistancesToRigidBodies
    ros::Timer timer_dlo_state_pub_; // renderMinDistancesToRigidBodies

    // ROS Parameters
    bool p_active_;
    bool p_reset_;

    Real gravity_x_;
    Real gravity_y_;
    Real gravity_z_;
    
    Real dt_;
    bool set_sim_rate_auto_; //param to set the simulation rate and dt automatically

    int num_substeps_;
    int num_steps_;

    bool is_collision_handling_enabled_;
    bool visualize_min_distances_;

    Real contact_tolerance_;
    Real contact_sdf_domain_offset_;

    int dlo_visualization_mode_;
    int rb_visualization_mode_;

    // Dlo mesh properties: (Assuming dlo is a cylinder)
    Real dlo_l_;
    Real dlo_r_;
    Real dlo_density_;
    int dlo_num_segments_;

    Real dlo_zero_stretch_stiffness_;
    Real dlo_young_modulus_;
    Real dlo_torsion_modulus_;

    bool use_zero_stretch_stiffness_;

    Real global_damp_coeff_v_; 
    Real global_damp_coeff_w_;

    Real initial_height_; // initial dlo height from ground (m)
    std::vector<Real> dlo_translation_;
    std::vector<Real> dlo_rotationAxis_;
    Real dlo_rotationAngle_;
    std::vector<Real> dlo_scale_;

    int num_hang_corners_; // num of corners to hang dlo from (options: 0,1,2)
    std::vector<int> custom_static_particles_; // particle ids to set as static

    std::string custom_static_particles_odom_topic_prefix_;
    std::string custom_static_particles_cmd_vel_topic_prefix_;

    Real simulation_rate_;
    Real rendering_rate_;
    /*
    Real wrench_pub_rate_;
    */
    Real rendering_rb_rate_;
    Real min_dist_to_rb_pub_rate_;
    Real dlo_state_pub_rate_;

    std::string rb_scene_config_path_;

    std::string dlo_state_topic_name_;
    std::string dlo_markers_topic_name_;
    std::string dlo_frame_id_;
    
    /*
    std::string wrench_01_topic_name_;
    std::string wrench_02_topic_name_;
    std::string wrench_03_topic_name_;
    std::string wrench_04_topic_name_;

    std::string wrench_01_frame_id_;
    std::string wrench_02_frame_id_;
    std::string wrench_03_frame_id_;
    std::string wrench_04_frame_id_;
    */
    std::string min_dist_to_rb_topic_name_;
    std::string min_dist_markers_topic_name_;

    std::string change_particle_dynamicity_topic_name_;
    std::string set_particle_dynamicity_service_name_;

    std::string set_dlo_young_modulus_service_name_;
    std::string set_dlo_torsion_modulus_service_name_;
    std::string get_dlo_young_modulus_service_name_;
    std::string get_dlo_torsion_modulus_service_name_;
    std::string change_dlo_young_modulus_topic_name_;
    std::string change_dlo_torsion_modulus_topic_name_;

    std::string enable_collision_handling_service_name_;

    // Dlo visualization parameters 
    Real point_marker_scale_;

    std::vector<Real> point_marker_color_rgba_;

    Real line_marker_scale_multiplier_;

    std::vector<Real> line_marker_color_rgba_;

    Real frame_marker_scale_;
    Real frame_marker_axis_length_;

    // Rigid Body visualization parameters
    std::string rb_markers_topic_name_;

    Real rb_line_marker_scale_multiplier_;

    std::vector<Real> rb_line_marker_color_rgba_;

    std::vector<Real> rb_mesh_marker_color_rgba_;

    Real min_dist_line_marker_scale_multiplier_;

    std::vector<Real> min_dist_line_marker_color_rgba_;

    // Other variables
    Eigen::Matrix<Real,3,1> gravity_;

    bool is_auto_sim_rate_set_;
    
    bool is_rob_01_attached_;
    bool is_rob_02_attached_;
    bool is_rob_03_attached_;
    bool is_rob_04_attached_;

    int rob_01_attached_id_;
    int rob_02_attached_id_;
    int rob_03_attached_id_;
    int rob_04_attached_id_;

    /*
    Eigen::Matrix<Real,3,1> rob_01_attached_force_;
    Eigen::Matrix<Real,3,1> rob_02_attached_force_;
    Eigen::Matrix<Real,3,1> rob_03_attached_force_;
    Eigen::Matrix<Real,3,1> rob_04_attached_force_;
    */
    int time_frames_; //to report the performance per frame
    Real time_sum_; //to report the performance per frame

    int marker_id_; // for bookkeeping of the visualization markers of rigid bodies used in the code.

    // TODO: In the future, implement dlo and other deformable/rigid objects as PBDObjects
    // std::vector<PBDObject> sim_objects_;
    pbd_object::Dlo dlo_;
    std::vector<utilities::RigidBodySceneLoader::RigidBodyData> rigid_bodies_;

    // utilities::CollisionHandler collision_handler_;
    // std::unique_ptr<utilities::CollisionHandler> collision_handler_;
    utilities::CollisionHandler* collision_handler_;

    static void contactCallbackFunction(const unsigned int contactType,
                                        const unsigned int bodyIndex1, 
                                        const unsigned int bodyIndex2,
                                        const Eigen::Matrix<Real, 3, 1> &cp1, 
                                        const Eigen::Matrix<Real, 3, 1> &cp2,
                                        const Eigen::Matrix<Real, 3, 1> &normal, 
                                        const Real dist,
                                        const Real restitutionCoeff, 
                                        const Real frictionCoeffStatic,
                                        const Real frictionCoeffDynamic,
                                        void *userData);
    
}; // class DloSimulator
} // namespace dlo_simulator
#endif /* !DLO_SIMULATOR_H */