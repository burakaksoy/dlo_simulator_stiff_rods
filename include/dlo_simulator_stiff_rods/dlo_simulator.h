/*
 * Author: Burak Aksoy
 */

#ifndef DLO_SIMULATOR_H
#define DLO_SIMULATOR_H

#include <ros/ros.h>
#include <visualization_msgs/Marker.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/WrenchStamped.h>
#include <std_msgs/Float32.h>
#include <nav_msgs/Odometry.h>

#include <std_srvs/Empty.h>

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

#include "dlo_simulator_stiff_rods/utilities/dlo.h"

namespace dlo_simulator
{

class DloSimulator
{
public:
    DloSimulator(ros::NodeHandle &nh, ros::NodeHandle &nh_local, boost::recursive_mutex &mtx);
    ~DloSimulator();
private:
    // Functions ---------------------------------
    
    pbd_object::MeshDLO createMeshDLO(const std::string &name, const Real &dlo_l, const Real &dlo_z, const int &dlo_num_segments);

    /*
    void readAttachedRobotForces();
    */

    // Helper functions to publish created markers
    void publishRvizPoints(const std::vector<geometry_msgs::Point> &points);
    void publishRvizLines(const std::vector<geometry_msgs::Point> &points);

    // Creates the markers to publish
    void drawRviz(const std::vector<Eigen::Matrix<Real,3,1>> *poses);
    void drawRvizWireframe(const std::vector<Eigen::Matrix<Real,3,1>> *poses, const Eigen::Matrix2Xi *ids);
    void drawRvizRod(const std::vector<Eigen::Matrix<Real,3,1>> *poses,
                                         const std::vector<Eigen::Quaternion<Real>> *orients,
                                         const std::vector<Real> *lengths);

    // Timer callback functions
    void simulate(const ros::TimerEvent& e);
    void render(const ros::TimerEvent& e);
    /*
    void publishWrenches(const ros::TimerEvent& e);
    */

    // Odometry callback functions
    void odometryCb_01(const nav_msgs::Odometry::ConstPtr odom_msg);
    void odometryCb_02(const nav_msgs::Odometry::ConstPtr odom_msg);
    void odometryCb_03(const nav_msgs::Odometry::ConstPtr odom_msg);
    void odometryCb_04(const nav_msgs::Odometry::ConstPtr odom_msg);
    
    // Service functions
    bool updateParams(std_srvs::Empty::Request& req, std_srvs::Empty::Response& res);
    void reset();
    void initialize() { std_srvs::Empty empt; updateParams(empt.request, empt.response); }

    // ROS variables---------------------------------
    ros::NodeHandle nh_;
    ros::NodeHandle nh_local_;
    boost::recursive_mutex &mtx_;

    ros::Publisher pub_dlo_points_;
    /*
    ros::Publisher pub_wrench_stamped_01_;
    ros::Publisher pub_wrench_stamped_02_;
    ros::Publisher pub_wrench_stamped_03_;
    ros::Publisher pub_wrench_stamped_04_;
    */
    ros::ServiceServer params_srv_;
    
    ros::Subscriber sub_odom_01_;
    ros::Subscriber sub_odom_02_;
    ros::Subscriber sub_odom_03_;
    ros::Subscriber sub_odom_04_;
    
    ros::Timer timer_render_;
    ros::Timer timer_simulate_;
    /*
    ros::Timer timer_wrench_pub_;
    */
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

    int num_hang_corners_; // num of corners to hang dlo from (options: 0,1,2)

    bool use_direct_kkt_solver_;

    Real simulation_rate_;
    Real rendering_rate_;
    /*
    Real wrench_pub_rate_;
    */
    std::string dlo_points_topic_name_;
    std::string dlo_points_frame_id_;
    
    std::string odom_01_topic_name_;
    std::string odom_02_topic_name_;
    std::string odom_03_topic_name_;
    std::string odom_04_topic_name_;

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
    Real dlo_rob_z_offset_; // Additional  attachment height to robots

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

    // TODO: In the future, implement dlo and other deformable/rigid objects as PBDObjects
    // std::vector<PBDObject> sim_objects_;
    pbd_object::Dlo dlo_;
    
}; // class DloSimulator
} // namespace dlo_simulator
#endif /* !DLO_SIMULATOR_H */