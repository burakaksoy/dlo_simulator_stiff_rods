#ifndef ENERGY_BASED_MODEL_CERES_H
#define ENERGY_BASED_MODEL_CERES_H

#include <ros/ros.h> // For ROS logging
#include "dlo_simulator_stiff_rods/utilities/dlo.h"

#include "Eigen/Core"
#include "Eigen/Geometry"

#include <ceres/ceres.h>
#include <ceres/rotation.h>

namespace utilities
{
    class StretchBendTwistCost 
    {
    public:
        StretchBendTwistCost(Real avg_segment_length,
                            Eigen::Matrix<Real, 3, 1> connector0,
                            Eigen::Matrix<Real, 3, 1> connector1,
                            Eigen::Matrix<Real,3,1> rest_darboux_vector,
                            Eigen::Matrix<Real, 6, 6> sqrt_stiffness_mat):
                        avg_segment_length_(std::move(avg_segment_length)),
                        connector0_(std::move(connector0)),
                        connector1_(std::move(connector1)),
                        rest_darboux_vector_(std::move(rest_darboux_vector)),
                        sqrt_stiffness_mat_(std::move(sqrt_stiffness_mat)) {}

        template <typename T>
        bool operator()(const T* const p0_ptr,
                        const T* const q0_ptr,
                        const T* const p1_ptr,
                        const T* const q1_ptr,
                        T* residuals_ptr) const {
            Eigen::Map<const Eigen::Matrix<T, 3, 1>> p0(p0_ptr);
            Eigen::Map<const Eigen::Quaternion<T>> q0(q0_ptr);

            Eigen::Map<const Eigen::Matrix<T, 3, 1>> p1(p1_ptr);
            Eigen::Map<const Eigen::Quaternion<T>> q1(q1_ptr);


            // Represent the displacement between the two connectors in the Global frame.
            Eigen::Matrix<T, 3, 1> connector0 = q0 * connector0_.template cast<T>() + p0;
            Eigen::Matrix<T, 3, 1> connector1 = q1 * connector1_.template cast<T>() + p1;


            // Compute the relative transformation between the two frames.
            Eigen::Quaternion<T> q0_inv = q0.conjugate();
            // Current darboux vector (imaginary part of it)
            Eigen::Matrix<T,3,1> omega = T(2./avg_segment_length_) * (q0_inv * q1).vec(); //darboux vector


            // Compute the residuals.
            Eigen::Map<Eigen::Matrix<T, 6, 1>> residuals(residuals_ptr);
            
            // [ position         ]   [ delta_p          ]
            // [ orientation (3x1)] = [ 2 * delta_q(0:2) ]

            residuals.template block<3, 1>(0, 0) = connector0 - connector1;
            residuals.template block<3, 1>(3, 0) = omega - rest_darboux_vector_.template cast<T>();

            // Scale the residuals by the measurement uncertainty.
            residuals.applyOnTheLeft(sqrt_stiffness_mat_.template cast<T>());

            return true;
        }

        static ceres::CostFunction* Create(
            const Real& avg_segment_length,
            const Eigen::Matrix<Real, 3, 1>& connector0,
            const Eigen::Matrix<Real, 3, 1>& connector1,
            const Eigen::Matrix<Real,3,1>& rest_darboux_vector ,
            const Eigen::Matrix<Real, 6, 6>& sqrt_stiffness_mat) {
            return new ceres::AutoDiffCostFunction<StretchBendTwistCost, 6, 3, 4, 3, 4>(
                new StretchBendTwistCost(avg_segment_length, 
                                        connector0,
                                        connector1,
                                        rest_darboux_vector, 
                                        sqrt_stiffness_mat));
        }

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    private:
        // average segment length
        const Real avg_segment_length_;

        // Each segment's connector position to the next segment in the local frame of each segment.
        const Eigen::Matrix<Real, 3, 1> connector0_;
        const Eigen::Matrix<Real, 3, 1> connector1_;

        // rest darboux vector (imaginary part of it)
        const Eigen::Matrix<Real,3,1> rest_darboux_vector_;

        // The square root of the stiffness matrix. Diagonal matrix.
        const Eigen::Matrix<Real, 6, 6> sqrt_stiffness_mat_;
    };

    class GravityEnergyCost
    {
    public:
        GravityEnergyCost(Real mass, 
                        Eigen::Matrix<Real, 3, 1> gravity,
                        Eigen::Matrix<Real, 3, 1> ref_point): 
                            mass_(std::move(mass)),
                            gravity_(std::move(gravity)),
                            ref_point_(std::move(ref_point)) {}

        template <typename T>
        bool operator()(const T *const p0_ptr, 
                        T *residual) const
        {
            Eigen::Map<const Eigen::Matrix<T, 3, 1>> p0(p0_ptr);

            // Compute the height vector
            Eigen::Matrix<T, 3, 1> height = p0 - ref_point_.template cast<T>();

            // Compute the dot product of gravity and height
            T gravity_dot_height = -gravity_.template cast<T>().dot(height);

            // Ensure the dot product is non-negative to avoid sqrt of negative number
            T safe_gravity_dot_height = ceres::fmax(T(1e-15), gravity_dot_height);

            // Compute the residual
            residual[0] = ceres::sqrt(T(2.0 * mass_) * safe_gravity_dot_height);

            return true;
        }

        static ceres::CostFunction* Create(
            const Real& mass, 
            const Eigen::Matrix<Real, 3, 1>& gravity,
            const Eigen::Matrix<Real, 3, 1>& ref_point) {
            return new ceres::AutoDiffCostFunction<GravityEnergyCost, 1, 3>(
                new GravityEnergyCost(mass, gravity, ref_point));
        }

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    private:
        const Real mass_;
        const Eigen::Matrix<Real, 3, 1> gravity_;
        const Eigen::Matrix<Real, 3, 1> ref_point_;
    };

    /* TODO
    class GravityForceCost 
    {
    public:
        GravityForceCost(Real avg_segment_length0,
                            Real avg_segment_length1,
                            Eigen::Matrix<Real, 3, 1> connector00,
                            Eigen::Matrix<Real, 3, 1> connector01,
                            Eigen::Matrix<Real, 3, 1> connector10,
                            Eigen::Matrix<Real, 3, 1> connector11,
                            Eigen::Matrix<Real,3,1> rest_darboux_vector,
                            Eigen::Matrix<Real, 3, 3> sqrt_stiffness_mat):
                        avg_segment_length0_(std::move(avg_segment_length0)),
                        avg_segment_length1_(std::move(avg_segment_length1)),
                        connector00_(std::move(connector00)),
                        connector01_(std::move(connector01)),
                        connector10_(std::move(connector10)),
                        connector11_(std::move(connector11)),
                        rest_darboux_vector_(std::move(rest_darboux_vector)),
                        sqrt_stiffness_mat_(std::move(sqrt_stiffness_mat)) {}

        template <typename T>
        bool operator()(const T* const p1_ptr,
                        const T* const q1_ptr,
                        T* residuals_ptr) const {

            Eigen::Map<const Eigen::Matrix<T, 3, 1>> p1(p1_ptr);
            Eigen::Map<const Eigen::Quaternion<T>> q1(q1_ptr);


            // Represent the displacement between the two connectors in the Global frame.
            Eigen::Matrix<T, 3, 1> connector0 = q0 * connector0_.template cast<T>() + p0;
            Eigen::Matrix<T, 3, 1> connector1 = q1 * connector1_.template cast<T>() + p1;


            // Compute the relative transformation between the two frames.
            Eigen::Quaternion<T> q0_inv = q0.conjugate();
            // Current darboux vector (imaginary part of it)
            Eigen::Matrix<T,3,1> omega = T(2./avg_segment_length_) * (q0_inv * q1).vec(); //darboux vector


            // Compute the residuals as net force acting on the segment
            Eigen::Map<Eigen::Matrix<T, 3, 1>> residuals(residuals_ptr);
            
            residuals.template block<3, 1>(0, 0) = connector0 - connector1;

            // Scale the residuals by the measurement uncertainty.
            residuals.applyOnTheLeft(sqrt_stiffness_mat_.template cast<T>());

            return true;
        }

        static ceres::CostFunction* Create(
            const Real& avg_segment_length0,
            const Real& avg_segment_length1,
            const Eigen::Matrix<Real, 3, 1>& connector00,
            const Eigen::Matrix<Real, 3, 1>& connector01,
            const Eigen::Matrix<Real, 3, 1>& connector10,
            const Eigen::Matrix<Real, 3, 1>& connector11,
            const Eigen::Matrix<Real,3,1>& rest_darboux_vector ,
            const Eigen::Matrix<Real, 3, 3>& sqrt_stiffness_mat) {
            return new ceres::AutoDiffCostFunction<GravityForceCost, 3, 3, 4>(
                new GravityForceCost(avg_segment_length, 
                                        connector0,
                                        connector1,
                                        rest_darboux_vector, 
                                        sqrt_stiffness_mat));
        }

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    private:
        // average segment length
        const Real avg_segment_length0_;
        const Real avg_segment_length1_;

        // Each segment's connector position to the next segment in the local frame of each segment.
        const Eigen::Matrix<Real, 3, 1> connector00_;
        const Eigen::Matrix<Real, 3, 1> connector01_;
        const Eigen::Matrix<Real, 3, 1> connector10_;
        const Eigen::Matrix<Real, 3, 1> connector11_;

        // rest darboux vector (imaginary part of it)
        const Eigen::Matrix<Real,3,1> rest_darboux_vector_;

        // The square root of the stiffness matrix. Diagonal matrix.
        const Eigen::Matrix<Real, 3, 3> sqrt_stiffness_mat_;
    };
    */

} // namespace utilities

#endif // ENERGY_BASED_MODEL_CERES_H