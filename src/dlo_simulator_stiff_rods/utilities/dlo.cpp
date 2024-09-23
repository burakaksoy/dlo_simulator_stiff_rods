/*
 * Author: Burak Aksoy
 */

#include "dlo_simulator_stiff_rods/utilities/dlo.h"

using namespace pbd_object;

const Real eps_ = static_cast<Real>(1e-10);

Dlo::Dlo(){

}

Dlo::Dlo(const MeshDLO &mesh, 
        const Real &zero_stretch_stiffness, 
        const Real &young_modulus, 
        const Real &torsion_modulus, 
        const Real &density,
        const Real &radius,
        const bool &use_zero_stretch_stiffness,
        const Real &global_damp_coeff_v,
        const Real &global_damp_coeff_w):
    mesh_(mesh),
    zero_stretch_stiffness_(zero_stretch_stiffness),
    young_modulus_(young_modulus),
    torsion_modulus_(torsion_modulus),
    density_(density),
    radius_(radius),
    use_zero_stretch_stiffness_(use_zero_stretch_stiffness),
    global_damp_coeff_v_(global_damp_coeff_v),
    global_damp_coeff_w_(global_damp_coeff_w)
{
    num_particles_ = mesh_.vertices.size();
    num_quaternions_ = mesh_.quaternions.size();

    std::cout << "num particles: " << num_particles_ << std::endl;
    std::cout << "num quaternions: " << num_quaternions_ << std::endl;
    
    // Particle data
    pos_ = mesh_.vertices;
    prev_pos_ = mesh_.vertices;
    rest_pos_ = mesh_.vertices;
    vel_ = Eigen::Matrix<Real,3,Eigen::Dynamic>::Zero(3,num_particles_); // Create a velocity vector of particles filled with zeros 
    for_ = Eigen::Matrix<Real,3,Eigen::Dynamic>::Zero(3,num_particles_);
    inv_mass_.assign(num_particles_, 0.0);
    is_dynamic_.assign(num_particles_, true); // initially assume all particles are dynamic

    // std::cout << "pos_:\n" << pos_ << std::endl;
    // std::cout << "prev_pos_:\n" << prev_pos_ << std::endl;
    // std::cout << "vel_:\n" << vel_ << std::endl;
    // std::cout << "inv_mass_:\n" << inv_mass_ << std::endl;
    
    // Orientation Data
    ori_ = mesh_.quaternions;
    prev_ori_ = mesh_.quaternions;
    rest_ori_ = mesh_.quaternions;
    omega_ = Eigen::Matrix<Real,3,Eigen::Dynamic>::Zero(3,num_quaternions_); // Create a angular velocity vector of orientations filled with zeros 
    tor_ = Eigen::Matrix<Real,3,Eigen::Dynamic>::Zero(3,num_quaternions_);
    inv_iner_.assign(num_quaternions_,Eigen::Matrix<Real,3,3>::Zero());
    iner_.assign(num_quaternions_,Eigen::Matrix<Real,3,3>::Zero());

    // std::cout << "ori_:\n" << ori_ << std::endl;
    // std::cout << "prev_ori_:\n" << prev_ori_ << std::endl;
    // std::cout << "omega_:\n" << omega_ << std::endl;
    // std::cout << "inv_iner_:\n" << inv_iner_ << std::endl;
    // std::cout << "iner_:\n" << iner_ << std::endl;

    // Initialize constraint values for each constraint
    RHS_.assign(num_quaternions_-1,Eigen::Matrix<Real,6,1>::Zero());

    // Initialize 2 jacobians for each constraint
    bendingAndTorsionJacobians_.resize(num_quaternions_-1);
	std::vector<Eigen::Matrix<Real,3,3>> sampleJacobians(2);
	sampleJacobians[0].setZero();
	sampleJacobians[1].setZero();
	std::fill(bendingAndTorsionJacobians_.begin(), bendingAndTorsionJacobians_.end(), sampleJacobians);

    setMasses();

    setStretchBendTwistConstraints();

    // attached_ids_  //Initially empty vector of integers to store the ids of attached (fixed) particles.
}

Dlo::~Dlo(){

}


void Dlo::setStretchBendTwistConstraints(){
    /* 
    Specifies:
    - stretchBendTwist_ids_ (from DLO mesh)
    - stretchBendTwist_restDarbouxVectors_
    - stretchBendTwist_constraintPosInfo_:
        col 0:	connector in segment 0 (local)
        col 1:	connector in segment 1 (local)
        col 2:	connector in segment 0 (global)
        col 3:	connector in segment 1 (global)    
    - average_segment_lengths_
    */

    std::vector<Eigen::Vector2i> quaternion_ids;
    for(unsigned int id = 0; id < num_quaternions_-1; id++){
        Eigen::Vector2i ids(id, id+1); // quaternion ids
        quaternion_ids.push_back(ids);
    }

    Eigen::MatrixXi quaternion_ids_mat(2,quaternion_ids.size());
    for (int i = 0; i < quaternion_ids.size(); i++) {
        quaternion_ids_mat.col(i) = quaternion_ids[i];
    }

    stretchBendTwist_ids_ = quaternion_ids_mat;

    for (int i = 0; i < stretchBendTwist_ids_.cols(); i++){
        int id0 = stretchBendTwist_ids_(0,i);
        int id1 = stretchBendTwist_ids_(1,i);
        
        // compute rest Darboux vector based on eqn 7    
        Eigen::Quaternion<Real> rest_darboux_vect = ori_[id0].conjugate() * ori_[id1]; 
        Real averageSegmentLength = 0.5*(mesh_.segment_lengths[id0] + mesh_.segment_lengths[id1]);
        stretchBendTwist_restDarbouxVectors_.push_back((2./averageSegmentLength)*rest_darboux_vect.vec());
        average_segment_lengths_.push_back(averageSegmentLength);

        // set initial constraint position info 
        const Eigen::Matrix<Real,3,3> rot0 = ori_[id0].toRotationMatrix();
        const Eigen::Matrix<Real,3,3> rot1 = ori_[id1].toRotationMatrix();

        Eigen::Matrix<Real, 3, 4> constraintPosInfo;
        // locally each constraint position is along the z axis at the half lenght of the segment
        constraintPosInfo.col(0) = Eigen::Matrix<Real,3,1>(0,0,0.5*mesh_.segment_lengths[id0]);
        constraintPosInfo.col(1) = Eigen::Matrix<Real,3,1>(0,0,-0.5*mesh_.segment_lengths[id1]);
        constraintPosInfo.col(2) = rot0 * constraintPosInfo.col(0) + mesh_.vertices[id0];
        constraintPosInfo.col(3) = rot1 * constraintPosInfo.col(1) + mesh_.vertices[id1];

        stretchBendTwist_constraintPosInfo_.push_back(constraintPosInfo);
    }
}


void Dlo::setMasses(){
    /* 
    Specifies:
    - inv_mass_, using the edge rest lengths, DLO length, radius and density information.
    - inv_iner_,
    - iner_
    */

    for (int i = 0; i < num_particles_; i++) {
        const Real &l = mesh_.segment_lengths[i]; // segment length
        Real V = M_PI*(radius_*radius_)*l;  // segment volume
        Real mass = V * density_; // segment mass

        if (mass > 0.0) {
            Real p_0_mass = ( inv_mass_[i] > 0.0 ) ? 1.0/inv_mass_[i] : 0.0;

            Real q_0_inertia_xx = ( inv_iner_[i](0,0) > 0.0 ) ? 1.0/inv_iner_[i](0,0) : 0.0;
            Real q_0_inertia_yy = ( inv_iner_[i](1,1) > 0.0 ) ? 1.0/inv_iner_[i](1,1) : 0.0;
            Real q_0_inertia_zz = ( inv_iner_[i](2,2) > 0.0 ) ? 1.0/inv_iner_[i](2,2) : 0.0;

            p_0_mass += mass;
            
            q_0_inertia_xx += (0.25)*mass*(radius_*radius_) + (1./12.)*mass*(l*l);
            q_0_inertia_yy += (0.25)*mass*(radius_*radius_) + (1./12.)*mass*(l*l);
            q_0_inertia_zz += 0.5*mass*(radius_*radius_);

            inv_mass_[i] = 1.0/p_0_mass;
            
            inv_iner_[i](0,0) = 1.0/q_0_inertia_xx;
            inv_iner_[i](1,1) = 1.0/q_0_inertia_yy;
            inv_iner_[i](2,2) = 1.0/q_0_inertia_zz;

            iner_[i](0,0) = q_0_inertia_xx;
            iner_[i](1,1) = q_0_inertia_yy;
            iner_[i](2,2) = q_0_inertia_zz;
        }
    }

    // To debug
    Eigen::Matrix<Real,1,Eigen::Dynamic> inv_mass_eigen(num_particles_);
    for (int i = 0; i < inv_mass_.size(); i++) {
        inv_mass_eigen(i) = inv_mass_[i];
    }
    // std::cout << "inv_mass_:\n" << inv_mass_eigen << std::endl;
    std::cout << "particle masses:\n" << inv_mass_eigen.cwiseInverse() << " kg." << std::endl;
    std::cout << "Total dlo mass:\n" << inv_mass_eigen.cwiseInverse().sum() << " kg." << std::endl;
}


void Dlo::changeParticleDynamicity(const int &particle, const bool &is_dynamic){
    if(particle < 0 || particle >= is_dynamic_.size()) {
        throw std::out_of_range("Index out of bounds");
    }
    if (is_dynamic_[particle] != is_dynamic){
        is_dynamic_[particle] = is_dynamic;

        if (!is_dynamic){
            // Update Velocity to Zero when making the particle static
            updateAttachedVelocity(particle, 
                                    Eigen::Matrix<Real,3,1>::Zero(), 
                                    Eigen::Matrix<Real,3,1>::Zero());
        }
    }
}

void Dlo::setStaticParticles(const std::vector<int> &particles){
    for (const int& i : particles)
    {
        changeParticleDynamicity(i, false);
    }
}

void Dlo::setDynamicParticles(const std::vector<int> &particles){
    for (const int& i : particles)
    {
        changeParticleDynamicity(i, true);
    }
}

const bool Dlo::isStaticParticle(const int &particle){
    if(particle < 0 || particle >= is_dynamic_.size()) {
        throw std::out_of_range("Index out of bounds");
    }
    return !is_dynamic_[particle];
}

const bool Dlo::isDynamicParticle(const int &particle){
    if(particle < 0 || particle >= is_dynamic_.size()) {
        throw std::out_of_range("Index out of bounds");
    }
    return is_dynamic_[particle];
}

void Dlo::setYoungModulus(const Real &young_modulus){
    young_modulus_ = young_modulus;
}

void Dlo::setTorsionModulus(const Real &torsion_modulus){
    torsion_modulus_ = torsion_modulus;
}

const Real Dlo::getYoungModulus(){
    return young_modulus_;
}

const Real Dlo::getTorsionModulus(){
    return torsion_modulus_;
}

void Dlo::preSolve(const Real &dt, const Eigen::Matrix<Real,3,1> &gravity){
    #pragma omp parallel default(shared)
    {
        // Semi implicit euler (position)
        #pragma omp for schedule(static) 
        // #pragma omp parallel for 
        for (int i = 0; i< num_particles_; i++){
            if (is_dynamic_[i]){
                vel_.col(i) += gravity*dt;
            }
            prev_pos_[i] = pos_[i];
            pos_[i] += vel_.col(i)*dt;
        }

        // Semi implicit euler (rotation)
        #pragma omp for schedule(static) 
        // #pragma omp parallel for
        for (int i = 0; i< num_quaternions_; i++){
            if (is_dynamic_[i]){
                //assume zero external torque.
                Eigen::Matrix<Real,3,1> torque = Eigen::Matrix<Real,3,1>::Zero();    
                // integration 
                omega_.col(i) += dt * inv_iner_[i] * (torque - (omega_.col(i).cross(iner_[i]*omega_.col(i))));
            }
            prev_ori_[i] = ori_[i];
            Eigen::Quaternion<Real> angVelQ(0.0, omega_.col(i)(0), omega_.col(i)(1), omega_.col(i)(2));
            ori_[i].coeffs() += dt * 0.5 * (angVelQ * ori_[i]).coeffs();
            ori_[i].normalize();
        }
    }
}

void Dlo::solve(const Real &dt){
    solveStretchBendTwistConstraints(dt);
}

void Dlo::solveStretchBendTwistConstraints(const Real &dt){
    // Inits before projection (alpha tilde, eqn 24)
    Real dt_sqr = (dt*dt);
    Real inv_dt_sqr = static_cast<Real>(1.0)/dt_sqr;
    // compute compliance parameter of the stretch constraint part
    Eigen::Matrix<Real,3,1> stretch_compliance; // upper diagonal of eqn 24

    if (use_zero_stretch_stiffness_){
        stretch_compliance <<
            inv_dt_sqr / (zero_stretch_stiffness_),
            inv_dt_sqr / (zero_stretch_stiffness_),
            inv_dt_sqr / (zero_stretch_stiffness_);
    }
    else{
        // Calculate Stretch compliance based on E,G and dlo_r
        Real area(static_cast<Real>(M_PI) * std::pow(radius_, static_cast<Real>(2.0)));
        stretch_compliance <<
            inv_dt_sqr / (zero_stretch_stiffness_*torsion_modulus_*area),
            inv_dt_sqr / (zero_stretch_stiffness_*torsion_modulus_*area),
            inv_dt_sqr / (zero_stretch_stiffness_*young_modulus_*area);
    }

    // compute compliance parameter of the bending and torsion constraint part
    Eigen::Matrix<Real,3,1> bending_and_torsion_compliance_fixed; //lower diagonal of eqn 24
    Eigen::Matrix<Real,3,1> bending_and_torsion_compliance; // 

    Real secondMomentOfArea(static_cast<Real>(M_PI_4) * std::pow(radius_, static_cast<Real>(4.0)));
	Real bendingStiffness(young_modulus_ * secondMomentOfArea);
	Real torsionStiffness(static_cast<Real>(2.0) * torsion_modulus_ * secondMomentOfArea);

    // assumption: the rod axis follows the z-axis of the local frame
    bending_and_torsion_compliance_fixed << 
        inv_dt_sqr / bendingStiffness,
        inv_dt_sqr / bendingStiffness,
        inv_dt_sqr / torsionStiffness;

    // // Damping coeffs.
    // Real damping_shear_   = 0.0; // 5.0e-5; // 1.0e-9;
    // Real damping_stretch_ = 0.0; // 5.0e-5; // 1.0e-9;
    // Real damping_bending_ = 5.0e-4;
    // Real damping_torsion_ = 0.0; //5.0e-5;
    
    // Eigen::Matrix<Real,3,1> beta_tilde_shear_stretch;

    // beta_tilde_shear_stretch << 
    //     dt_sqr * damping_shear_,
    //     dt_sqr * damping_shear_,
    //     dt_sqr * damping_stretch_;
    
    // Eigen::Matrix<Real,3,1> beta_tilde_bending_torsion;

    // beta_tilde_bending_torsion <<
    //     dt_sqr * damping_bending_,
    //     dt_sqr * damping_bending_,
    //     dt_sqr * damping_torsion_;


    max_error_ = 0.0;

    // updates on the constraints
    // For each constraint
    int i = 0;
    const int n = stretchBendTwist_restDarbouxVectors_.size(); // num. of constraints
    for (int itr = 0; itr < n; itr++){
        // Bilateral interleaving order
        if (n % 2 == 0) {
            (itr % 2 == 0) ? i = itr :  i = n-itr;
        } else {
            (itr % 2 == 0) ? i = itr :  i = n-itr-1;
        }

        // IDs 
        const int& id0 = stretchBendTwist_ids_(0,i);
        const int& id1 = stretchBendTwist_ids_(1,i);

        const Real& averageSegmentLength = average_segment_lengths_[i];

        bending_and_torsion_compliance = bending_and_torsion_compliance_fixed * (static_cast<Real>(1.0)/averageSegmentLength); // why?
        // Because to remove the dependency of the stiffness to the segment length.

        // TODO: Do the same thing above for stretch compliance when not use_zero_stretch_stiffness_?

        // inverse masses of these segments
        const Real& invMass0 = inv_mass_[id0];
        const Real& invMass1 = inv_mass_[id1];

        const bool &isDynamic0 = is_dynamic_[id0];
        const bool &isDynamic1 = is_dynamic_[id1];

        // Current positions
        Eigen::Matrix<Real,3,1>& p0 = pos_[id0];
        Eigen::Matrix<Real,3,1>& p1 = pos_[id1];

        // Current orientations
        Eigen::Quaternion<Real>& q0 = ori_[id0];
        Eigen::Quaternion<Real>& q1 = ori_[id1];
        const Eigen::Matrix<Real,3,3> rot0 = q0.toRotationMatrix();
        const Eigen::Matrix<Real,3,3> rot1 = q1.toRotationMatrix();

        // // Previous positions for damping
        // const Eigen::Matrix<Real,3,1>& p0_prev = prev_pos_[id0];
        // const Eigen::Matrix<Real,3,1>& p1_prev = prev_pos_[id1];

        // // Previous orientations for damping
        // const Eigen::Quaternion<Real>& q0_prev = prev_ori_[id0];
        // const Eigen::Quaternion<Real>& q1_prev = prev_ori_[id1];

        // inverse inertia matrices of these segments (TODO: these needs to be rotated to world frame)
        Eigen::Matrix<Real,3,3> inertiaInverseW0 = inv_iner_[id0]; // local
        Eigen::Matrix<Real,3,3> inertiaInverseW1 = inv_iner_[id1]; // local

        inertiaInverseW0 = rot0*inertiaInverseW0 * rot0.transpose(); //world
        inertiaInverseW1 = rot1*inertiaInverseW1 * rot1.transpose(); //world

        // Current constraint pos info needs to be updated
        Eigen::Matrix<Real, 3, 4>& constraintPosInfo = stretchBendTwist_constraintPosInfo_[i];

        // update constraint (for eqn 23, upper part)
        constraintPosInfo.col(2) = rot0 * constraintPosInfo.col(0) + p0;
        constraintPosInfo.col(3) = rot1 * constraintPosInfo.col(1) + p1;

        const Eigen::Matrix<Real,3,1>& connector0 = constraintPosInfo.col(2);
        const Eigen::Matrix<Real,3,1>& connector1 = constraintPosInfo.col(3);

        // rest darboux vector (imaginary part of it)
        Eigen::Matrix<Real,3,1>& restDarbouxVector = stretchBendTwist_restDarbouxVectors_[i]; 

        // Current darboux vector (imaginary part of it)
        Eigen::Matrix<Real,3,1> omega = (2./averageSegmentLength)*(q0.conjugate() * q1).vec();   //darboux vector

        // Compute zero-stretch part of constraint violation (eqn 23, upper part)
        // Compute bending and torsion part of constraint violation (eqn 23, lower part)

        // fill right hand side of the linear equation system (Equation (19))
        Eigen::Matrix<Real, 6, 1>& rhs= RHS_[i];
        rhs.block<3, 1>(0, 0) = - (connector0 - connector1); // stretchViolation;
        rhs.block<3, 1>(3, 0) = - (omega - restDarbouxVector); // bendingAndTorsionViolation;

        // compute max error
		for (unsigned char j(0); j < 6; ++j)
		{
			max_error_ = std::max(max_error_, std::abs(rhs[j]));
		}
        
        // compute G matrices (Equation (27)), but first the bottom part of the matrix
        Eigen::Matrix<Real, 4, 3> G0, G1;
        // w component at index 3
        G0 <<
             q0.w(),  q0.z(), -q0.y(),
            -q0.z(),  q0.w(),  q0.x(),
             q0.y(), -q0.x(),  q0.w(),
            -q0.x(), -q0.y(), -q0.z();
        G0 *= static_cast<Real>(0.5);
        // w component at index 3
        G1 <<
             q1.w(),  q1.z(), -q1.y(),
            -q1.z(),  q1.w(),  q1.x(),
             q1.y(), -q1.x(),  q1.w(),
            -q1.x(), -q1.y(), -q1.z();
        G1 *= static_cast<Real>(0.5);

        // compute bending and torsion Jacobians (Equation (10) and Equation (11)), but first the right part of the matrix
	    Eigen::Matrix<Real, 3, 4> jOmega0, jOmega1;
        // w component at index 3, Equation (11)
        jOmega0 <<
            -q1.w(), -q1.z(),  q1.y(), q1.x(),
             q1.z(), -q1.w(), -q1.x(), q1.y(),
            -q1.y(),  q1.x(), -q1.w(), q1.z();
        // w component at index 3, Equation (10)
        jOmega1 <<
             q0.w(),  q0.z(), -q0.y(), -q0.x(),
            -q0.z(),  q0.w(),  q0.x(), -q0.y(),
             q0.y(), -q0.x(),  q0.w(), -q0.z();
        jOmega0 *= static_cast<Real>(2.0) / averageSegmentLength;
        jOmega1 *= static_cast<Real>(2.0) / averageSegmentLength;

        // Lower right part of eqn 25
        Eigen::Matrix<Real, 3, 3> & jOmegaG0 = bendingAndTorsionJacobians_[i][0];
        jOmegaG0 = jOmega0*G0;
        

        // Lower right part of eqn 26
        Eigen::Matrix<Real, 3, 3> & jOmegaG1 = bendingAndTorsionJacobians_[i][1];
        jOmegaG1 = jOmega1*G1;

        // Start actual solving from here -------
        // compute matrix of the linear equation system (using Equations (25), (26), and (28) in Equation (19))
        Eigen::Matrix<Real, 6, 6> JMJT; // = Eigen::Matrix<Real, 6, 6>::Zero(); // Initialize place holder for J*M^-1*J^T

        // compute stretch block
        Eigen::Matrix<Real,3,3> K1, K2;
        if (isDynamic0){ computeMatrixK(connector0, invMass0, p0, inertiaInverseW0, K1);}
        else{K1.setZero();}
        if (isDynamic1){ computeMatrixK(connector1, invMass1, p1, inertiaInverseW1, K2);}
        else{K2.setZero();}
        JMJT.block<3, 3>(0, 0) = K1 + K2;

        // compute coupling blocks
        const Eigen::Matrix<Real,3,1> ra = connector0 - p0;
        const Eigen::Matrix<Real,3,1> rb = connector1 - p1;

        Eigen::Matrix<Real,3,3> ra_crossT, rb_crossT;
        crossProductMatrix(-ra, ra_crossT); // use -ra to get the transpose
        crossProductMatrix(-rb, rb_crossT); // use -rb to get the transpose

        Eigen::Matrix<Real,3,3> offdiag(Eigen::Matrix<Real,3,3>::Zero());
        if (isDynamic0){offdiag = jOmegaG0 * inertiaInverseW0 * ra_crossT * (-1);}
        if (isDynamic1){offdiag += jOmegaG1 * inertiaInverseW1 * rb_crossT;}
        JMJT.block<3, 3>(3, 0) = offdiag;
        JMJT.block<3, 3>(0, 3) = offdiag.transpose();

        // compute bending and torsion block
        Eigen::Matrix<Real,3,3> MInvJT0 = inertiaInverseW0 * jOmegaG0.transpose();
        Eigen::Matrix<Real,3,3> MInvJT1 = inertiaInverseW1 * jOmegaG1.transpose();

        Eigen::Matrix<Real,3,3> JMJTOmega(Eigen::Matrix<Real,3,3>::Zero());
        if (isDynamic0){ JMJTOmega = jOmegaG0*MInvJT0;}
        if (isDynamic1){ JMJTOmega += jOmegaG1*MInvJT1;}
        JMJT.block<3, 3>(3, 3) = JMJTOmega;

        // // Multiply each row with compliance times damping divided by dt
        // JMJT.row(0) *= (1.0 + (stretch_compliance(0)*beta_tilde_shear_stretch(0)/dt)); 
        // JMJT.row(1) *= (1.0 + (stretch_compliance(1)*beta_tilde_shear_stretch(1)/dt)); 
        // JMJT.row(2) *= (1.0 + (stretch_compliance(2)*beta_tilde_shear_stretch(2)/dt)); 
        // JMJT.row(3) *= (1.0 + (bending_and_torsion_compliance(0)*beta_tilde_bending_torsion(0)/dt)); 
        // JMJT.row(4) *= (1.0 + (bending_and_torsion_compliance(1)*beta_tilde_bending_torsion(1)/dt)); 
        // JMJT.row(5) *= (1.0 + (bending_and_torsion_compliance(2)*beta_tilde_bending_torsion(2)/dt)); 

        // add compliance
        JMJT(0, 0) += stretch_compliance(0);
        JMJT(1, 1) += stretch_compliance(1);
        JMJT(2, 2) += stretch_compliance(2);
        JMJT(3, 3) += bending_and_torsion_compliance(0);
        JMJT(4, 4) += bending_and_torsion_compliance(1);
        JMJT(5, 5) += bending_and_torsion_compliance(2);

        // // Calculate Jv = J1*v1 + J2*v2 for damping
        // Eigen::Matrix<Real,6,1> Jv;

        // // Current velocities
        // const Eigen::Matrix<Real,3,1> v0 = (p0 - p0_prev)/dt;
        // const Eigen::Matrix<Real,3,1> v1 = (p1 - p1_prev)/dt;

        // // Current angular velocities
        // const Eigen::Matrix<Real,3,1> w0 = 2.0*(q0 * q0_prev.conjugate()).vec()/dt;
        // const Eigen::Matrix<Real,3,1> w1 = 2.0*(q1 * q1_prev.conjugate()).vec()/dt;

        // Jv.block<3,1>(0,0) = v0 - v1 + ra_crossT*w0 - rb_crossT*w1;
        // Jv.block<3,1>(3,0) = jOmegaG0 * w0 + jOmegaG1 * w1;

        // //Update rhs with -alpha_tilde*beta_tilde*Jv
        // rhs(0) -= stretch_compliance(0)*beta_tilde_shear_stretch(0) * Jv(0); 
        // rhs(1) -= stretch_compliance(1)*beta_tilde_shear_stretch(1) * Jv(1); 
        // rhs(2) -= stretch_compliance(2)*beta_tilde_shear_stretch(2) * Jv(2); 
        // rhs(3) -= bending_and_torsion_compliance(0)*beta_tilde_bending_torsion(0) * Jv(3); 
        // rhs(4) -= bending_and_torsion_compliance(1)*beta_tilde_bending_torsion(1) * Jv(4); 
        // rhs(5) -= bending_and_torsion_compliance(2)*beta_tilde_bending_torsion(2) * Jv(5); 

        // solve linear equation system (Equation 19)
        Eigen::Matrix<Real,6,1> deltaLambda = JMJT.ldlt().solve(rhs);

        // compute position and orientation updates (using Equations (25), (26), and (28) in Equation (20))
        const Eigen::Matrix<Real,3,1> & deltaLambdaStretch = deltaLambda.block<3, 1>(0, 0);
        const Eigen::Matrix<Real,3,1> & deltaLambdaBendingAndTorsion = deltaLambda.block<3, 1>(3, 0);

        // Now apply the corrections -----------------------------
        if (isDynamic0)
        {
            p0 += invMass0 * deltaLambdaStretch;
            q0.coeffs() += G0 * (inertiaInverseW0 * -ra_crossT * deltaLambdaStretch + MInvJT0 * deltaLambdaBendingAndTorsion);
            q0.normalize();
        }

        if (isDynamic1)
        {
            p1 += -invMass1 * deltaLambdaStretch;
            q1.coeffs() += G1 * (inertiaInverseW1 * rb_crossT * deltaLambdaStretch + MInvJT1 * deltaLambdaBendingAndTorsion);
            q1.normalize();
        }

        // Calculate the Force/Torque (Wrench) at each segment ---------------
        for_.col(id0) += inv_dt_sqr * deltaLambdaStretch;
        for_.col(id1) += inv_dt_sqr * -deltaLambdaStretch;

        tor_.col(id0) += inv_dt_sqr * (-ra_crossT * deltaLambdaStretch + jOmegaG0.transpose() * deltaLambdaBendingAndTorsion);
        tor_.col(id1) += inv_dt_sqr * ( rb_crossT * deltaLambdaStretch + jOmegaG1.transpose() * deltaLambdaBendingAndTorsion);
    }
}

void Dlo::computeMatrixK(const Eigen::Matrix<Real,3,1> &connector, 
                        const Real invMass, 
                        const Eigen::Matrix<Real,3,1> &x, 
                        const Eigen::Matrix<Real,3,3> &inertiaInverseW, 
                        Eigen::Matrix<Real,3,3> &K) {
    // This function computes the upper left block of J*M^-1*J.T directly.
	if (invMass != 0.0)
	{
        // vector from center of mass to conneting point in world frame
		const Eigen::Matrix<Real,3,1> v = connector - x;
		const Real a = v[0];
		const Real b = v[1];
		const Real c = v[2];

		// J is symmetric
		const Real j11 = inertiaInverseW(0, 0);
		const Real j12 = inertiaInverseW(0, 1);
		const Real j13 = inertiaInverseW(0, 2);
		const Real j22 = inertiaInverseW(1, 1);
		const Real j23 = inertiaInverseW(1, 2);
		const Real j33 = inertiaInverseW(2, 2);

		K(0, 0) = c*c*j22 - b*c*(j23 + j23) + b*b*j33 + invMass;
		K(0, 1) = -(c*c*j12) + a*c*j23 + b*c*j13 - a*b*j33;
		K(0, 2) = b*c*j12 - a*c*j22 - b*b*j13 + a*b*j23;
		K(1, 0) = K(0, 1);
		K(1, 1) = c*c*j11 - a*c*(j13 + j13) + a*a*j33 + invMass;
		K(1, 2) = -(b*c*j11) + a*c*j12 + a*b*j13 - a*a*j23;
		K(2, 0) = K(0, 2);
		K(2, 1) = K(1, 2);
		K(2, 2) = b*b*j11 - a*b*(j12 + j12) + a*a*j22 + invMass;
	}
	else
		K.setZero();
}

void Dlo::crossProductMatrix(const Eigen::Matrix<Real,3,1> &v, Eigen::Matrix<Real,3,3> &v_hat){
	v_hat <<     0, -v(2),  v(1),
              v(2),     0, -v(0),
             -v(1),  v(0),     0;
    }

void Dlo::postSolve(const Real &dt){
    // Update velocities
    #pragma omp parallel default(shared)
    {
        // Update linear velocities
        #pragma omp for schedule(static) 
        // #pragma omp parallel for 
        for (int i = 0; i< num_particles_; i++){
            if (is_dynamic_[i]){
                vel_.col(i) = (pos_[i] - prev_pos_[i])/dt;
            }
        }
        // Update angular velocities
        #pragma omp for schedule(static) 
        // #pragma omp parallel for 
        for (int i = 0; i< num_quaternions_; i++){
            if (is_dynamic_[i]){
                const Eigen::Quaternion<Real> relRot = (ori_[i] * prev_ori_[i].conjugate());
                omega_.col(i) = (relRot.w() >= 0) ? (2.0*relRot.vec()/dt) : (-2.0*relRot.vec()/dt);
            }
        }
    }

    // Create an artificial global damping
    #pragma omp parallel default(shared)
    {
        // Damp linear velocities
        #pragma omp for schedule(static) 
        // #pragma omp parallel for 
        for (int i = 0; i< num_particles_; i++){
            if (is_dynamic_[i]){
                vel_.col(i) -= std::min(1.0, (global_damp_coeff_v_/num_particles_)*dt*inv_mass_[i]) * vel_.col(i);
                // divide damping coeff by num_particles_ to get rid of the segment number dependent damping response
            }
        }
        // Damp angular velocities
        #pragma omp for schedule(static) 
        // #pragma omp parallel for 
        for (int i = 0; i< num_quaternions_; i++){
            if (is_dynamic_[i]){
                Eigen::Matrix<Real,3,1> dw = (global_damp_coeff_w_/num_quaternions_)*dt*inv_iner_[i]*omega_.col(i);
                // divide damping coeff by num_quaternions_ to get rid of the segment number dependent damping response
                
                if(dw.norm() >= omega_.col(i).norm()){
                    omega_.col(i).setZero();
                } else {
                    omega_.col(i) -= dw;
                }
            }
        }
    }
}

void Dlo::resetForces(){
    // Also Reset accumulated forces for the next iteration
    for_.setZero();
    tor_.setZero();
}

void Dlo::normalizeForces(const int &num_substeps){
    // Find the average forces of the timestep during the substep iterations
    for_ /= num_substeps;
    tor_ /= num_substeps;
}

// Find the nearest 3D position vector col id in the given matrix
int Dlo::findNearestPositionVectorId(const std::vector<Eigen::Matrix<Real,3,1>>& matrix, 
                                     const Eigen::Matrix<Real,3,1>& pos) {
  int nearestId = -1;
  Real minDistance = std::numeric_limits<Real>::max();
  for (int i = 0; i < matrix.size(); ++i) {
    const Eigen::Matrix<Real,3,1>& currentPos = matrix[i];
    Real currentDistance = (currentPos - pos).norm();
    if (currentDistance < minDistance) {
      nearestId = i;
      minDistance = currentDistance;
    }
  }
  return nearestId;
}

int Dlo::attachNearest(const Eigen::Matrix<Real,3,1> &pos){
    int id = findNearestPositionVectorId(pos_,pos);
    // Make that particle stationary
    if (id >= 0){
        changeParticleDynamicity(id,false);
    }
    return id;
}

void Dlo::updateAttachedPose(const int &id, 
                             const Eigen::Matrix<Real,3,1> &pos, 
                             const Eigen::Quaternion<Real> &ori ){
    pos_[id] = pos;
    ori_[id] = ori;
}

void Dlo::updateAttachedVelocity(const int &id, 
                                const Eigen::Matrix<Real,3,1> &vel, 
                                const Eigen::Matrix<Real,3,1> &omega){
    vel_.col(id) = vel;
    omega_.col(id) = omega;
}

Eigen::Matrix2Xi *Dlo::getStretchBendTwistIdsPtr(){
    return &stretchBendTwist_ids_;
}

std::vector<Eigen::Matrix<Real,3,1>> *Dlo::getPosPtr(){
    return &pos_;
}

std::vector<Eigen::Matrix<Real,3,1>> *Dlo::getPrevPosPtr(){
    return &prev_pos_;
}

Eigen::Matrix<Real,3,Eigen::Dynamic> *Dlo::getVelPtr(){
    return &vel_;
}

Eigen::Matrix<Real,3,Eigen::Dynamic> *Dlo::getForPtr(){
    return &for_;
}

std::vector<Real> *Dlo::getInvMassPtr(){
    return &inv_mass_;
}

std::vector<Eigen::Quaternion<Real>> *Dlo::getOriPtr(){
    return &ori_;
}

std::vector<Eigen::Quaternion<Real>> *Dlo::getPrevOriPtr(){
    return &prev_ori_;
}

Eigen::Matrix<Real,3,Eigen::Dynamic> *Dlo::getAngVelPtr(){
    return &omega_;
}

Eigen::Matrix<Real,3,Eigen::Dynamic> *Dlo::getTorPtr(){
    return &tor_;
}

std::vector<Eigen::Matrix<Real,3,3>> *Dlo::getInvInerPtr(){
    return &inv_iner_;
}

std::vector<Eigen::Matrix<Real,3,3>> *Dlo::getInerPtr(){
    return &iner_;
}

std::vector<Real> *Dlo::getSegmentLengthsPtr(){
    return &mesh_.segment_lengths;
}

Real Dlo::getMaxError(){
    return max_error_;
}

std::vector<int> *Dlo::getAttachedIdsPtr(){
    return &attached_ids_;
}

