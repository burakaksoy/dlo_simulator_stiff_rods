/*
 * Author: Burak Aksoy
 */

#ifndef DLO_H
#define DLO_H

#include <math.h>
#include <cmath>
#include <numeric>
#include <vector>
#include <list>
#include <unordered_map> 
#include <algorithm>    // std::sort, std::min

#include <Eigen/Dense>
#include <Eigen/Geometry>
// #include <scipy/spatial.h>

#include <string>
#include <iostream> // std::cout

#include <float.h>

#include <omp.h>

#define USE_DOUBLE // comment out if you would like to use float.

#ifdef USE_DOUBLE
typedef double Real;
#else
typedef float Real;
#endif

namespace pbd_object
{

struct MeshDLO
{
    std::string name;
    std::vector<Eigen::Matrix<Real,3,1>> vertices;
    std::vector<Eigen::Quaternion<Real>> quaternions;
    std::vector<Real> segment_lengths;
};

/** Node in the simulated tree structure */
struct Node {
    Node() {
        object = NULL; D = Dinv = J = Eigen::Matrix<Real, 6, 6>::Zero(); parent = NULL;
        soln.setZero(); index = 0;
    };
    Node *parent;
    bool isconstraint;
    int index;

    void *object;
    
    Eigen::Matrix<Real, 6, 6> D, Dinv, J;
    
    std::vector <Node*> children;
    
    Eigen::Matrix<Real, 6, 1> soln;
    Eigen::LDLT<Eigen::Matrix<Real, 6, 6>> DLDLT;
};

class Dlo
{
public:
    Dlo();
    Dlo(const MeshDLO &mesh, 
        const Real &zero_stretch_stiffness, 
        const Real &young_modulus, 
        const Real &torsion_modulus, 
        const Real &density,
        const Real &radius,
        const bool &use_direct_kkt_solver,
        const bool &use_zero_stretch_stiffness,
        const Real &global_damp_coeff_v,
        const Real &global_damp_coeff_w);
    ~Dlo();

    void preSolve(const Real &dt, const Eigen::Matrix<Real,3,1> &gravity);
    void solve(const Real &dt);
    void postSolve(const Real &dt);

    void hangFromCorners(const int &num_corners);

    
    int attachNearest(const Eigen::Matrix<Real,3,1> &pos);
    void updateAttachedPose(const int &id, 
                            const Eigen::Matrix<Real,3,1> &pos, 
                            const Eigen::Quaternion<Real> &ori);

    /*
    void resetForces();
    */

    Eigen::Matrix2Xi *getStretchBendTwistIdsPtr();

    std::vector<Eigen::Matrix<Real,3,1>> *getPosPtr();
    Eigen::Matrix<Real,3,Eigen::Dynamic> *getVelPtr();
    Eigen::Matrix<Real,3,Eigen::Dynamic> *getForPtr();

    std::vector<Eigen::Quaternion<Real>> *getOriPtr();
    std::vector<Real> *getSegmentLengthsPtr();

    Real getMaxError();
    

    std::vector<int> *getAttachedIdsPtr();    

private:
    // Functions
    void setStretchBendTwistConstraints();
    void setMasses();


    int findNearestPositionVectorId(const std::vector<Eigen::Matrix<Real,3,1>>& matrix, 
                                    const Eigen::Matrix<Real,3,1>& pos);
    

    void solveStretchBendTwistConstraints(const Real &dt);
    void computeMatrixK(const Eigen::Matrix<Real,3,1> &connector, 
                        const Real invMass, 
                        const Eigen::Matrix<Real,3,1> &x, 
                        const Eigen::Matrix<Real,3,3> &inertiaInverseW, 
                        Eigen::Matrix<Real,3,3> &K);
    void crossProductMatrix(const Eigen::Matrix<Real,3,1> &v, 
                            Eigen::Matrix<Real,3,3> &v_hat);

    // Variables
    MeshDLO mesh_;
    
    Real density_; // dlo mass per meter cube (kg/m^3)
    Real radius_; // radius of dlo assuming it's cylndrical

    int num_particles_;
    int num_quaternions_;

    // Particle Data
    std::vector<Eigen::Matrix<Real,3,1>> pos_;
    std::vector<Eigen::Matrix<Real,3,1>> prev_pos_;
    std::vector<Eigen::Matrix<Real,3,1>> rest_pos_;
    Eigen::Matrix<Real,3,Eigen::Dynamic> vel_;
    Eigen::Matrix<Real,3,Eigen::Dynamic> for_;
    std::vector<Real> inv_mass_;

    // Orientation Data
    std::vector<Eigen::Quaternion<Real>> ori_;
    std::vector<Eigen::Quaternion<Real>> prev_ori_;
    std::vector<Eigen::Quaternion<Real>> rest_ori_;
    Eigen::Matrix<Real,3,Eigen::Dynamic> omega_;
    Eigen::Matrix<Real,3,Eigen::Dynamic> tor_;
    std::vector<Eigen::Matrix<Real,3,3>> inv_iner_;
    std::vector<Eigen::Matrix<Real,3,3>> iner_;

    // Zero-stretch bending twisting constraint info
    Eigen::Matrix2Xi stretchBendTwist_ids_; // stores segment id pair of each constraint
    std::vector<Eigen::Matrix<Real,3,1>> stretchBendTwist_restDarbouxVectors_;
    std::vector<Eigen::Matrix<Real,3,4>> stretchBendTwist_constraintPosInfo_;
    std::vector<Real> average_segment_lengths_;

    std::vector<Eigen::Matrix<Real,6,1>> RHS_;
    std::vector<std::vector<Eigen::Matrix<Real,3,3>>> bendingAndTorsionJacobians_;
    

    Real zero_stretch_stiffness_;
    Real young_modulus_;
    Real torsion_modulus_;


    // -----------------------------------------------------------------
    // BEGIN: Parameters and functions for Tree structured linear time direct solver
    Node* root_;
    std::list<Node*> * forward_; // from leave to the root
    std::list<Node*> * backward_; // from root to the leaves

    Real max_error_;

    void initTree();
    void initNodes();
    void initSegmentNode(Node *n);

    void orderMatrixH(Node *n);

    void factor(const Eigen::Matrix<Real,3,1> &stretch_compliance,
                const Eigen::Matrix<Real,3,1> &bending_and_torsion_compliance);

    void solver();

    bool use_direct_kkt_solver_;
    bool use_zero_stretch_stiffness_;

    Real global_damp_coeff_v_; 
    Real global_damp_coeff_w_;
    // END: Parameters and functions for Tree structured linear time direct solver
    // -----------------------------------------------------------------

    std::vector<int> attached_ids_; // ids of robot attached particles
    // int grab_id_;
    // Real grab_inv_mass_;
};

} // namespace pbd_object

#endif /* !DLO_H */