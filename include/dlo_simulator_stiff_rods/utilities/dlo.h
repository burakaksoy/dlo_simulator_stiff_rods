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

struct Mesh
{
    // EXPLANATIONS OF THE DATA ROWS OF AN .OBJ FILE
    
    // vertices (vertex data):
    // The numbers following "v" are the x, y, and z coordinates of a vertex in 3D space.

    // face_tri_ids (face data):
    // Lines beginning with "f "define the faces of the geometry, using the vertices defined earlier.
    // Each set of numbers in a face definition refers to a vertex index and its corresponding texture coordinate index, separated by a slash (/). The indices start from 1.
    // For instance, "f 1/1 2/2 3/3" defines a face using vertices 1, 2, and 3, with their corresponding texture coordinates 1, 2, and 3.
    // These faces are typically triangles in OBJ files, but can also be other polygons.

    // tex_coords (texture coordinate data):
    // Lines beginning with "vt" represent texture coordinates, used to map a 2D texture onto a 3D model.
    // The numbers following "vt" are the u and v coordinates (similar to x and y) in the texture space, typically ranging from 0.0 to 1.0.

    // normals (vertex normal data):
    // stores the normal vectors for each vertex.
    // Lines beginning with vn represents the normal vector of a vertex, with x, y, and z components.
    // Normal vectors are essential for lighting calculations, as they indicate the direction a surface is facing.

    std::string name; // Name of the mesh, used for identification.

    // Eigen::MatrixX3d vertices;
    Eigen::Matrix<Real,Eigen::Dynamic,3> vertices; // for vertex data (v)
    Eigen::MatrixX3i face_tri_ids; // for face data (f)
    Eigen::Matrix<Real, Eigen::Dynamic, 2> tex_coords; // for texture coordinates (vt), 
    Eigen::Matrix<Real, Eigen::Dynamic, 3> normals; // for normal vectors on each vertex (vn) 
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
        const bool &use_zero_stretch_stiffness,
        const Real &global_damp_coeff_v,
        const Real &global_damp_coeff_w);
    ~Dlo();

    void preSolve(const Real &dt, const Eigen::Matrix<Real,3,1> &gravity);
    void solve(const Real &dt);
    void postSolve(const Real &dt);    

    void changeParticleDynamicity(const int &particle, const bool &is_dynamic);
    
    void setStaticParticles(const std::vector<int> &particles);
    void setDynamicParticles(const std::vector<int> &particles);
    
    const bool isStaticParticle(const int &particle);
    const bool isDynamicParticle(const int &particle);

    void setYoungModulus(const Real &young_modulus);
    void setTorsionModulus(const Real &torsion_modulus);

    const Real getYoungModulus();
    const Real getTorsionModulus();
    
    int attachNearest(const Eigen::Matrix<Real,3,1> &pos);
    void updateAttachedPose(const int &id, 
                            const Eigen::Matrix<Real,3,1> &pos, 
                            const Eigen::Quaternion<Real> &ori);

    void updateAttachedVelocity(const int &id, 
                                const Eigen::Matrix<Real,3,1> &vel, 
                                const Eigen::Matrix<Real,3,1> &omega);

    void resetForces();
    void normalizeForces(const int &num_substeps);

    Eigen::Matrix2Xi *getStretchBendTwistIdsPtr();

    std::vector<Eigen::Matrix<Real,3,1>> *getPosPtr();
    std::vector<Eigen::Matrix<Real,3,1>> *getPrevPosPtr();
    Eigen::Matrix<Real,3,Eigen::Dynamic> *getVelPtr();
    Eigen::Matrix<Real,3,Eigen::Dynamic> *getForPtr();
    std::vector<Real> *getInvMassPtr();

    std::vector<Eigen::Quaternion<Real>> *getOriPtr();
    std::vector<Eigen::Quaternion<Real>> *getPrevOriPtr();
    Eigen::Matrix<Real,3,Eigen::Dynamic> *getAngVelPtr();
    Eigen::Matrix<Real,3,Eigen::Dynamic> *getTorPtr();
    std::vector<Eigen::Matrix<Real,3,3>> *getInvInerPtr();
    std::vector<Eigen::Matrix<Real,3,3>> *getInerPtr();

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
    std::vector<bool> is_dynamic_; // vector that holds the data whether the particle is dynamic or not

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

    Real max_error_;

    bool use_zero_stretch_stiffness_;

    Real global_damp_coeff_v_; 
    Real global_damp_coeff_w_;
    
    std::vector<int> attached_ids_; // ids of robot attached particles
};

} // namespace pbd_object

#endif /* !DLO_H */