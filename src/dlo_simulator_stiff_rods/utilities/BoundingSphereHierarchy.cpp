#include "dlo_simulator_stiff_rods/utilities/BoundingSphereHierarchy.h"

#include <iostream>
#include <unordered_set>
#include <set>

using pool_set = std::set<unsigned int>;
using namespace utilities;

PointCloudBSH::PointCloudBSH() : super(0, 10) {
    m_vertices = nullptr;
    m_vertices_v = nullptr;
}

void PointCloudBSH::init(const Eigen::Matrix<Real, Eigen::Dynamic, 3> *vertices, 
                         const unsigned int numVertices) {
    m_lst.resize(numVertices);
    m_vertices = vertices;
}

void PointCloudBSH::init(const Eigen::Matrix<Real, 3, 1> *vertices, 
                        const unsigned int numVertices)
{
	m_lst.resize(numVertices);
	m_vertices_v = vertices;
}

const Eigen::Matrix<Real, 3, 1>& PointCloudBSH::entity_position(unsigned int i) const {
    if (m_vertices){
        m_tempPosition = m_vertices->row(i).transpose();
        return m_tempPosition;
    }
    else{
        return m_vertices_v[i];
    }
}

void PointCloudBSH::compute_hull(unsigned int b, 
                                 unsigned int n, 
                                 BoundingSphere& hull) const {
    auto vertices_subset = std::vector<Eigen::Matrix<Real, 3, 1>>(n);

    if (m_vertices){
        for (unsigned int i = b; i < b + n; ++i){
            vertices_subset[i - b] = m_vertices->row(m_lst[i]).transpose();
        }
    }
    else{
        for (unsigned int i = b; i < b + n; ++i){
            vertices_subset[i - b] = m_vertices_v[m_lst[i]];
        }
    }
    
    const BoundingSphere s(vertices_subset);

    hull.x() = s.x();
    hull.r() = s.r();
}

void PointCloudBSH::compute_hull_approx(unsigned int b, 
                                        unsigned int n, 
                                        BoundingSphere& hull) const {
	// compute center											
    Eigen::Matrix<Real, 3, 1> x;
    x.setZero();
    Real radius2 = 0.0;


    if (m_vertices){
        for (unsigned int i = b; i < b + n; i++) {
            x += m_vertices->row(m_lst[i]).transpose();
        }
        x /= static_cast<Real>(n);
        for (unsigned int i = b; i < b + n; i++) {
            radius2 = std::max(radius2, (x - m_vertices->row(m_lst[i]).transpose()).squaredNorm());
        }
    }
    else{
        for (unsigned int i = b; i < b + n; i++) {
            x += m_vertices_v[m_lst[i]];
        }
        x /= static_cast<Real>(n);
        for (unsigned int i = b; i < b + n; i++) {
            radius2 = std::max(radius2, (x - m_vertices_v[m_lst[i]]).squaredNorm());
        }
    }

    hull.x() = x;
    hull.r() = sqrt(radius2);
}
