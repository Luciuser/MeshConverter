#include "remesh_parameter.h"

#include "meshAlgorithm.h"

MESHIO::RemeshParameter::RemeshParameter()
{
    size_function_ = [](double x, double y, double z)->double {
        return 0.0;
    };
}

MESHIO::RemeshParameter::~RemeshParameter()
{
    if (aabbTree_ != nullptr) {
        delete aabbTree_;
        aabbTree_ = nullptr;
    }
}

int MESHIO::RemeshParameter::buildAABBTree(Mesh& mesh)
{
    if (aabbTree_ != nullptr) {
        delete aabbTree_;
        aabbTree_ = nullptr;
    }

    // build AABB tree
    double hmax, hmin, average;
    MESHIO::calculateEdgesLength(mesh, hmax, hmin, average);
    aabbTree_ = new aabb::Tree(3, hmin * 0.0001);

    // set data in
    int nFacet = mesh.Topo.rows();

    for (int i = 0; i < nFacet; i++) {
        std::vector<double> lower;
        std::vector<double> upper;
        MESHIO::calculateBoundingBoxForOneElement(mesh, i, lower, upper);
        aabbTree_->insertParticle(i, lower, upper);
    }

    b_use_aabb_tree_ = true;

    return 0;
}

int MESHIO::RemeshParameter::setAABBTreeActive()
{
    if (aabbTree_ == nullptr) {
        std::cout << "Warnning: no AABB tree which has bound while setting aabb tree active" << std::endl;
        return -1;
    }
    else {
        b_use_aabb_tree_ = true;
    }

    return 0;
}

int MESHIO::RemeshParameter::setAABBTreeDead()
{
    b_use_aabb_tree_ = false;
    return 0;
}

double MESHIO::RemeshParameter::getVertexSizeOrTargetLow(const Eigen::Vector3d& vertex)
{
    return this->b_use_size_function_ ? low_ratio_ * getVertexSize(vertex) : low_ratio_ * target_length_;
}

double MESHIO::RemeshParameter::getVertexSizeOrTargetHigh(const Eigen::Vector3d& vertex)
{
    return this->b_use_size_function_ ? high_ratio_ * getVertexSize(vertex) : high_ratio_ * target_length_;
}


double MESHIO::RemeshParameter::getVertexSize(const Eigen::Vector3d &vertex)
{
    double size = size_function_(vertex.x(), vertex.y(), vertex.z());
    if (size < hmin_) {
        return hmin_;
    }
    else if (size > hmax_) {
        return hmax_;
    }
    
    return size;
}