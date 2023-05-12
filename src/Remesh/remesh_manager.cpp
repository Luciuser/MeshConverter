#include "remesh_manager.h"

#include<algorithm>

#include "meshAlgorithm.h"
#include "tiger_sizingfunction.h"

#include "intersection.h"

#define AABB_TREE // TODO(jin): still have many bugs

using namespace MESHIO::polymesh;


MESHIO::RemeshManager::RemeshManager()
{
}

MESHIO::RemeshManager::~RemeshManager()
{
    //mesh_ = nullptr;
    if (abtree_ != nullptr) {
        delete abtree_;
        abtree_ = nullptr;
    }
}

int MESHIO::RemeshManager::checkInput(const std::vector<double>& points, const std::vector<int>& triangles, const std::vector<int>& surfaceID)
{
    if (points.size() % 3 != 0) {
        std::cout << "Error in input: The point coordinates array length is not n*3." << std::endl;
        return -1;
    }
    if (triangles.size() % 3 != 0) {
        std::cout << "Error in input: The triangle topology array length is not n*3." << std::endl;
        return -1;
    }

    int nPoints = points.size() / 3;
    int nFacets = triangles.size() / 3;
    int nSurfaceID = surfaceID.size();
    //std::cout << "nPoints " << nPoints << " nFacets " << nFacets << " nSurfaceID " << nSurfaceID << std::endl;
    //std::cout << "Points: " << std::endl;
    //for (int i = 0; i < nPoints; i++) {
    //    for (int j = 0; j < 3; j++) {
    //        std::cout << points[3 * i + j] << " ";
    //    }
    //    std::cout << std::endl;
    //}
    //std::cout << "Facets: " << std::endl;
    //for (int i = 0; i < nFacets; i++) {
    //    for (int j = 0; j < 3; j++) {
    //        std::cout << triangles[3 * i + j] << " ";
    //    }
    //    std::cout << std::endl;
    //}
    //std::cout << "SurfaceID: " << std::endl;
    //for (int i = 0; i < nSurfaceID; i++) {
    //    std::cout << surfaceID[i] << " " << std::endl;
    //}

    if (nPoints == 0) {
        std::cout << "Error in input: The point coordinates array length is 0." << std::endl;
        return -1;
    }

    if (nFacets == 0) {
        std::cout << "Error in input: The triangle topology array length is 0." << std::endl;
        return -1;
    }

    if (nFacets != nSurfaceID) {
        std::cout << "Error in input: The triangle topology array length is different with surfaceID array length." << std::endl;
        return -1;
    }
    
    for (int i = 0; i < triangles.size(); i++) {
        if (triangles[i] < 0) {
            std::cout << "Error in input: There is one index which is less than 0 in triangle topology array." << std::endl;
            return -1;
        }
        else if (triangles[i] >= nPoints) {
            std::cout << "Error in input: There is one index which is more than points num in triangle topology array." << std::endl;
            return -1;
        }
    }

    return 0;
}

int MESHIO::RemeshManager::getSizeFunction(const std::string filename, std::function<double(double, double, double)>& size_function)
{
    int sfObjID;
    int success = API_Create_SurfBKG_SF(
        filename.c_str(),
        1.2,
        &sfObjID
    );

    size_function = API_Sizing_Query;

    return 0;
}

int MESHIO::RemeshManager::initial(const Mesh& mesh)
{
    // mesh initial
    half_mesh_.clear();
    for (int i = 0; i < mesh.Vertex.rows(); ++i) {
        half_mesh_.addVertex(mesh.Vertex.row(i), Eigen::Vector2d(0, 0));
    }
    for (int i = 0; i < mesh.Topo.rows(); ++i) {
        auto t = std::vector<size_t>{ size_t(mesh.Topo(i, 0)), size_t(mesh.Topo(i, 1)), size_t(mesh.Topo(i, 2)) };
        half_mesh_.addPolyFace(t);
    }

    // AABB_tree initial
    if (abtree_ == nullptr) {
        delete abtree_;
    }
    abtree_ = nullptr;
    get_AABB_tree(&half_mesh_, abtree_);

    // parameter initial
    parameter_.setTargetLength(calculateTargetEdgeLength(&half_mesh_) / 2.0);
    return 0;
}

int MESHIO::RemeshManager::addSizeFunction(std::function<double(double, double, double)> size_function)
{
    parameter_.setSizeFunction(size_function);
    return 0;
}

int MESHIO::RemeshManager::deleteSizeFunction()
{
    parameter_.deleteSizeFunction();
    return 0;
}

double MESHIO::RemeshManager::calculateTargetEdgeLength(PolyMesh* mesh)
{
    double target_edge_length = 0.0;
    std::cout << mesh->numEdges() << std::endl;
    for (EdgeIter e_it = mesh->edges_begin(); e_it != mesh->edges_end(); ++e_it)
    {
        target_edge_length += (*e_it)->length();
    }
    target_edge_length /= mesh->numEdges();
    return target_edge_length;
}

int MESHIO::RemeshManager::getMesh(Mesh& mesh)
{
    mesh.Vertex.resize(half_mesh_.numVertices(), 3);
    int v_id = 0, f_id = 0;
    for (VertexIter iter = half_mesh_.vertices_begin(); iter != half_mesh_.vertices_end(); iter++) {
        mesh.Vertex.row(v_id++) = (*iter)->position();
    }
    mesh.Topo.resize(half_mesh_.numPolygons(), 3);
    mesh.Masks.resize(half_mesh_.numPolygons(), 1);
    for (FaceIter iter = half_mesh_.polyfaces_begin(); iter != half_mesh_.polyfaces_end(); iter++) {
        std::vector<int> index;
        for (FaceVertexIter fv_it = half_mesh_.fv_iter(*iter); fv_it.isValid(); ++fv_it)
        {
            MVert* fv = *fv_it;
            index.push_back(fv->index());
        }
        mesh.Topo.row(f_id) = Eigen::Vector3i(index[0], index[1], index[2]);
        mesh.Masks(f_id++, 0) = 0;
    }

    return 0;
}

int MESHIO::RemeshManager::repair(Mesh& mesh, double eps)
{
    MESHIO::removeDulplicatePoint(mesh.Vertex, mesh.Topo, eps);
    MESHIO::removeDegradationTopo(mesh);
    MESHIO::removeDulplicateTopo(mesh);
    MESHIO::removeHangingPoint(mesh);
    return 0;
}

int MESHIO::RemeshManager::buildAABBTree(Mesh& mesh)
{
    return parameter_.buildAABBTree(mesh);
}

int MESHIO::RemeshManager::buildAABBTree(Mesh& mesh, aabb::Tree** tree)
{
    // build AABB tree
    double hmax, hmin, average;
    MESHIO::calculateEdgesLength(mesh, hmax, hmin, average);
    (*tree) = new aabb::Tree(3, hmin * 0.0001);

    // set data in
    int nFacet = mesh.Topo.rows();

    for (int i = 0; i < nFacet; i++) {
        std::vector<double> lower;
        std::vector<double> upper;
        MESHIO::calculateBoundingBoxForOneElement(mesh, i, lower, upper);
        (*tree)->insertParticle(i, lower, upper);
    }

    return 0;
}

int MESHIO::RemeshManager::updateAABBTree(MPolyFace* polyface, aabb::Tree* tree)
{
    MHalfedge* he_begin = polyface->halfEdge();
    MHalfedge* he = he_begin;

    std::vector<Eigen::Vector3d> vertexes;

    //temp.push_back(he->fromVertex()->position());
    do {
        vertexes.push_back(he->toVertex()->position());
        he = he->next();
    } while (he != he_begin);

    std::vector<double> lower;
    std::vector<double> upper;
    calculateBoundingBoxForOneElement(vertexes, lower, upper);

    if (tree->findParticle(polyface->index())) {
        tree->updateParticle(polyface->index(), lower, upper);
    }
    else {
        tree->insertParticle(polyface->index(), lower, upper);
    }

    return 0;
}

int MESHIO::RemeshManager::calculateBoundingBoxForOneElement(const std::vector<Eigen::Vector3d>& point, std::vector<double>& lower, std::vector<double>& upper)
{
    lower.clear();
    upper.clear();

    lower.resize(3, std::numeric_limits<double>::max());
    upper.resize(3, -std::numeric_limits<double>::max());
    for (int i = 0; i < point.size(); i++) {
        if (lower[0] > point[i].x()) {
            lower[0] = point[i].x();
        }
        if (lower[1] > point[i].y()) {
            lower[1] = point[i].y();
        }
        if (lower[2] > point[i].z()) {
            lower[2] = point[i].z();
        }

        if (upper[0] < point[i].x()) {
            upper[0] = point[i].x();
        }
        if (upper[1] < point[i].y()) {
            upper[1] = point[i].y();
        }
        if (upper[2] < point[i].z()) {
            upper[2] = point[i].z();
        }
    }

    return 0;
}

int MESHIO::RemeshManager::triangleIntersectionCheck(Eigen::Vector3d& a, Eigen::Vector3d& b, Eigen::Vector3d& c, Eigen::Vector3d& o, Eigen::Vector3d& p, Eigen::Vector3d& q)
{
    return tri_tri_inter(a.data(), b.data(), c.data(), o.data(), p.data(), q.data());
}

int MESHIO::RemeshManager::remesh(const Mesh& mesh, Mesh& out)
{
    initial(mesh);
    remesh();
    getMesh(out);

    return 0;
}

int MESHIO::RemeshManager::remesh(const Mesh& mesh, std::function<double(double, double, double)>& size_function, Mesh& out)
{
    initial(mesh);
    addSizeFunction(size_function);
    remesh();
    getMesh(out);

    return 0;
}

int MESHIO::RemeshManager::remesh()
{
    for (int i = 0; i < 10; i++)
    {
        std::cout << "    Remesh in " << i << "th" << std::endl;
        split_long_edges(&half_mesh_, parameter_);
        collapse_short_edges(&half_mesh_, parameter_);
        equalize_valences(&half_mesh_, parameter_);
        tangential_relaxation(&half_mesh_);
        project_to_surface(&half_mesh_, abtree_);
    }
    return 0;
}

void MESHIO::RemeshManager::split_long_edges(PolyMesh* mesh, RemeshParameter& parameter)
{
    int NE = mesh->numEdges();
    for (int i = 0; i < NE; i++)
    {
        MEdge* e = mesh->edge(i);
        double len = e->length();
        Eigen::Vector3d middle_point = e->getCenter();
        double size = parameter.getVertexSizeOrTargetHigh(middle_point);
        if (len > size)
        {
            auto vt = mesh->splitEdgeTriangle(e);

            // use aabb tree to avoid intersection of triangles
            if (parameter.getAABBTreeActive()) {
                // find all polygen new
                MHalfedge* he_begin = vt->halfEdge();
                if (he_begin->toVertex() == vt) {
                    he_begin = he_begin->next();
                }
                if (he_begin->fromVertex() != vt) {
                    std::cout << "Error: error in split because of wrong topology." << std::endl;
                    return;
                }
                MHalfedge* he = he_begin;
                std::vector<MPolyFace*> adj_faces;
                do {
                    adj_faces.push_back(he->polygon());
                    do {
                        he = he->next();
                    } while (he->toVertex() != vt);
                    he = he->pair();

                } while (he != he_begin);

                // update all polygen which may be changed
                for (int i = 0; i < adj_faces.size(); i++) {
                    updateAABBTree(adj_faces[i], parameter_.aabbTree());
                }
            }

        }
    }
}

void MESHIO::RemeshManager::collapse_short_edges(PolyMesh* mesh, RemeshParameter& parameter)
{
    int NE = mesh->numEdges();
    for (int i = NE - 1; i >= 0; i--)
    {
        if (i > mesh->numEdges() - 1)
            continue;
        MEdge* e = mesh->edge(i);
        MHalfedge* he = e->halfEdge();
        MVert* p0 = he->fromVertex();
        MVert* p1 = he->toVertex();
        if (!mesh->is_collapse_ok(he))
            continue;
        if (mesh->isBoundary(p0) || mesh->isBoundary(p1))
            continue;
        double len = e->length();
        Eigen::Vector3d middle_point = e->getCenter();
        if (len < parameter.getVertexSizeOrTargetLow(middle_point))
        {
            bool is_collapse = true;
            for (VertexVertexIter vv_it = mesh->vv_iter(p0); vv_it.isValid(); ++vv_it)
            {
                MVert* vv = *vv_it;
                double len = (p1->position() - vv->position()).norm();
                Eigen::Vector3d middle_point = p1->position() * 0.5 + vv->position() * 0.5;
                if (len > parameter.getVertexSizeOrTargetHigh(middle_point))
                {
                    is_collapse = false;
                    break;
                }
            }
            if (is_collapse)
                mesh->collapseTriangle(he);
        }
    }
}

void MESHIO::RemeshManager::delete_lowdegree(PolyMesh* mesh)
{
    int NE = mesh->numEdges();
    for (int i = NE - 1; i >= 0; i--)
    {
        if (i > mesh->numEdges() - 1)
            continue;

        MEdge* e = mesh->edge(i);
        MHalfedge* he = e->halfEdge();
        MVert* p0 = he->fromVertex();
        MVert* p1 = he->toVertex();

        if (
            (!mesh->isBoundary(e)) && (!mesh->isBoundary(p0)) && (!mesh->isBoundary(p1)) &&
            ((mesh->valence(p0) == 3) || (mesh->valence(p1) == 3)))
        {
            // std::cout << "delete non-boundary 3 degree point" << std::endl;
            // std::cout << "collapse edge p0-p1£º "
            //<< p0->index() << " " << p1->index() << std::endl;

            mesh->collapseTriangle(he);
        }
    }
}

void MESHIO::RemeshManager::equalize_valences(PolyMesh* mesh, RemeshParameter& parameter)
{
    std::vector<int> target_valence;
    int deviation_pre, deviation_post;
    for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it)
    {
        if (mesh->isBoundary(*v_it))
            target_valence.push_back(4);
        else
            target_valence.push_back(6);
    }

    for (EdgeIter e_it = mesh->edges_begin(); e_it != mesh->edges_end(); ++e_it)
    {
        if (mesh->isBoundary(*e_it) || !mesh->is_flip_ok_Triangle(*e_it))
            continue;

        MHalfedge* he1 = (*e_it)->halfEdge(); // edge 1->2
        MVert* v1 = he1->fromVertex();
        MVert* v2 = he1->toVertex();

        MHalfedge* he2 = (*e_it)->halfEdge()->next(); // edge 2->3
        MVert* v3 = he2->toVertex();
        MHalfedge* he3 = (*e_it)->halfEdge()->pair()->next(); // edge 1->4
        MVert* v4 = he3->toVertex();

        /*angle judge*/
        // *** notation(original configuration)
        //
        //                        1            4
        //                       *---------*
        //						/ \		        /
        //					   /    \  E2  /
        //					  /       \    /
        //				     /	E1	    \ /
        //				    *---------*
        //                  3         2
        //
        //      Diagonal swapping. '
        /**
         *        a
         *      /
         *    /
         * b ----- c
         * calculate the angle of (ba and bc)
         */

        auto vectorAngle = [](const Eigen::Vector3d& a,
            const Eigen::Vector3d& b,
            const Eigen::Vector3d& c)
        {
            Eigen::Vector3d ba = (a - b).normalized();
            Eigen::Vector3d bc = (c - b).normalized();
            double cos_theta = ba.dot(bc);
            double ang = std::acos(cos_theta) * 180.0 / M_PI;
            return ang;
        };

        double a123 = vectorAngle(v1->position(), v2->position(), v3->position());
        double a231 = vectorAngle(v2->position(), v3->position(), v1->position());
        double a142 = vectorAngle(v1->position(), v4->position(), v2->position());
        double a421 = vectorAngle(v4->position(), v2->position(), v1->position());
        double a312 = 180 - a123 - a231;
        double a214 = 180 - a142 - a421;
        std::vector<double> BeforeSwap{ a123, a231, a142, a421, a312, a214 };

        auto MinAngle_before = std::min_element(BeforeSwap.begin(), BeforeSwap.end());

        double a234 = vectorAngle(v2->position(), v3->position(), v4->position());
        double a423 = vectorAngle(v4->position(), v2->position(), v3->position());
        double a143 = vectorAngle(v1->position(), v4->position(), v3->position());
        double a314 = vectorAngle(v3->position(), v1->position(), v4->position());
        double a342 = 180 - a234 - a423;
        double a431 = 180 - a143 - a314;
        std::vector<double> AfterSwap{ a234, a423, a143, a314, a342, a431 };

        auto MinAngle_after = std::min_element(AfterSwap.begin(), AfterSwap.end());

        // if the topology will make intersection like this, this operation shouble be banned
        //
        //                    1            
        //                    *
        //			        /  |  \ 
        //			      /    |    \
        //		        /     *      \
        //            /    / 2 \    \
        //           /  /           \ \
        //         * 3               4*
        if (a123 + a421 > 180 || a312 + a214 > 180) {
            continue;
        }

        // if the min angle get smaller after flip, this operation should be banned!
        if ((*MinAngle_after) < (*MinAngle_before))
            continue;

        /*degree judge*/
        deviation_pre = abs(int(mesh->valence(v1) - target_valence[v1->index()])) + abs(int(mesh->valence(v2) - target_valence[v2->index()])) + abs(int(mesh->valence(v3) - target_valence[v3->index()])) + abs(int(mesh->valence(v4) - target_valence[v4->index()]));
        deviation_post = abs(int(mesh->valence(v1) - 1 - target_valence[v1->index()])) + abs(int(mesh->valence(v2) - 1 - target_valence[v2->index()])) + abs(int(mesh->valence(v3) + 1 - target_valence[v3->index()])) + abs(int(mesh->valence(v4) + 1 - target_valence[v4->index()]));

        // non-boundary point shouldn't have degree less than 3
        if ((!mesh->isBoundary(v1)) && (mesh->valence(v1) - 1 <= 3))
            continue;
        if ((!mesh->isBoundary(v2)) && (mesh->valence(v2) - 1 <= 3))
            continue;
        if ((!mesh->isBoundary(v3)) && (mesh->valence(v3) + 1 <= 3))
            continue;
        if ((!mesh->isBoundary(v4)) && (mesh->valence(v4) + 1 <= 3))
            continue;

        if (deviation_pre < deviation_post) {
            continue;
        }

        // use aabb tree to judge whether exists intersection
        std::vector<double> lower;
        std::vector<double> upper;

        std::vector<Eigen::Vector3d> triangle143 = { v1->position(), v4->position(), v3->position() };
        std::vector<Eigen::Vector3d> triangle234 = { v2->position(), v3->position(), v4->position() };
        if (parameter_.aabbTree()) {
            /* aabb tree judge */

            // get simple aabb box intersection by aabb tree
            calculateBoundingBoxForOneElement(triangle143, lower, upper);
            parameter_.aabbTree()->insertParticle(UINT_MAX, lower, upper);
            auto intersection143 = parameter_.aabbTree()->query(UINT_MAX);
            parameter_.aabbTree()->removeParticle(UINT_MAX);
            std::sort(intersection143.begin(), intersection143.end());
            if (intersection143.size() >= 100) {
                std::cout << "Warning: the triangles intersection num of aabb box is " << intersection143.size() << std::endl;
            }

            calculateBoundingBoxForOneElement(triangle234, lower, upper);
            parameter_.aabbTree()->insertParticle(UINT_MAX, lower, upper);
            auto intersection234 = parameter_.aabbTree()->query(UINT_MAX);
            parameter_.aabbTree()->removeParticle(UINT_MAX);
            if (intersection234.size() >= 100) {
                std::cout << "Warning: the triangles intersection num of aabb box is " << intersection234.size() << std::endl;
            }

            // try to find intersection robust
            bool b_intersect = false;
            for (auto iter : intersection143) {
                if (b_intersect) {
                    break;
                }
                if (iter != UINT_MAX && iter != he1->polygon()->index() && iter != he1->pair()->polygon()->index()) {
                    MHalfedge* iter_he_begin = getPolyFace(iter)->halfEdge();
                    MHalfedge* iter_he = iter_he_begin;
                    std::vector<MVert*> iter_poly_vertexes;
                    do
                    {
                        MVert* iter_fv = iter_he->toVertex();
                        iter_poly_vertexes.push_back(iter_fv);

                        iter_he = iter_he->next();
                    } while (iter_he != iter_he_begin);

                    if (iter_poly_vertexes.size() < 3) {
                        std::cout << "Warning: there exist one poly which only has two points while swapping." << std::endl;
                        b_intersect = true;
                        continue;
                    }

                    int success = triangleIntersectionCheck(v1->position(), v4->position(), v3->position(), iter_poly_vertexes[0]->position(), iter_poly_vertexes[1]->position(), iter_poly_vertexes[2]->position());
                    if (success != interresult::DISJOINT) {
                        b_intersect = true;
                    }
                }
            }
            for (auto iter : intersection234) {
                if (b_intersect) {
                    break;
                }
                if (iter != UINT_MAX && iter != he1->polygon()->index() && iter != he1->pair()->polygon()->index()) {
                    MHalfedge* iter_he_begin = getPolyFace(iter)->halfEdge();
                    MHalfedge* iter_he = iter_he_begin;
                    std::vector<MVert*> iter_poly_vertexes;
                    do
                    {
                        MVert* iter_fv = iter_he->toVertex();
                        iter_poly_vertexes.push_back(iter_fv);

                        iter_he = iter_he->next();
                    } while (iter_he != iter_he_begin);

                    if (iter_poly_vertexes.size() < 3) {
                        std::cout << "Warning: there exist one poly which only has two points while swapping." << std::endl;
                        b_intersect = true;
                        continue;
                    }

                    int success = triangleIntersectionCheck(v1->position(), v4->position(), v3->position(), iter_poly_vertexes[0]->position(), iter_poly_vertexes[1]->position(), iter_poly_vertexes[2]->position());
                    if (success != interresult::DISJOINT) {
                        b_intersect = true;
                    }
                }
            }

            // if there exist intersection after flip operation, this operation should be banned!
            if (b_intersect) {
                continue;
            }
        }


        // finish all test
        if (parameter_.aabbTree()) {
            // update aabb tree: the same with the flip edge topology
            calculateBoundingBoxForOneElement(triangle234, lower, upper);
            parameter_.aabbTree()->updateParticle(he1->polygon()->index(), lower, upper);
            calculateBoundingBoxForOneElement(triangle143, lower, upper);
            parameter_.aabbTree()->updateParticle(he1->pair()->polygon()->index(), lower, upper);
        }

        // face 1 2 3 ==> face 2 3 4
        // face 1 4 2 ==> face 1 4 3
        mesh->flipEdgeTriangle(*e_it);
    }

    /*boundary judge*/
    // if two points are boundary but the edge is none boundary,then must flip!
    // for (EdgeIter e_it = mesh->edges_begin(); e_it != mesh->edges_end(); ++e_it)
    //{
    //	if (mesh->isBoundary(*e_it) || !mesh->is_flip_ok_Triangle(*e_it)) continue;

    //	MHalfedge* he1 = (*e_it)->halfEdge(); //edge 1->2
    //	MVert* v1 = he1->fromVertex();
    //	MVert* v2 = he1->toVertex();

    //	MHalfedge* he2 = (*e_it)->halfEdge()->next();//edge 2->3
    //	MVert* v3 = he2->toVertex();
    //	MHalfedge* he3 = (*e_it)->halfEdge()->pair()->next();//edge 1->4
    //	MVert* v4 = he3->toVertex();

    //	if ((mesh->isBoundary(v1) && mesh->isBoundary(v2))
    //		&& (mesh->isBoundary(v3) && mesh->isBoundary(v4)))
    //		continue;

    //	if ((!mesh->isBoundary(*e_it)) && (mesh->isBoundary(v1) && mesh->isBoundary(v2)))

    //		mesh->flipEdgeTriangle(*e_it);
    //}
}

void MESHIO::RemeshManager::tangential_relaxation(PolyMesh* mesh)
{
    mesh->updateMeshNormal();

    auto get_TriFace_Area = [](const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2)
    {
        return sqrt(0.5 * ((v1 - v0).cross((v2 - v0))).squaredNorm());
    };

    for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it)
    {
        if (mesh->isBoundary(*v_it))
            continue;
        double count = 0.0;
        Eigen::Vector3d p = (*v_it)->position();
        Eigen::Vector3d q(0.0, 0.0, 0.0);

        // CPT
        double area = 0, sum_area = 0;
        for (VertexFaceIter vf_it = mesh->vf_iter(*v_it); vf_it.isValid(); ++vf_it)
        {
            std::vector<Eigen::Vector3d> vertexlist;
            for (FaceVertexIter fv_it = mesh->fv_iter(*vf_it); fv_it.isValid(); ++fv_it)
            {
                vertexlist.push_back((*fv_it)->position());
            }

            area = get_TriFace_Area(vertexlist[0], vertexlist[1], vertexlist[2]);
            sum_area += area;

            q += ((*vf_it)->getFaceCenter()) * area;
        }
        q /= sum_area;

        // Laplace
        // for (VertexVertexIter vv_it = mesh->vv_iter(*v_it); vv_it.isValid(); ++vv_it)
        //{
        //	q += (*vv_it)->position();
        //	++count;
        // }
        // q /= count;
        Eigen::Vector3d n = (*v_it)->normal();
        n.normalize();

        Eigen::Vector3d newPos = q + (n.dot(p - q)) * n;

        //// For continue surface.
        // Eigen::Vector2d res_uv = {(*v_it)->u(), (*v_it)->v()};
        // double dis;
        // surface->OrthogonalProjection(newPos, 1e-5, res_uv, dis);
        // surface->param_to_coord(res_uv, newPos);

        (*v_it)->setPosition(newPos);
    }
}

void MESHIO::RemeshManager::get_AABB_tree(PolyMesh* mesh, AABB_Tree*& abtree)
{
    std::vector<Vector3f> point_set;
    point_set.clear();
    for (FaceIter f_it = mesh->polyfaces_begin(); f_it != mesh->polyfaces_end(); ++f_it)
    {
        for (FaceVertexIter fv_it = mesh->fv_iter(*f_it); fv_it.isValid(); ++fv_it)
        {
            MVert* fv = *fv_it;
            Vector3f p;
            p[0] = float(fv->x());
            p[1] = float(fv->y());
            p[2] = float(fv->z());
            point_set.push_back(p);
        }
    }
    abtree = new AABB_Tree(point_set);
}

void MESHIO::RemeshManager::project_to_surface(PolyMesh* mesh, AABB_Tree* abtree)
{
    for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it)
    {
        Vector3f p;
        p[0] = float((*v_it)->x());
        p[1] = float((*v_it)->y());
        p[2] = float((*v_it)->z());
        Vector3f ab_nearst_point;
        abtree->findNearstPoint(p, ab_nearst_point);
        Eigen::Vector3d new_point;
        new_point[0] = double(ab_nearst_point[0]);
        new_point[1] = double(ab_nearst_point[1]);
        new_point[2] = double(ab_nearst_point[2]);
        (*v_it)->setPosition(new_point);
    }
}
