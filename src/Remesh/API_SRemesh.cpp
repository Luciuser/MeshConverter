#include "API_SRemesh.h"

#include "meshAlgorithm.h"
#include "tiger_sizingfunction.h"
#include "remesh_manager.h"

int API_remesh_non_manifold_file(
    // ----- input -----//
    const std::vector<double>& points,
    const std::vector<int>& triangles,
    const std::vector<int>& surfaceID,
    // ----- output -----//
    std::vector<double>& points_out,
    std::vector<int>& triangles_out,
    std::vector<int>& surfaceID_out,
    // ----- parameter -----//
    const int iterations,
    const std::string size_vtk_file_name,
    const double size_target,
    const bool b_useAABBtree,
    const bool b_split,
    const double dulplicate_point_eps,
    const bool b_write_vtk,
    const std::string write_vtk_file_name,
    const bool b_write_split_vtk)
{
    std::function<double(double, double, double)> size_function;

    if (size_vtk_file_name == "NO_USE") {
        size_function = nullptr;
    }
    else {
        MESHIO::RemeshManager RM;
        RM.getSizeFunction(size_vtk_file_name, size_function);
    }

    return API_remesh_non_manifold(
        points,
        triangles,
        surfaceID,
        points_out,
        triangles_out,
        surfaceID_out,
        iterations,
        size_function,
        size_target,
        b_useAABBtree,
        b_split,
        dulplicate_point_eps,
        b_write_vtk,
        write_vtk_file_name,
        b_write_split_vtk);
}

int API_remesh_non_manifold(
    // ----- input -----//
    const std::vector<double>& points,
    const std::vector<int>& triangles,
    const std::vector<int>& surfaceID,
    // ----- output -----//
    std::vector<double>& points_out,
    std::vector<int>& triangles_out,
    std::vector<int>& surfaceID_out,
    // ----- parameter -----//
    const int iterations,
    const std::function<double(double, double, double)>& size_function,
    const double size_target,
    const bool b_useAABBtree,
    const bool b_split,
    const double dulplicate_point_eps,
    const bool b_write_vtk,
    const std::string write_vtk_file_name,
    const bool b_write_split_vtk)
{
    MESHIO::RemeshManager RM;

    // mesh check
    if (RM.checkInput(points, triangles, surfaceID) != 0) {
        return -1;
    }
    std::cout << "    check input end" << std::endl;

    // mesh input
    Mesh mesh;
    MESHIO::setData(points, triangles, surfaceID, mesh);

    // mesh repair
    RM.repair(mesh, dulplicate_point_eps);
    MESHIO::resetOrientation(mesh);

    // parameter
    // iteration
    RM.setIteration(iterations);
    // size function
    if (size_function == nullptr) {
        if (size_target > 0) {
            RM.setUserTargetLength(size_target);
        }
    }
    else {
        RM.addSizeFunction(size_function);
    }

    // intesecetion aabb treee
    if (b_useAABBtree) {
        RM.buildAABBTree(mesh);
    }

    // core remesh
    if (!b_split) {
        // not split
        RM.remesh(mesh, mesh);
    }
    else {
        if (size_function == nullptr) {
            if (size_target <= 0) {
                RM.initialUserTargetLength(mesh);
            }
        }

        // split mesh to several small part
        std::vector<Mesh> split_mesh_list;
        MESHIO::splitDifferentFaces(mesh, split_mesh_list);

        for (int i = 0; i < split_mesh_list.size(); i++) {

            //std::string before_write_vtk_file_name = write_vtk_file_name;
            //before_write_vtk_file_name += "before_10.vtk";
            //MESHIO::writeVTK(before_write_vtk_file_name, split_mesh_list[i], "surface_id");

            // local parameter

            // core remesh
            std::cout << "remesh" << i << "th part " << std::endl;
            RM.remesh(split_mesh_list[i], split_mesh_list[i]);

            if (b_write_split_vtk) {
                std::cout << "write file to " << write_vtk_file_name << std::endl;
                std::string new_write_vtk_file_name = write_vtk_file_name;
                new_write_vtk_file_name += std::to_string(i) + ".vtk";
                MESHIO::writeVTK(new_write_vtk_file_name, split_mesh_list[i], "surface_id");
            }
        }

        MESHIO::addMesh(split_mesh_list, mesh);
    }

    std::cout << "remesh end" << std::endl;

    // mesh repair
    RM.repair(mesh, dulplicate_point_eps);

    // mesh output
    MESHIO::getData(mesh, points_out, triangles_out, surfaceID_out);

    if (b_write_vtk) {
        std::cout << "write file to " << write_vtk_file_name << std::endl;
        MESHIO::writeVTK(write_vtk_file_name, mesh, "surface_id");
    }

    return 0;
}