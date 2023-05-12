#include "API_SRemesh.h"

#include "meshAlgorithm.h"
#include "tiger_sizingfunction.h"
#include "remesh_manager.h"

int API_remesh_non_manifold(
    const std::vector<double>& points,
    const std::vector<int>& triangles,
    const std::vector<int>& surfaceID,
    const std::string size_vtk_file_name,
    std::vector<double>& points_out,
    std::vector<int>& triangles_out,
    std::vector<int>& surfaceID_out,
    const bool b_useAABBtree,
    const bool b_split,
    const bool b_write_vtk,
    const std::string write_vtk_file_name)
{
    MESHIO::RemeshManager RM;
    std::function<double(double, double, double)> size_function;

    RM.getSizeFunction(size_vtk_file_name, size_function);

    return API_remesh_non_manifold(
        points,
        triangles,
        surfaceID,
        size_function,
        points_out,
        triangles_out,
        surfaceID_out,
        b_useAABBtree,
        b_split,
        b_write_vtk,
        write_vtk_file_name);
}

int API_remesh_non_manifold(
    const std::vector<double>& points,
    const std::vector<int>& triangles,
    const std::vector<int>& surfaceID,
    const std::function<double(double, double, double)>& size_function,
    std::vector<double>& points_out,
    std::vector<int>& triangles_out,
    std::vector<int>& surfaceID_out,
    const bool b_useAABBtree,
    const bool b_split,
    const bool b_write_vtk,
    const std::string write_vtk_file_name)
{
    if (!b_split) {
        // NOT split
        MESHIO::RemeshManager RM;

        // mesh check
        if (RM.checkInput(points, triangles, surfaceID) != 0) {
            return -1;
        }
        std::cout << "    check input end" << std::endl;

        // mesh input
        Mesh mesh;
        MESHIO::setData(points, triangles, surfaceID, mesh);

        // parameter
        RM.addSizeFunction(size_function);

        // intesecetion aabb treee
        if (b_useAABBtree) {
            RM.buildAABBTree(mesh);
        }

        // mesh repair
        RM.repair(mesh, 1e-6);
        MESHIO::resetOrientation(mesh);

        // core remesh
        RM.remesh(mesh, mesh);
        std::cout << "remesh end" << std::endl;

        // mesh repair
        RM.repair(mesh);

        // mesh output
        MESHIO::getData(mesh, points_out, triangles_out, surfaceID_out);

        if (b_write_vtk) {
            std::cout << "write file to " << write_vtk_file_name << std::endl;
            MESHIO::writeVTK(write_vtk_file_name, mesh, "surface_id");
        }

    }
    else {
        // split mesh to several small part

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
        RM.repair(mesh, 1e-6);
        MESHIO::resetOrientation(mesh);

        // intesecetion aabb treee
        if (b_useAABBtree) {
            RM.buildAABBTree(mesh);
        }

        // split mesh to several small part
        std::vector<Mesh> split_mesh_list;
        MESHIO::splitDifferentFaces(mesh, split_mesh_list);

        // global parameter
        RM.addSizeFunction(size_function);

        for (int i = 0; i < split_mesh_list.size(); i++) {
            //if (i != 21) {
            //    continue;
            //}
            //std::string before_write_vtk_file_name = write_vtk_file_name;
            //before_write_vtk_file_name += "before_10.vtk";
            //MESHIO::writeVTK(before_write_vtk_file_name, split_mesh_list[i], "surface_id");

            // local parameter

            // core remesh
            std::cout << "remesh" << i << "th part " << std::endl;
            RM.remesh(split_mesh_list[i], split_mesh_list[i]);
            std::cout << "remesh end" << std::endl;

            //if (b_write_vtk) {
            //    std::cout << "write file to " << write_vtk_file_name << std::endl;
            //    std::string new_write_vtk_file_name = write_vtk_file_name;
            //    new_write_vtk_file_name += std::to_string(i) + ".vtk";
            //    MESHIO::writeVTK(new_write_vtk_file_name, split_mesh_list[i], "surface_id");
            //}

        }

        // mesh repair
        MESHIO::addMesh(split_mesh_list, mesh);
        RM.repair(mesh, 1e-4);

        // mesh output
        MESHIO::getData(mesh, points_out, triangles_out, surfaceID_out);

        if (b_write_vtk) {
            std::cout << "write file to " << write_vtk_file_name << std::endl;
            MESHIO::writeVTK(write_vtk_file_name, mesh, "surface_id");
        }
    }

    return 0;
}