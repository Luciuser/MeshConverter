cmake_minimum_required (VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
project(MeshConverter)
set(CMAKE_CXX_FLAGS "-std=c++11" ${CMAKE_CXX_FLAGS})
# options
option(MESHORIENT_BUILD_EXECUTABLE "Build executable." OFF)
option(MESHORIENT_DOWNLOAD "Download extern dependent library." OFF)
if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()
 

if(DOWNLOAD_WITH_CNPM)
    set (GITHUB_REPOSITE "github.com.cnpmjs.org")
elseif(DOWNLOAD_WITH_GITCLONE)
    set (GITHUB_REPOSITE "gitclone.com/github.com")
else()
    set (GITHUB_REPOSITE "github.com")
endif()

if(NOT EXTERNAL_DIR)
    set(EXTERNAL_DIR ${CMAKE_SOURCE_DIR}/extern)
endif()

if(MESHORIENT_DOWNLOAD)
    include(MeshConverterDependencies)
else()
endif()

set(EIGEN_BASE ${EXTERNAL_DIR}/eigen)
set(IGL_BASE ${EXTERNAL_DIR}/libigl)

add_compile_options(-fPIC)
IF(MSVC)
 add_compile_options(-bigobj)
endif()

add_subdirectory(src/MeshOrient)

# remesh
set(SREMESH_SOURCES
    src/meshIO.cpp
    src/meshAlgorithm.cpp
    src/Remesh/API_SRemesh.cpp
    src/Remesh/tiger_sizingfunction.h
    
    src/Halfedge/PolyMesh_Base.cpp
    src/Halfedge/PolyMesh.cpp
    src/Halfedge/PolyMeshIterators.cpp
    src/Halfedge/AABB_Tree.cpp
)

add_library(SRmesh_lib SHARED ${SREMESH_SOURCES})
target_link_libraries(SRmesh_lib MeshOrient)
target_link_libraries(SRmesh_lib MeshOrient)

set(SOURCES
    src/meshIO.cpp
    src/meshAlgorithm.cpp
    src/MeshConverter.cpp
    src/remesh.cpp
    
    src/Halfedge/PolyMesh_Base.cpp
    src/Halfedge/PolyMesh.cpp
    src/Halfedge/PolyMeshIterators.cpp
    src/Halfedge/AABB_Tree.cpp
    
    )

include_directories(./src)
include_directories(./extern/cli11)
include_directories(${EIGEN_BASE})
include_directories(${IGL_BASE}/include)
include_directories(./src/MeshOrient)

add_executable(MeshConverter ${SOURCES})
target_link_libraries(MeshConverter MeshOrient SRmesh_lib)



#add_executable(SRmesh_exe src/Remesh/SRemesh_test)
#target_link_libraries(SRmesh_exe SRmesh_lib)