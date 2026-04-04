#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <vector>
#include <cmath>
#include <algorithm>

#include <chrono/geometry/ChTriangleMeshConnected.h>

#include <openvdb/openvdb.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/io/File.h>

using namespace chrono;
namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_mesh.obj> [resolution_divider=300]\n";
        std::cerr << "resolution_divider determines adaptive voxel size: voxel_size = max_extent / divider\n";
        return 1;
    }

    std::string input_mesh_path = argv[1];
    float resolution_divider = 280.0f; // approx 0.02mm voxel size for 5.6mm extent
    float half_band = 100.0f; // matches 100 in Chrono

    if (argc >= 3) {
        resolution_divider = std::stof(argv[2]);
    }
    if (argc >= 4) {
        half_band = std::stof(argv[3]);
    }

    openvdb::initialize();

    // 1. Load the mesh
    ChTriangleMeshConnected trimesh;
    if (!trimesh.LoadWavefrontMesh(input_mesh_path, true, true)) {
        std::cerr << "Failed to load input mesh: " << input_mesh_path << "\n";
        return 1;
    }
    
    // 2. Compute adaptive voxel size
    ChVector3d min_pt = trimesh.GetBoundingBox().min;
    ChVector3d max_pt = trimesh.GetBoundingBox().max;
    ChVector3d extent = max_pt - min_pt;
    double max_extent = std::max({extent.x(), extent.y(), extent.z()});
    double voxel_size = max_extent / resolution_divider;
    
    std::cout << "Mesh loaded: " << trimesh.GetNumTriangles() << " triangles.\n";
    std::cout << "Bounding box extent: " << extent.x() << " x " << extent.y() << " x " << extent.z() << "\n";
    std::cout << "Adaptive voxel size (max_extent / " << resolution_divider << "): " << voxel_size << "\n";

    // Prepare point and triangle lists for OpenVDB
    const auto& vertices = trimesh.GetCoordsVertices();
    const auto& triangles = trimesh.GetIndicesVertexes();
    
    std::vector<openvdb::Vec3s> points;
    points.reserve(vertices.size());
    for (const auto& v : vertices) {
        points.emplace_back(static_cast<float>(v.x()), static_cast<float>(v.y()), static_cast<float>(v.z()));
    }
    std::vector<openvdb::Vec3I> tri_faces;
    tri_faces.reserve(triangles.size());
    for (const auto& tri : triangles) {
        tri_faces.emplace_back(tri.x(), tri.y(), tri.z());
    }
    std::vector<openvdb::Vec4I> quad_faces; // empty
    
    auto transform = openvdb::math::Transform::createLinearTransform(voxel_size);
    std::cout << "Using half band width (voxels): " << half_band << "\n";
    std::cout << "Converting mesh to VDB volume...\n";
    openvdb::FloatGrid::Ptr grid = openvdb::tools::meshToSignedDistanceField<openvdb::FloatGrid>(
        *transform, points, tri_faces, quad_faces, half_band, half_band
    );
    grid->setName("master_sdf");
    grid->setGridClass(openvdb::GRID_LEVEL_SET);

    // 3. Prepare output directories
    fs::path in_path(input_mesh_path);
    std::string base_name = in_path.stem().string();
    fs::path base_out_dir = fs::path("data") / "openvdb";
    fs::path out_vdb_dir = base_out_dir / "vdb";
    fs::path out_obj_dir = base_out_dir / "obj";
    
    fs::create_directories(out_vdb_dir);
    fs::create_directories(out_obj_dir);
    
    fs::path out_vdb_path = out_vdb_dir / (base_name + ".vdb");
    fs::path out_obj_path = out_obj_dir / (base_name + ".obj");

    // 4. Save VDB
    std::cout << "Saving VDB file to: " << out_vdb_path.string() << "\n";
    openvdb::io::File vdb_file(out_vdb_path.string());
    openvdb::GridPtrVec grids;
    grids.push_back(grid);
    vdb_file.write(grids);
    vdb_file.close();

    // 5. Volume to Mesh (Zero-isosurface)
    std::cout << "Extracting zero-isosurface to mesh...\n";
    openvdb::tools::VolumeToMesh mesher(0.0 /* isovalue */);
    mesher(*grid);
    
    // 6. Save zero-isosurface OBJ
    std::cout << "Saving visualization OBJ file to: " << out_obj_path.string() << "\n";
    std::ofstream obj_out(out_obj_path.string());
    if (!obj_out.is_open()) {
        std::cerr << "Failed to open output OBJ file for writing!\n";
        return 1;
    }
    // Write vertices
    for (size_t i = 0; i < mesher.pointListSize(); ++i) {
        openvdb::Vec3s pt = mesher.pointList()[i];
        obj_out << "v " << pt.x() << " " << pt.y() << " " << pt.z() << "\n";
    }
    // Write polygons
    for (size_t i = 0; i < mesher.polygonPoolListSize(); ++i) {
        const openvdb::tools::PolygonPool& pool = mesher.polygonPoolList()[i];
        for (size_t j = 0; j < pool.numTriangles(); ++j) {
            openvdb::Vec3I tri = pool.triangle(j);
            // OBJ indices are 1-based
            obj_out << "f " << (tri[0] + 1) << " " << (tri[1] + 1) << " " << (tri[2] + 1) << "\n";
        }
        for (size_t j = 0; j < pool.numQuads(); ++j) {
            openvdb::Vec4I quad = pool.quad(j);
            obj_out << "f " << (quad[0] + 1) << " " << (quad[1] + 1) << " " << (quad[2] + 1) << " " << (quad[3] + 1) << "\n";
        }
    }
    obj_out.close();

    std::cout << "Done! Output files generated successfully.\n";
    return 0;
}