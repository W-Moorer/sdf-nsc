#include <iostream>
#include <string>
#include <chrono/geometry/ChTriangleMeshConnected.h>
#include <platform/backend/spcc/VDBSDFField.h>

using namespace platform::backend::spcc;
using namespace chrono;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_obj_file>\n";
        return 1;
    }

    std::string obj_path = argv[1];
    std::cout << "Loading mesh from: " << obj_path << std::endl;

    ChTriangleMeshConnected trimesh;
    if (!trimesh.LoadWavefrontMesh(obj_path, true, true)) {
        std::cerr << "Failed to load OBJ file!\n";
        return 1;
    }

    std::cout << "Mesh loaded: " << trimesh.GetNumTriangles() << " triangles.\n";

    VDBSDFField sdf_field;
    VDBSDFField::BuildOptions options;
    options.voxel_size = 0.01;
    std::cout << "Converting to SDF (dX = " << options.voxel_size << ")...\n";
    
    // Convert!
    bool result = sdf_field.BuildFromTriangleMesh(trimesh, options);

    std::cout << "Conversion successful! Ready? " << sdf_field.IsReady() << "\n";

    // Query tests
    ChVector3d p_out(1.0, 1.0, 1.0);  // presumably outside
    ChVector3d p_in(0.0, 0.0, 0.0);   // presumably inside depending on mesh center
    
    double dist_out;
    ChVector3d grad_out;
    sdf_field.QueryPhiGradM(p_out, dist_out, grad_out);
    
    double dist_in;
    ChVector3d grad_in;
    sdf_field.QueryPhiGradM(p_in, dist_in, grad_in);

    std::cout << "Query at (1,1,1): dist = " << dist_out << ", grad = (" << grad_out.x() << ", " << grad_out.y() << ", " << grad_out.z() << ")\n";
    std::cout << "Query at (0,0,0): dist = " << dist_in << ", grad = (" << grad_in.x() << ", " << grad_in.y() << ", " << grad_in.z() << ")\n";

    // Also do a grid check based on bounding box
    ChVector3d bb_min = trimesh.GetBoundingBox().min;
    ChVector3d bb_max = trimesh.GetBoundingBox().max;
    std::cout << "Bounding box: min = " << bb_min.x() << "," << bb_min.y() << "," << bb_min.z() 
              << " max = " << bb_max.x() << "," << bb_max.y() << "," << bb_max.z() << "\n";

    // Find a point actually inside or near the surface by stepping through the volume
    ChVector3d center = (bb_min + bb_max) * 0.5;
    bool found_surface = false;
    for (int ix = 0; ix < 10 && !found_surface; ++ix) {
        for (int iy = 0; iy < 10 && !found_surface; ++iy) {
            for (int iz = 0; iz < 10; ++iz) {
                ChVector3d p(
                    bb_min.x() + (bb_max.x()-bb_min.x()) * (ix/9.0),
                    bb_min.y() + (bb_max.y()-bb_min.y()) * (iy/9.0),
                    bb_min.z() + (bb_max.z()-bb_min.z()) * (iz/9.0)
                );
                double d; ChVector3d g;
                sdf_field.QueryPhiGradM(p, d, g);
                if (d < 0.02) { // Found something interacting with the narrow band
                    std::cout << "Found interesting point at " << p.x() << "," << p.y() << "," << p.z() 
                              << " => dist: " << d << " grad: " << g.x() << "," << g.y() << "," << g.z() << "\n";
                    found_surface = true;
                    break;
                }
            }
        }
    }

    return 0;
}
