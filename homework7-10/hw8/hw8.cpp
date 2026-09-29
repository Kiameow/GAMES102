#include <iostream>
#include <easy3d/viewer/viewer.h>
#include <easy3d/core/surface_mesh.h>

using namespace easy3d;

void main() {
    std::cout << "Hello World" << std::endl;

    // do not use temporal object here, the viewer would do delete to the mesh
    SurfaceMesh* mesh = new SurfaceMesh();
    auto v0 = mesh->add_vertex(vec3(0, 0, 0));
    auto v1 = mesh->add_vertex(vec3(1, 0, 0));
    auto v2 = mesh->add_vertex(vec3(0, 1, 0));
    mesh->add_triangle(v0, v1, v2);

    Viewer viewer("Easy3D test");
    viewer.add_model(mesh);
    viewer.run();
}
