#include <iostream>
#include <easy3d/util/initializer.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/viewer/viewer.h>
#include <easy3d/renderer/renderer.h>
#include <easy3d/renderer/drawable_lines.h>
#include <easy3d/renderer/drawable_points.h>
#include <easy3d/algo/delaunay_2d.h>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/core/random.h>
#include <vector>
#include <cstdlib>

using namespace easy3d;

std::vector<vec2> GenerateRandomPoints2d(int n, float min_val = 0.0f, float max_val = 1.0f) {
    std::vector<vec2> points;
    points.reserve(n);

    for (int i = 0; i < n; i++) {
        points.emplace_back(random_float(), random_float());
    }

    return points;
} 

static vec2 intersect(const vec2& a, const vec2& b, int edge, float val) {
    // edge: 0=left(x=val), 1=right(x=val), 2=bottom(y=val), 3=top(y=val)
    if (edge == 0 || edge == 1) {
        // 与竖直线 x = val 求交
        float t = (val - a.x) / (b.x - a.x);
        return vec2(val, a.y + t * (b.y - a.y));
    } else {
        // 与水平线 y = val 求交
        float t = (val - a.y) / (b.y - a.y);
        return vec2(a.x + t * (b.x - a.x), val);
    }
}

std::vector<vec2> clipEdge(const std::vector<vec2>& poly, int edge, float val) {
    std::vector<vec2> out;
    for (size_t i = 0; i < poly.size(); i++) {
        vec2 cur = poly[i];
        vec2 prev = poly[(i + poly.size() - 1) % poly.size()];
        auto inside = [&](vec2 p) {
            if (edge == 0) return p.x >= val; // left
            if (edge == 1) return p.x <= val; // right
            if (edge == 2) return p.y >= val; // bottom
            return p.y <= val;                // top
        };
        bool curIn = inside(cur), prevIn = inside(prev);
        if (curIn) {
            if (!prevIn) out.push_back(intersect(prev, cur, edge, val));
            out.push_back(cur);
        } else if (prevIn) {
            out.push_back(intersect(prev, cur, edge, val));
        }
    }
    return out;
}

std::vector<vec2> clipToRect(std::vector<vec2> poly, float xmin, float ymin, float xmax, float ymax) {
    poly = clipEdge(poly, 0, xmin);
    poly = clipEdge(poly, 1, xmax);
    poly = clipEdge(poly, 2, ymin);
    poly = clipEdge(poly, 3, ymax);
    return poly;
}

vec2 circumcenter(const vec2& a, const vec2& b, const vec2& c) {
    vec2 ab = b - a;
    vec2 ac = c - a;

    vec2 ab_ortho(-ab.y, ab.x);
    vec2 ac_ortho(-ac.y, ac.x);

    float cross_2d = ab.x * ac.y - ab.y * ac.x;

    vec2 O = (ab_ortho * ac.length2() - ac_ortho * ab.length2())
             / (cross_2d * 2.0f);

    return O + a;
}

float clamp(float val, float min_val, float max_val) {
    if (val < min_val) return min_val;
    else if (val > max_val) return max_val;

    return val;
}

int main() {
    // 0.initialize easy3d
    initialize();
    
    // 1.random points generation
    srand(42);

    int N = 50;
    std::vector<vec2> points = GenerateRandomPoints2d(N);

    auto cloud = new PointCloud;
    std::vector<vec3> points3d;

    for (const auto& p : points) {
        cloud->add_vertex(vec3(p.x, p.y, 0.0f));
    }

    // 2.Delaunay triangulation
    Delaunay2 delaunay;
    delaunay.set_vertices(points);

    const int* tri_indices = delaunay.tri_to_v();
    unsigned int num_triangles = delaunay.nb_triangles();

    auto mesh = new SurfaceMesh;
    std::vector<SurfaceMesh::Vertex> mesh_vertices;
    for (const auto& p : points) {
        mesh_vertices.push_back(mesh->add_vertex(vec3(p.x, p.y, 0.0f)));
    }

    std::vector<vec2> tri_circumcenter(num_triangles);

    for (unsigned int i = 0; i < num_triangles; i++) {
        int v0 = tri_indices[i * 3 + 0];
        int v1 = tri_indices[i * 3 + 1];
        int v2 = tri_indices[i * 3 + 2];
        mesh->add_triangle(mesh_vertices[v0], mesh_vertices[v1], mesh_vertices[v2]);

        tri_circumcenter[i] = circumcenter(points[v0], points[v1], points[v2]);
    }

    // 3.iterate through vertices and get incident faces and their circumcenters
    // 先算一次网格重心，用于判断外向方向
    vec2 mesh_centroid(0.0f, 0.0f);
    for (auto v : mesh->vertices()) {
        vec3 p = mesh->position(v);
        mesh_centroid.x += p.x;
        mesh_centroid.y += p.y;
    }
    mesh_centroid /= (float)mesh->n_vertices();

    std::vector<std::vector<vec3>> voronoiVertices(points.size());

    for (auto v : mesh->vertices()) {
        vec3 vp = mesh->position(v);
        vec2 v2(vp.x, vp.y);

        std::vector<SurfaceMesh::Face> incident_faces;
        // 每条凸包边记录：(相邻三角形外心 cc, 凸包边另一端点 w)
        std::vector<std::pair<vec2, vec2>> boundary_infos;
        bool isBoundary = false;

        for (auto h : mesh->halfedges(v)) {
            auto f = mesh->face(h);
            if (f.is_valid()) {
                incident_faces.push_back(f);
            } else {
                isBoundary = true;
                // 对向半边所在的面 —— 就是与这条凸包边相邻的那个三角形
                auto fopp = mesh->face(mesh->opposite(h));
                if (fopp.is_valid()) {
                    vec2 cc = tri_circumcenter[fopp.idx()];
                    auto w = mesh->target(h);
                    vec3 wp = mesh->position(w);
                    boundary_infos.emplace_back(cc, vec2(wp.x, wp.y));
                }
            }
        }
        if (incident_faces.empty()) continue;

        // 按外心相对 v 的角度排序
        std::sort(incident_faces.begin(), incident_faces.end(),
            [&](const SurfaceMesh::Face& fa, const SurfaceMesh::Face& fb) {
                vec2 ca = tri_circumcenter[fa.idx()];
                vec2 cb = tri_circumcenter[fb.idx()];
                return std::atan2(ca.y - v2.y, ca.x - v2.x)
                    < std::atan2(cb.y - v2.y, cb.x - v2.x);
            });

        std::vector<vec2> poly;

        if (!isBoundary || boundary_infos.size() < 2) {
            // 内部点：直接连成一圈
            for (auto f : incident_faces)
                poly.push_back(tri_circumcenter[f.idx()]);
        } else {
            // 边界点：在缺口两端各插入一个远端辅助点
            vec2 cc_first = tri_circumcenter[incident_faces.front().idx()];
            vec2 cc_last  = tri_circumcenter[incident_faces.back().idx()];

            // 根据外心匹配出对应的凸包边另一端 w
            auto find_w = [&](const vec2& cc) -> vec2 {
                for (auto& bi : boundary_infos) {
                    if (std::abs(bi.first.x - cc.x) < 1e-9f &&
                        std::abs(bi.first.y - cc.y) < 1e-9f)
                        return bi.second;
                }
                return cc;
            };

            auto aux_point = [&](const vec2& cc, const vec2& w) -> vec2 {
                // 凸包边方向
                vec2 e = w - v2;
                // 垂直平分线方向
                vec2 d(-e.y, e.x);
                float len = std::sqrt(d.x * d.x + d.y * d.y);
                if (len > 1e-12f) { d.x /= len; d.y /= len; }
                // 让 d 指向"远离网格重心"的一侧，即朝外
                vec2 u = mesh_centroid - cc;
                if (d.x * u.x + d.y * u.y > 0) { d.x = -d.x; d.y = -d.y; }
                const float L = 10.0f;   // 足够远，超过 [0,1]^2 的范围即可
                return vec2(cc.x + L * d.x, cc.y + L * d.y);
            };

            vec2 p_first = aux_point(cc_first, find_w(cc_first));
            vec2 p_last  = aux_point(cc_last,  find_w(cc_last));

            // 顺序：p_first → cc_first → ... → cc_last → p_last
            poly.push_back(p_first);
            for (auto f : incident_faces)
                poly.push_back(tri_circumcenter[f.idx()]);
            poly.push_back(p_last);
        }

        // 统一走 Sutherland–Hodgman 裁剪到 [0,1]^2
        std::vector<vec2> clipped = clipToRect(poly, 0.0f, 0.0f, 1.0f, 1.0f);

        // 输出线段
        for (size_t i = 0; i < clipped.size(); ++i) {
            vec2 a = clipped[i];
            vec2 b = clipped[(i + 1) % clipped.size()];
            voronoiVertices[v.idx()].push_back(vec3(a.x, a.y, 0.0f));
            voronoiVertices[v.idx()].push_back(vec3(b.x, b.y, 0.0f));
        }
    }

    Viewer viewer("hw8");
    viewer.add_model(cloud);
    viewer.add_model(mesh);

    auto mesh_renderer = mesh->renderer();
    auto edges = mesh_renderer->get_lines_drawable("edges");
    if (edges) {
        edges->set_visible(true);
    }

    auto voronoi_drawable = mesh_renderer->add_lines_drawable("voronoi");
    if (voronoi_drawable) {
        std::vector<vec3> allVoronoiLines;
        for (const auto& cell : voronoiVertices) {
            allVoronoiLines.insert(allVoronoiLines.end(), cell.begin(), cell.end());
        }

        voronoi_drawable->update_vertex_buffer(allVoronoiLines);
        voronoi_drawable->set_uniform_coloring(vec4(1.0f, 0.0f, 0.0f, 1.0f));
        voronoi_drawable->set_line_width(2.0f);
    }

    viewer.run();

    return 0;
}
