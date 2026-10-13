#define NOMINMAX

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <easy3d/util/initializer.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/viewer/viewer.h>
#include <easy3d/renderer/renderer.h>
#include <easy3d/renderer/drawable_lines.h>
#include <easy3d/renderer/drawable_points.h>
#include <easy3d/algo/delaunay_2d.h>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/core/random.h>

using namespace easy3d;

// ============================================================
// 工具函数（与你之前一致）
// ============================================================
std::vector<vec2> GenerateRandomPoints2d(int n) {
    std::vector<vec2> points;
    points.reserve(n);
    for (int i = 0; i < n; ++i)
        points.emplace_back(random_float(), random_float());
    return points;
}

static vec2 intersect(const vec2& a, const vec2& b, int edge, float val) {
    if (edge == 0 || edge == 1) {
        float t = (val - a.x) / (b.x - a.x);
        return vec2(val, a.y + t * (b.y - a.y));
    } else {
        float t = (val - a.y) / (b.y - a.y);
        return vec2(a.x + t * (b.x - a.x), val);
    }
}

std::vector<vec2> clipEdge(const std::vector<vec2>& poly, int edge, float val) {
    std::vector<vec2> out;
    for (size_t i = 0; i < poly.size(); ++i) {
        vec2 cur = poly[i];
        vec2 prev = poly[(i + poly.size() - 1) % poly.size()];
        auto inside = [&](vec2 p) {
            if (edge == 0) return p.x >= val;
            if (edge == 1) return p.x <= val;
            if (edge == 2) return p.y >= val;
            return p.y <= val;
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
    float cross_2d = ab.x * ac.y - ab.y * ac.x;

    // 面积过小 → 退化，返回三角形重心代替（比 NaN 安全）
    if (std::abs(cross_2d) < 1e-12f) {
        return vec2((a.x + b.x + c.x) / 3.0f,
                    (a.y + b.y + c.y) / 3.0f);
    }

    vec2 ab_ortho(-ab.y, ab.x);
    vec2 ac_ortho(-ac.y, ac.x);
    vec2 O = (ab_ortho * ac.length2() - ac_ortho * ab.length2())
             / (cross_2d * 2.0f);
    return O + a;
}

namespace {

double signedArea(const std::vector<vec2>& poly) {
    const size_t n = poly.size();
    if (n < 3) return 0.0;
    double A = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const vec2& p = poly[i];
        const vec2& q = poly[(i + 1) % n];
        A += (double)p.x * q.y - (double)q.x * p.y;
    }
    return 0.5 * A;
}

vec2 polygonCentroid(const std::vector<vec2>& poly) {
    const size_t n = poly.size();
    if (n == 0) return vec2(0.0f, 0.0f);
    if (n < 3) {
        vec2 s(0.0f, 0.0f);
        for (auto& p : poly) s += p;
        return vec2(s.x / (float)n, s.y / (float)n);
    }
    double A = 0.0, cx = 0.0, cy = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const vec2& p = poly[i];
        const vec2& q = poly[(i + 1) % n];
        const double cr = (double)p.x * q.y - (double)q.x * p.y;
        A += cr;
        cx += ((double)p.x + q.x) * cr;
        cy += ((double)p.y + q.y) * cr;
    }
    A *= 0.5;
    if (std::abs(A) < 1e-12) {
        vec2 s(0.0f, 0.0f);
        for (auto& p : poly) s += p;
        return vec2(s.x / (float)n, s.y / (float)n);
    }
    return vec2((float)(cx / (6.0 * A)), (float)(cy / (6.0 * A)));
}

std::vector<vec2> extractPolygon(const std::vector<vec3>& line_pts) {
    std::vector<vec2> poly;
    poly.reserve(line_pts.size());
    constexpr float eps2 = 1e-14f;
    for (const auto& p : line_pts) {
        vec2 q(p.x, p.y);
        if (!poly.empty()) {
            const vec2& last = poly.back();
            float dx = q.x - last.x, dy = q.y - last.y;
            if (dx * dx + dy * dy < eps2) continue;
        }
        poly.push_back(q);
    }
    if (poly.size() > 1) {
        const vec2& a = poly.front();
        const vec2& b = poly.back();
        float dx = a.x - b.x, dy = a.y - b.y;
        if (dx * dx + dy * dy < eps2) poly.pop_back();
    }
    return poly;
}

void makeCCW(std::vector<vec2>& poly) {
    if (signedArea(poly) < 0.0) std::reverse(poly.begin(), poly.end());
}

double sceneDiagonal(const std::vector<vec2>& pts) {
    if (pts.empty()) return 1.0;
    double minx = pts[0].x, maxx = pts[0].x;
    double miny = pts[0].y, maxy = pts[0].y;
    for (auto& p : pts) {
        minx = std::min(minx, (double)p.x); maxx = std::max(maxx, (double)p.x);
        miny = std::min(miny, (double)p.y); maxy = std::max(maxy, (double)p.y);
    }
    return std::hypot(maxx - minx, maxy - miny);
}

} // namespace

// ============================================================
// CVTViewer：交互式 CVT
// ============================================================
class CVTViewer : public Viewer {
public:
    explicit CVTViewer(const std::string& title) : Viewer(title) {}

    void run_cvt() {
        points     = GenerateRandomPoints2d(N);
        old_points = points;
        iter_count = 0;
        cached_voronoi_lines.assign(N, {});

        rebuild_display();

        std::cout << "\n==== CVT Viewer ====\n"
                  << "  'n' : 执行 20 次 Lloyd 迭代\n"
                  << "  'r' : 重置为随机点\n"
                  << "  'q' / ESC : 退出\n"
                  << "=====================\n\n";

        run();   // 只调用一次
    }

protected:
    bool key_press_event(int key, int modifiers) override {
        if (key == 'N') {
            bool converged = false;
            int done = 0;
            for (int i = 0; i < 20; ++i) {
                converged = iterate_once();
                ++iter_count;
                ++done;
                if (converged) break;
            }
            rebuild_display();

            std::cout << "[CVT] this batch: " << done
                      << "  total: " << iter_count;
            if (converged) std::cout << " (converged)";
            std::cout << "\n";
            return true;
        }
        if (key == 'R') {
            points     = GenerateRandomPoints2d(N);
            old_points = points;
            iter_count = 0;
            cached_voronoi_lines.assign(N, {});
            rebuild_display();
            std::cout << "[CVT] 重置\n";
            return true;
        }
        return Viewer::key_press_event(key, modifiers);
    }

private:
    static constexpr int N = 50;

    std::vector<vec2> points;
    std::vector<vec2> old_points;
    int iter_count = 0;

    SurfaceMesh* mesh  = nullptr;
    PointCloud*  cloud = nullptr;

    // 上一次迭代算出的 Voronoi cell 线段（用于显示）
    std::vector<std::vector<vec3>> cached_voronoi_lines;

    // ---------- 单次 Lloyd 迭代 ----------
    // 返回 true 表示达到收敛阈值
    bool iterate_once() {
        old_points = points;

        // --- Delaunay ---
        Delaunay2 delaunay;
        delaunay.set_vertices(points);
        const int* tri_indices = delaunay.tri_to_v();
        const unsigned int num_triangles = delaunay.nb_triangles();

        // 用局部 mesh 做拓扑计算，避免污染显示用 mesh
        SurfaceMesh tmp;
        std::vector<SurfaceMesh::Vertex> mvs;
        mvs.reserve(points.size());
        for (auto& p : points)
            mvs.push_back(tmp.add_vertex(vec3(p.x, p.y, 0.0f)));

        std::vector<vec2> tri_cc(num_triangles);
        for (unsigned int i = 0; i < num_triangles; ++i) {
            int v0 = tri_indices[i * 3 + 0];
            int v1 = tri_indices[i * 3 + 1];
            int v2 = tri_indices[i * 3 + 2];
            tmp.add_triangle(mvs[v0], mvs[v1], mvs[v2]);
            tri_cc[i] = circumcenter(points[v0], points[v1], points[v2]);
        }

        // --- 网格重心 ---
        vec2 mesh_centroid(0.0f, 0.0f);
        for (auto v : tmp.vertices()) {
            vec3 p = tmp.position(v);
            mesh_centroid.x += p.x;
            mesh_centroid.y += p.y;
        }
        if (tmp.n_vertices() > 0)
            mesh_centroid /= (float)tmp.n_vertices();

        // --- 每个顶点的 Voronoi cell ---
        std::vector<std::vector<vec3>> voronoiLines(points.size());
        const double diag = sceneDiagonal(points);
        const float  L    = (float)(4.0 * (std::max)(diag, 1e-6));

        for (auto v : tmp.vertices()) {
            vec3 vp = tmp.position(v);
            vec2 v2(vp.x, vp.y);

            std::vector<SurfaceMesh::Face> incident_faces;
            std::vector<std::pair<int, vec2>> boundary_infos;
            bool isBoundary = false;

            for (auto h : tmp.halfedges(v)) {
                auto f = tmp.face(h);
                if (f.is_valid()) {
                    incident_faces.push_back(f);
                } else {
                    isBoundary = true;
                    auto fopp = tmp.face(tmp.opposite(h));
                    if (fopp.is_valid()) {
                        auto w = tmp.target(h);
                        vec3 wp = tmp.position(w);
                        boundary_infos.emplace_back((int)fopp.idx(),
                                                    vec2(wp.x, wp.y));
                    }
                }
            }
            if (incident_faces.empty()) continue;

            std::sort(incident_faces.begin(), incident_faces.end(),
                [&](const SurfaceMesh::Face& fa, const SurfaceMesh::Face& fb) {
                    vec2 ca = tri_cc[fa.idx()];
                    vec2 cb = tri_cc[fb.idx()];
                    return std::atan2(ca.y - v2.y, ca.x - v2.x)
                         < std::atan2(cb.y - v2.y, cb.x - v2.x);
                });

            std::vector<vec2> poly;

            if (!isBoundary || boundary_infos.size() < 2) {
                poly.reserve(incident_faces.size());
                for (auto f : incident_faces)
                    poly.push_back(tri_cc[f.idx()]);
            } else {
                const vec2 cc_first = tri_cc[incident_faces.front().idx()];
                const vec2 cc_last  = tri_cc[incident_faces.back().idx()];

                auto find_w = [&](const SurfaceMesh::Face& f) -> vec2 {
                    for (auto& bi : boundary_infos)
                        if (bi.first == (int)f.idx()) return bi.second;
                    return vec2(0.0f, 0.0f);
                };

                auto aux_point = [&](const vec2& cc, const vec2& w) -> vec2 {
                    vec2 e = w - v2;
                    vec2 d(-e.y, e.x);
                    float len = std::sqrt(d.x * d.x + d.y * d.y);
                    if (len > 1e-12f) { d.x /= len; d.y /= len; }
                    vec2 u = mesh_centroid - cc;
                    if (d.x * u.x + d.y * u.y > 0) { d.x = -d.x; d.y = -d.y; }
                    return vec2(cc.x + L * d.x, cc.y + L * d.y);
                };

                const vec2 p_first = aux_point(cc_first, find_w(incident_faces.front()));
                const vec2 p_last  = aux_point(cc_last,  find_w(incident_faces.back()));

                poly.reserve(incident_faces.size() + 2);
                poly.push_back(p_first);
                for (auto f : incident_faces)
                    poly.push_back(tri_cc[f.idx()]);
                poly.push_back(p_last);
            }

            makeCCW(poly);
            std::vector<vec2> clipped = clipToRect(poly, 0.0f, 0.0f, 1.0f, 1.0f);
            if (clipped.size() < 3) continue;
            makeCCW(clipped);

            auto& lines = voronoiLines[v.idx()];
            lines.reserve(clipped.size() * 2);
            for (size_t i = 0; i < clipped.size(); ++i) {
                const vec2& a = clipped[i];
                const vec2& b = clipped[(i + 1) % clipped.size()];
                lines.emplace_back(a.x, a.y, 0.0f);
                lines.emplace_back(b.x, b.y, 0.0f);
            }
        }

        // --- 更新 points 为面积质心 ---
        for (size_t i = 0; i < points.size(); ++i) {
            std::vector<vec2> poly = extractPolygon(voronoiLines[i]);
            if (poly.size() < 3) continue;   // 退化，保持原位
            makeCCW(poly);
            points[i] = polygonCentroid(poly);
        }

        // --- 收敛判断 ---
        double maxDisp = 0.0, sumDisp = 0.0;
        for (size_t i = 0; i < points.size(); ++i) {
            double dx = points[i].x - old_points[i].x;
            double dy = points[i].y - old_points[i].y;
            double d  = std::hypot(dx, dy);
            maxDisp = (std::max)(maxDisp, d);
            sumDisp += d;
        }
        double meanDisp = sumDisp / (std::max<size_t>)(points.size(), (size_t)1);
        double tol      = 1e-5 * (std::max)(diag, 1e-6);

        std::cout << "    maxDisp="  << maxDisp
                  << "  meanDisp=" << meanDisp
                  << "  tol="      << tol << "\n";

        // 缓存给显示用
        cached_voronoi_lines = std::move(voronoiLines);

        return maxDisp < tol && meanDisp < 0.5 * tol;
    }

    // ---------- 用当前 points 重建显示 ----------
    void rebuild_display() {
        // 删旧模型（第二个参数 true 表示一并释放内存）
        if (mesh)  { delete_model(mesh); mesh  = nullptr; }
        if (cloud) { delete_model(cloud); cloud = nullptr; }

        // 点云
        cloud = new PointCloud;
        for (auto& p : points)
            cloud->add_vertex(vec3(p.x, p.y, 0.0f));
        add_model(cloud);

        // 三角网格
        mesh = new SurfaceMesh;
        Delaunay2 delaunay;
        delaunay.set_vertices(points);
        const int* tri_indices = delaunay.tri_to_v();
        const unsigned int num_tri = delaunay.nb_triangles();

        std::vector<SurfaceMesh::Vertex> mvs;
        mvs.reserve(points.size());
        for (auto& p : points)
            mvs.push_back(mesh->add_vertex(vec3(p.x, p.y, 0.0f)));
        for (unsigned int i = 0; i < num_tri; ++i) {
            int v0 = tri_indices[i * 3 + 0];
            int v1 = tri_indices[i * 3 + 1];
            int v2 = tri_indices[i * 3 + 2];
            mesh->add_triangle(mvs[v0], mvs[v1], mvs[v2]);
        }
        add_model(mesh);

        // 打开边显示 + 加 Voronoi 线条
        if (auto* mr = mesh->renderer()) {
            if (auto* edges = mr->get_lines_drawable("edges"))
                edges->set_visible(true);

            if (auto* vd = mr->add_lines_drawable("voronoi")) {
                std::vector<vec3> all;
                for (auto& cell : cached_voronoi_lines)
                    all.insert(all.end(), cell.begin(), cell.end());
                vd->update_vertex_buffer(all);
                vd->set_uniform_coloring(vec4(1.0f, 0.0f, 0.0f, 1.0f));
                vd->set_line_width(2.0f);
            }
        }

        update();   // 触发一次重绘
    }
};

// ============================================================
// main
// ============================================================
int main() {
    initialize();
    CVTViewer viewer("hw8 CVT");
    viewer.run_cvt();
    return 0;
}