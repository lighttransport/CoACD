#include "clip.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <sstream>
#include <stdexcept>

#if COACD_USE_CDT_TRIANGULATION
#include "CDT.h"
#include "CDTUtils.h"
#endif

namespace coacd
{
    namespace
    {
        using vec2d = std::array<double, 2>;
        constexpr double kTriangulationEps = 1e-10;

        bool SamePoint2D(const vec2d &a, const vec2d &b, double eps = kTriangulationEps)
        {
            return fabs(a[0] - b[0]) <= eps && fabs(a[1] - b[1]) <= eps;
        }

        double SignedArea2D(const vector<vec2d> &points, const vector<int> &cycle)
        {
            double area = 0.0;
            for (size_t i = 0; i < cycle.size(); ++i)
            {
                const vec2d &p = points[cycle[i] - 1];
                const vec2d &q = points[cycle[(i + 1) % cycle.size()] - 1];
                area += p[0] * q[1] - q[0] * p[1];
            }
            return area * 0.5;
        }

        double Cross2D(const vec2d &a, const vec2d &b, const vec2d &c)
        {
            return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        }

        bool PointOnSegment2D(const vec2d &p, const vec2d &a, const vec2d &b, double eps = kTriangulationEps)
        {
            if (fabs(Cross2D(a, b, p)) > eps)
                return false;

            return std::min(a[0], b[0]) - eps <= p[0] && p[0] <= std::max(a[0], b[0]) + eps &&
                   std::min(a[1], b[1]) - eps <= p[1] && p[1] <= std::max(a[1], b[1]) + eps;
        }

        int Orientation2D(const vec2d &a, const vec2d &b, const vec2d &c, double eps = kTriangulationEps)
        {
            const double cross = Cross2D(a, b, c);
            if (cross > eps)
                return 1;
            if (cross < -eps)
                return -1;
            return 0;
        }

        bool SegmentsIntersect2D(const vec2d &a0, const vec2d &a1, const vec2d &b0, const vec2d &b1,
                                 bool allow_shared_endpoints = true)
        {
            if (allow_shared_endpoints)
            {
                if (SamePoint2D(a0, b0) || SamePoint2D(a0, b1) || SamePoint2D(a1, b0) || SamePoint2D(a1, b1))
                    return false;
            }

            const int o1 = Orientation2D(a0, a1, b0);
            const int o2 = Orientation2D(a0, a1, b1);
            const int o3 = Orientation2D(b0, b1, a0);
            const int o4 = Orientation2D(b0, b1, a1);

            if (o1 != o2 && o3 != o4)
                return true;

            if (o1 == 0 && PointOnSegment2D(b0, a0, a1))
                return true;
            if (o2 == 0 && PointOnSegment2D(b1, a0, a1))
                return true;
            if (o3 == 0 && PointOnSegment2D(a0, b0, b1))
                return true;
            if (o4 == 0 && PointOnSegment2D(a1, b0, b1))
                return true;

            return false;
        }

        bool PointInPolygon2D(const vec2d &point, const vector<vec2d> &points, const vector<int> &cycle)
        {
            bool inside = false;
            for (size_t i = 0, j = cycle.size() - 1; i < cycle.size(); j = i++)
            {
                const vec2d &pi = points[cycle[i] - 1];
                const vec2d &pj = points[cycle[j] - 1];

                if (PointOnSegment2D(point, pj, pi))
                    return true;

                const bool intersect = ((pi[1] > point[1]) != (pj[1] > point[1])) &&
                                       (point[0] < (pj[0] - pi[0]) * (point[1] - pi[1]) / ((pj[1] - pi[1]) + 1e-30) + pi[0]);
                if (intersect)
                    inside = !inside;
            }
            return inside;
        }

        bool PointInTriangle2D(const vec2d &p, const vec2d &a, const vec2d &b, const vec2d &c)
        {
            const double c0 = Cross2D(a, b, p);
            const double c1 = Cross2D(b, c, p);
            const double c2 = Cross2D(c, a, p);

            const bool all_pos = c0 > kTriangulationEps && c1 > kTriangulationEps &&
                                 c2 > kTriangulationEps;
            const bool all_neg = c0 < -kTriangulationEps && c1 < -kTriangulationEps &&
                                 c2 < -kTriangulationEps;
            return all_pos || all_neg;
        }

        double SquaredDistance2D(const vec2d &a, const vec2d &b)
        {
            const double dx = a[0] - b[0];
            const double dy = a[1] - b[1];
            return dx * dx + dy * dy;
        }

        vector<pair<int, int>> BuildBoundaryEdges(const vector<pair<int, int>> &edges)
        {
            std::set<pair<int, int>> unique_edges;
            for (const auto &edge : edges)
            {
                if (edge.first != edge.second)
                {
                    const pair<int, int> canonical = edge.first < edge.second
                                                         ? edge
                                                         : pair<int, int>(edge.second, edge.first);
                    unique_edges.insert(canonical);
                }
            }

            return vector<pair<int, int>>(unique_edges.begin(), unique_edges.end());
        }

        bool TriangulationDebugEnabled()
        {
            return std::getenv("COACD_TRIANGULATION_DEBUG") != nullptr;
        }

        vector<int> CleanupPolygon(const vector<int> &polygon,
                                   const vector<vec2d> &points,
                                   bool remove_collinear = true);

        std::string CycleKey(const vector<int> &cycle)
        {
            vector<pair<int, int>> edges;
            edges.reserve(cycle.size());
            for (size_t i = 0; i < cycle.size(); ++i)
            {
                int a = cycle[i];
                int b = cycle[(i + 1) % cycle.size()];
                if (a > b)
                    std::swap(a, b);
                edges.push_back({a, b});
            }
            std::sort(edges.begin(), edges.end());

            std::ostringstream key;
            for (const auto &edge : edges)
                key << edge.first << '-' << edge.second << ';';
            return key.str();
        }

        vector<vector<int>> ExtractPlanarBoundaryCycles(const vector<vec2d> &points,
                                                        const vector<pair<int, int>> &edges)
        {
            struct HalfEdge
            {
                int from = -1;
                int to = -1;
                int twin = -1;
                int next = -1;
                double angle = 0.0;
                bool visited = false;
            };

            const vector<pair<int, int>> boundary_edges = BuildBoundaryEdges(edges);
            vector<HalfEdge> halfedges;
            halfedges.reserve(boundary_edges.size() * 2);
            map<int, vector<int>> outgoing;

            for (const auto &edge : boundary_edges)
            {
                const int halfedge_index = static_cast<int>(halfedges.size());
                const vec2d &a = points[edge.first - 1];
                const vec2d &b = points[edge.second - 1];
                halfedges.push_back({edge.first, edge.second, halfedge_index + 1, -1,
                                     std::atan2(b[1] - a[1], b[0] - a[0]), false});
                halfedges.push_back({edge.second, edge.first, halfedge_index, -1,
                                     std::atan2(a[1] - b[1], a[0] - b[0]), false});
                outgoing[edge.first].push_back(halfedge_index);
                outgoing[edge.second].push_back(halfedge_index + 1);
            }

            for (auto &[vertex, incident] : outgoing)
            {
                (void)vertex;
                std::sort(incident.begin(), incident.end(),
                          [&](int lhs, int rhs) { return halfedges[lhs].angle < halfedges[rhs].angle; });
            }

            for (HalfEdge &halfedge : halfedges)
            {
                auto &incident = outgoing[halfedge.to];
                auto it = std::find(incident.begin(), incident.end(), halfedge.twin);
                if (it == incident.end())
                    continue;

                const int twin_pos = static_cast<int>(std::distance(incident.begin(), it));
                const int next_pos = (twin_pos - 1 + static_cast<int>(incident.size())) %
                                     static_cast<int>(incident.size());
                halfedge.next = incident[next_pos];
            }

            vector<vector<int>> cycles;
            std::set<std::string> seen_cycles;
            for (HalfEdge &start_edge : halfedges)
            {
                if (start_edge.visited || start_edge.next < 0)
                    continue;

                vector<int> cycle;
                HalfEdge *edge = &start_edge;
                for (size_t guard = 0; guard <= halfedges.size(); ++guard)
                {
                    if (edge->visited)
                        break;

                    edge->visited = true;
                    cycle.push_back(edge->from);
                    if (edge->next < 0)
                        break;

                    edge = &halfedges[edge->next];
                    if (edge == &start_edge)
                        break;
                }

                cycle = CleanupPolygon(cycle, points, false);
                if (cycle.size() < 3)
                    continue;

                const std::string key = CycleKey(cycle);
                if (seen_cycles.insert(key).second)
                    cycles.push_back(std::move(cycle));
            }

            return cycles;
        }

        bool FindInteriorSample(const vector<vec2d> &points, const vector<int> &cycle, vec2d &sample)
        {
            if (cycle.size() < 3)
                return false;

            const double area = SignedArea2D(points, cycle);
            if (fabs(area) <= kTriangulationEps)
                return false;

            const double sign = area > 0.0 ? 1.0 : -1.0;
            for (size_t i = 0; i < cycle.size(); ++i)
            {
                const vec2d &prev = points[cycle[(i + cycle.size() - 1) % cycle.size()] - 1];
                const vec2d &curr = points[cycle[i] - 1];
                const vec2d &next = points[cycle[(i + 1) % cycle.size()] - 1];
                if (sign * Cross2D(prev, curr, next) > kTriangulationEps)
                {
                    sample = {(prev[0] + curr[0] + next[0]) / 3.0,
                              (prev[1] + curr[1] + next[1]) / 3.0};
                    return true;
                }
            }

            sample = points[cycle[0] - 1];
            return true;
        }

        bool IsBridgeVisible(const vector<vec2d> &points, const vector<int> &outer, const vector<int> &hole,
                             int outer_vertex, int hole_vertex, const vector<vector<int>> &all_holes,
                             const vector<int> &base_outer)
        {
            const vec2d &a = points[outer_vertex - 1];
            const vec2d &b = points[hole_vertex - 1];

            for (size_t i = 0; i < outer.size(); ++i)
            {
                const int u = outer[i];
                const int v = outer[(i + 1) % outer.size()];
                if (u == outer_vertex || v == outer_vertex || u == hole_vertex || v == hole_vertex)
                    continue;
                if (SegmentsIntersect2D(a, b, points[u - 1], points[v - 1]))
                    return false;
            }

            for (const auto &other_hole : all_holes)
            {
                for (size_t i = 0; i < other_hole.size(); ++i)
                {
                    const int u = other_hole[i];
                    const int v = other_hole[(i + 1) % other_hole.size()];
                    if (u == outer_vertex || v == outer_vertex || u == hole_vertex || v == hole_vertex)
                        continue;
                    if (SegmentsIntersect2D(a, b, points[u - 1], points[v - 1]))
                        return false;
                }
            }

            const vec2d midpoint = {(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5};
            if (!PointInPolygon2D(midpoint, points, base_outer))
                return false;
            for (const auto &other_hole : all_holes)
            {
                if (PointInPolygon2D(midpoint, points, other_hole))
                    return false;
            }

            return true;
        }

        bool BridgeHoleIntoPolygon(vector<int> &outer, vector<int> hole, const vector<vec2d> &points,
                                   const vector<vector<int>> &all_holes, const vector<int> &base_outer)
        {
            if (hole.empty())
                return true;

            int hole_start = 0;
            for (int i = 1; i < (int)hole.size(); ++i)
            {
                const vec2d &candidate = points[hole[i] - 1];
                const vec2d &best = points[hole[hole_start] - 1];
                if (candidate[0] > best[0] ||
                    (fabs(candidate[0] - best[0]) <= kTriangulationEps && candidate[1] < best[1]))
                    hole_start = i;
            }

            const int hole_vertex = hole[hole_start];
            int outer_start = -1;
            double best_distance = INF;
            for (int i = 0; i < (int)outer.size(); ++i)
            {
                const int outer_vertex = outer[i];
                if (!IsBridgeVisible(points, outer, hole, outer_vertex, hole_vertex, all_holes, base_outer))
                    continue;

                const double dist = SquaredDistance2D(points[outer_vertex - 1], points[hole_vertex - 1]);
                if (dist < best_distance)
                {
                    best_distance = dist;
                    outer_start = i;
                }
            }

            if (outer_start < 0)
                return false;

            vector<int> merged;
            merged.reserve(outer.size() + hole.size() + 2);
            merged.insert(merged.end(), outer.begin(), outer.begin() + outer_start + 1);
            merged.push_back(hole_vertex);
            for (int i = 1; i < (int)hole.size(); ++i)
                merged.push_back(hole[(hole_start + i) % hole.size()]);
            merged.push_back(hole_vertex);
            merged.push_back(outer[outer_start]);
            merged.insert(merged.end(), outer.begin() + outer_start + 1, outer.end());
            outer.swap(merged);
            return true;
        }

        vector<int> CleanupPolygon(const vector<int> &polygon, const vector<vec2d> &points,
                                   bool remove_collinear)
        {
            vector<int> cleaned;
            cleaned.reserve(polygon.size());
            for (int vertex : polygon)
            {
                if (cleaned.empty() || !SamePoint2D(points[cleaned.back() - 1], points[vertex - 1]))
                    cleaned.push_back(vertex);
            }
            if (cleaned.size() > 1 && SamePoint2D(points[cleaned.front() - 1], points[cleaned.back() - 1]))
                cleaned.pop_back();

            bool changed = remove_collinear;
            while (changed && cleaned.size() >= 3)
            {
                changed = false;
                for (size_t i = 0; i < cleaned.size(); ++i)
                {
                    const size_t prev = (i + cleaned.size() - 1) % cleaned.size();
                    const size_t next = (i + 1) % cleaned.size();
                    const vec2d &a = points[cleaned[prev] - 1];
                    const vec2d &b = points[cleaned[i] - 1];
                    const vec2d &c = points[cleaned[next] - 1];
                    if (fabs(Cross2D(a, b, c)) <= kTriangulationEps &&
                        PointOnSegment2D(b, a, c))
                    {
                        cleaned.erase(cleaned.begin() + i);
                        changed = true;
                        break;
                    }
                }
            }

            return cleaned;
        }

        bool DiagonalIsClear(const vector<int> &polygon, const vector<vec2d> &points, size_t ia, size_t ic)
        {
            const int a_idx = polygon[ia];
            const int c_idx = polygon[ic];
            const vec2d &a = points[a_idx - 1];
            const vec2d &c = points[c_idx - 1];

            for (size_t i = 0; i < polygon.size(); ++i)
            {
                const size_t j = (i + 1) % polygon.size();
                if (i == ia || i == ic || j == ia || j == ic)
                    continue;

                const vec2d &u = points[polygon[i] - 1];
                const vec2d &v = points[polygon[j] - 1];
                if (SegmentsIntersect2D(a, c, u, v))
                    return false;
            }
            return true;
        }

        bool TriangulateSimplePolygon(vector<int> polygon, const vector<vec2d> &points, vector<vec3i> &triangles)
        {
            polygon = CleanupPolygon(polygon, points, false);
            if (polygon.size() < 3)
            {
                if (TriangulationDebugEnabled())
                    std::cerr << "    polygon cleaned below triangle threshold size="
                              << polygon.size() << '\n';
                return false;
            }

            if (TriangulationDebugEnabled())
            {
                std::set<int> unique_vertices(polygon.begin(), polygon.end());
                std::cerr << "    polygon size=" << polygon.size()
                          << " unique_vertices=" << unique_vertices.size()
                          << " area=" << SignedArea2D(points, polygon) << '\n';
            }

            if (SignedArea2D(points, polygon) < 0.0)
                std::reverse(polygon.begin(), polygon.end());

            size_t guard = 0;
            const size_t max_guard = polygon.size() * polygon.size();
            while (polygon.size() > 3 && guard < max_guard)
            {
                bool ear_found = false;
                for (size_t i = 0; i < polygon.size(); ++i)
                {
                    const size_t prev = (i + polygon.size() - 1) % polygon.size();
                    const size_t next = (i + 1) % polygon.size();
                    const int a_idx = polygon[prev];
                    const int b_idx = polygon[i];
                    const int c_idx = polygon[next];
                    const vec2d &a = points[a_idx - 1];
                    const vec2d &b = points[b_idx - 1];
                    const vec2d &c = points[c_idx - 1];

                    if (Cross2D(a, b, c) <= kTriangulationEps)
                        continue;
                    if (!DiagonalIsClear(polygon, points, prev, next))
                        continue;

                    bool contains_point = false;
                    for (size_t j = 0; j < polygon.size(); ++j)
                    {
                        if (j == prev || j == i || j == next)
                            continue;
                        if (PointInTriangle2D(points[polygon[j] - 1], a, b, c))
                        {
                            contains_point = true;
                            break;
                        }
                    }
                    if (contains_point)
                        continue;

                    triangles.push_back({a_idx, b_idx, c_idx});
                    polygon.erase(polygon.begin() + i);
                    ear_found = true;
                    break;
                }

                if (!ear_found)
                {
                    if (TriangulationDebugEnabled())
                    {
                        int intersections = 0;
                        for (size_t a = 0; a < polygon.size(); ++a)
                        {
                            size_t b = (a + 1) % polygon.size();
                            for (size_t c = a + 1; c < polygon.size(); ++c)
                            {
                                size_t d = (c + 1) % polygon.size();
                                if (a == c || a == d || b == c || b == d)
                                    continue;
                                if (a == 0 && d == polygon.size() - 1)
                                    continue;
                                if (SegmentsIntersect2D(points[polygon[a] - 1],
                                                        points[polygon[b] - 1],
                                                        points[polygon[c] - 1],
                                                        points[polygon[d] - 1]))
                                    ++intersections;
                            }
                        }
                        std::cerr << "    no ear found for polygon size=" << polygon.size()
                                  << " self_intersections=" << intersections << '\n';
                    }
                    return false;
                }
                ++guard;
            }

            if (polygon.size() != 3)
                return false;

            triangles.push_back({polygon[0], polygon[1], polygon[2]});
            return true;
        }

        bool TriangulateCyclesWithEarClipping(const vector<vec2d> &points, vector<vector<int>> cycles,
                                              vector<vec3i> &triangles)
        {
            struct RingInfo
            {
                vector<int> indices;
                vec2d sample;
                double area = 0.0;
                int depth = 0;
                bool hole = false;
            };

            vector<RingInfo> rings;
            rings.reserve(cycles.size());
            for (auto &cycle : cycles)
            {
                cycle = CleanupPolygon(cycle, points, false);
                if (cycle.size() < 3)
                    continue;

                RingInfo ring;
                ring.indices = cycle;
                ring.area = SignedArea2D(points, ring.indices);
                if (fabs(ring.area) <= kTriangulationEps)
                    continue;
                if (!FindInteriorSample(points, ring.indices, ring.sample))
                    return false;
                rings.push_back(std::move(ring));
            }

            if (rings.empty())
                return false;

            for (size_t i = 0; i < rings.size(); ++i)
            {
                for (size_t j = 0; j < rings.size(); ++j)
                {
                    if (i == j)
                        continue;
                    if (PointInPolygon2D(rings[i].sample, points, rings[j].indices))
                        rings[i].depth++;
                }
                rings[i].hole = (rings[i].depth % 2) == 1;
                if (TriangulationDebugEnabled())
                {
                    std::cerr << "  ring[" << i << "] size=" << rings[i].indices.size()
                              << " depth=" << rings[i].depth
                              << " hole=" << (rings[i].hole ? 1 : 0) << '\n';
                }
            }

            for (size_t i = 0; i < rings.size(); ++i)
            {
                if (rings[i].hole)
                    continue;

                vector<int> outer = rings[i].indices;
                if (SignedArea2D(points, outer) < 0.0)
                    std::reverse(outer.begin(), outer.end());

                vector<vector<int>> holes;
                for (size_t j = 0; j < rings.size(); ++j)
                {
                    if (!rings[j].hole)
                        continue;
                    if (rings[j].depth != rings[i].depth + 1)
                        continue;
                    if (!PointInPolygon2D(rings[j].sample, points, rings[i].indices))
                        continue;

                    vector<int> hole = rings[j].indices;
                    if (SignedArea2D(points, hole) > 0.0)
                        std::reverse(hole.begin(), hole.end());
                    holes.push_back(std::move(hole));
                }

                std::sort(holes.begin(), holes.end(), [&](const vector<int> &lhs, const vector<int> &rhs)
                {
                    const vec2d &a = points[lhs[0] - 1];
                    const vec2d &b = points[rhs[0] - 1];
                    return a[0] > b[0];
                });

                for (size_t h = 0; h < holes.size(); ++h)
                {
                    vector<vector<int>> remaining_holes;
                    remaining_holes.reserve(holes.size());
                    for (size_t other = 0; other < holes.size(); ++other)
                    {
                        if (other != h)
                            remaining_holes.push_back(holes[other]);
                    }

                    if (!BridgeHoleIntoPolygon(outer, holes[h], points, remaining_holes, rings[i].indices))
                    {
                        if (TriangulationDebugEnabled())
                            std::cerr << "  bridge failed for outer ring " << i << " hole " << h << '\n';
                        return false;
                    }
                }

                if (!TriangulateSimplePolygon(outer, points, triangles))
                {
                    if (TriangulationDebugEnabled())
                        std::cerr << "  simple polygon triangulation failed for ring " << i
                                  << " merged_size=" << outer.size() << '\n';
                    return false;
                }
            }

            return !triangles.empty();
        }
    }

    void PrintEdgeSet(const vector<pair<int, int>> &edges)
    {
        ofstream of("../edge.txt");
        for (int i = 0; i < (int)edges.size(); i++)
        {
            of << i + 1 << ' ' << edges[i].first << ' ' << edges[i].second << endl;
        }
        of.close();
    }

    bool CreatePlaneRotationMatrix(const vector<vec3d> &border, const vector<pair<int, int>> &border_edges, vec3d &T,
                                   double R[3][3], Plane &plane)
    {
        (void)border_edges;
        int idx0 = 0;
        int idx1 = -1;
        int idx2 = -1;
        bool flag = false;

        for (int i = 1; i < (int)border.size(); i++)
        {
            double dist = sqrt(pow(border[idx0][0] - border[i][0], 2) + pow(border[idx0][1] - border[i][1], 2) + pow(border[idx0][2] - border[i][2], 2));
            if (dist > 0.01)
            {
                flag = true;
                idx1 = i;
                break;
            }
        }
        if (!flag)
            return false;
        flag = false;

        for (int i = 2; i < (int)border.size(); i++)
        {
            if (i == idx1)
                continue;
            vec3d p0 = border[idx0];
            vec3d p1 = border[idx1];
            vec3d p2 = border[i];
            vec3d AB, BC;
            AB[0] = p1[0] - p0[0];
            AB[1] = p1[1] - p0[1];
            AB[2] = p1[2] - p0[2];
            BC[0] = p2[0] - p1[0];
            BC[1] = p2[1] - p1[1];
            BC[2] = p2[2] - p1[2];

            double dot_product = AB[0] * BC[0] + AB[1] * BC[1] + AB[2] * BC[2];
            double res = dot_product / (sqrt(pow(AB[0], 2) + pow(AB[1], 2) + pow(AB[2], 2)) * sqrt(pow(BC[0], 2) + pow(BC[1], 2) + pow(BC[2], 2)));
            if (fabs(fabs(res) - 1) > 1e-6 && fabs(res) < INF)
            {
                flag = true;
                idx2 = i;
                break;
            }
        }
        if (!flag)
            return false;

        double t0, t1, t2;
        vec3d p0 = border[idx0], p1 = border[idx1], p2 = border[idx2];
        vec3d normal = CalFaceNormal(p0, p1, p2);

        if (normal[0] * plane.a > 0 || normal[1] * plane.b > 0 || normal[2] * plane.c > 0)
        {
            p0 = border[idx2];
            p1 = border[idx1];
            p2 = border[idx0];
        }
        plane.pFlag = true;
        plane.p0 = p2;
        plane.p1 = p1;
        plane.p2 = p0;

        T = p0;

        double eps = 0.0;
        R[0][0] = (p0[0] - p1[0]) / (sqrt(pow(p0[0] - p1[0], 2) + pow(p0[1] - p1[1], 2) + pow(p0[2] - p1[2], 2)) + eps);
        R[0][1] = (p0[1] - p1[1]) / (sqrt(pow(p0[0] - p1[0], 2) + pow(p0[1] - p1[1], 2) + pow(p0[2] - p1[2], 2)) + eps);
        R[0][2] = (p0[2] - p1[2]) / (sqrt(pow(p0[0] - p1[0], 2) + pow(p0[1] - p1[1], 2) + pow(p0[2] - p1[2], 2)) + eps);

        t0 = (p2[2] - p0[2]) * R[0][1] - (p2[1] - p0[1]) * R[0][2];
        t1 = (p2[0] - p0[0]) * R[0][2] - (p2[2] - p0[2]) * R[0][0];
        t2 = (p2[1] - p0[1]) * R[0][0] - (p2[0] - p0[0]) * R[0][1];
        R[2][0] = t0 / (sqrt(pow(t0, 2) + pow(t1, 2) + pow(t2, 2)) + eps);
        R[2][1] = t1 / (sqrt(pow(t0, 2) + pow(t1, 2) + pow(t2, 2)) + eps);
        R[2][2] = t2 / (sqrt(pow(t0, 2) + pow(t1, 2) + pow(t2, 2)) + eps);

        t0 = R[2][2] * R[0][1] - R[2][1] * R[0][2];
        t1 = R[2][0] * R[0][2] - R[2][2] * R[0][0];
        t2 = R[2][1] * R[0][0] - R[2][0] * R[0][1];
        R[1][0] = t0 / (sqrt(pow(t0, 2) + pow(t1, 2) + pow(t2, 2)) + eps);
        R[1][1] = t1 / (sqrt(pow(t0, 2) + pow(t1, 2) + pow(t2, 2)) + eps);
        R[1][2] = t2 / (sqrt(pow(t0, 2) + pow(t1, 2) + pow(t2, 2)) + eps);

        return true;
    }

    void SimpleCyclesFromEdges(const vector<pair<int, int>> &edges, vector<vector<int>> &simple_cycles)
    {
        (void)edges;
        simple_cycles.clear();
    }

    void FindCycleDirection(const vector<vec3d> &border, const vector<vector<int>> &cycles, const Plane &plane,
                            map<pair<int, int>, bool> &cycles_dir)
    {
        double R[3][3];
        vec3d T;
        Plane plane_copy = plane;
        if (!CreatePlaneRotationMatrix(border, {}, T, R, plane_copy))
            return;

        vector<vec2d> projected;
        projected.reserve(border.size());
        for (const auto &point : border)
        {
            const double x = point[0] - T[0];
            const double y = point[1] - T[1];
            const double z = point[2] - T[2];
            projected.push_back({
                R[0][0] * x + R[0][1] * y + R[0][2] * z,
                R[1][0] * x + R[1][1] * y + R[1][2] * z,
            });
        }

        for (const auto &cycle : cycles)
        {
            const bool ccw = SignedArea2D(projected, cycle) > 0.0;
            for (size_t i = 0; i < cycle.size(); ++i)
            {
                const int a = cycle[i];
                const int b = cycle[(i + 1) % cycle.size()];
                cycles_dir[{a, b}] = ccw;
            }
        }
    }

    short Triangulation(vector<vec3d> &border, const vector<pair<int, int>> &border_edges, vector<vec3i> &border_triangles, Plane &plane)
    {
        double R[3][3];
        vec3d T;

        bool flag = CreatePlaneRotationMatrix(border, border_edges, T, R, plane);
        if (!flag)
            return 1;

        vector<array<double, 2>> points;
        points.reserve(border.size());

        for (int i = 0; i < (int)border.size(); i++)
        {
            double x, y, z, px, py;
            x = border[i][0] - T[0];
            y = border[i][1] - T[1];
            z = border[i][2] - T[2];

            px = R[0][0] * x + R[0][1] * y + R[0][2] * z;
            py = R[1][0] * x + R[1][1] * y + R[1][2] * z;

            points.push_back({px, py});
        }

#if COACD_USE_CDT_TRIANGULATION
        int borderN = (int)points.size();
        CDT::Triangulation<double> cdt(
            CDT::detail::defaults::vertexInsertionOrder,
            CDT::IntersectingConstraintEdges::TryResolve,
            5e-17);

        try
        {
            cdt.insertVertices(
                points.begin(),
                points.end(),
                [](const std::array<double, 2> &p)
                { return p[0]; },
                [](const std::array<double, 2> &p)
                { return p[1]; });
            cdt.insertEdges(
                border_edges.begin(),
                border_edges.end(),
                [](const std::pair<int, int> &p)
                { return (int)p.first - 1; },
                [](const std::pair<int, int> &p)
                { return (int)p.second - 1; });
            cdt.eraseSuperTriangle();
        }
        catch (const std::runtime_error &e)
        {
            (void)e;
            return 2;
        }

        for (size_t i = 0; i < (size_t)cdt.triangles.size(); i++)
        {
            border_triangles.push_back({(int)cdt.triangles[i].vertices[0] + 1,
                                        (int)cdt.triangles[i].vertices[1] + 1,
                                        (int)cdt.triangles[i].vertices[2] + 1});
        }

        for (int i = (int)border.size(); i < borderN; i++)
        {
            double x, y, z;
            CDT::V2d<double> vertex = cdt.vertices[i];
            x = R[0][0] * vertex.x + R[1][0] * vertex.y + T[0];
            y = R[0][1] * vertex.x + R[1][1] * vertex.y + T[1];
            z = R[0][2] * vertex.x + R[1][2] * vertex.y + T[2];
            border.push_back({x, y, z});
        }

        return 0;
#else
        vector<vec2d> merged_points;
        vector<int> merged_to_border;
        vector<int> border_to_merged(points.size(), -1);
        for (size_t i = 0; i < points.size(); ++i)
        {
            int merged_index = -1;
            for (size_t j = 0; j < merged_points.size(); ++j)
            {
                if (SamePoint2D(points[i], merged_points[j]))
                {
                    merged_index = static_cast<int>(j);
                    break;
                }
            }

            if (merged_index < 0)
            {
                merged_index = static_cast<int>(merged_points.size());
                merged_points.push_back(points[i]);
                merged_to_border.push_back(static_cast<int>(i) + 1);
            }

            border_to_merged[i] = merged_index + 1;
        }

        vector<pair<int, int>> merged_edges;
        merged_edges.reserve(border_edges.size());
        for (const auto &edge : border_edges)
        {
            const int a = border_to_merged[edge.first - 1];
            const int b = border_to_merged[edge.second - 1];
            if (a != b)
                merged_edges.push_back({a, b});
        }

        vector<vector<int>> simple_cycles =
            ExtractPlanarBoundaryCycles(merged_points, merged_edges);
        if (TriangulationDebugEnabled())
        {
            std::cerr << "builtin triangulation: border=" << border.size()
                      << " merged_border=" << merged_points.size()
                      << " input_edges=" << border_edges.size()
                      << " merged_edges=" << merged_edges.size()
                      << " cycles=" << simple_cycles.size() << '\n';
            for (size_t i = 0; i < simple_cycles.size(); ++i)
                std::cerr << "  cycle[" << i << "] size=" << simple_cycles[i].size() << '\n';
        }
        if (simple_cycles.empty())
            return 2;

        vector<vec3i> merged_triangles;
        const bool success =
            TriangulateCyclesWithEarClipping(merged_points, simple_cycles, merged_triangles);
        if (!success && TriangulationDebugEnabled())
            std::cerr << "builtin triangulation: ear clipping failed\n";

        if (success)
        {
            border_triangles.reserve(merged_triangles.size());
            for (const vec3i &triangle : merged_triangles)
            {
                border_triangles.push_back({
                    merged_to_border[triangle[0] - 1],
                    merged_to_border[triangle[1] - 1],
                    merged_to_border[triangle[2] - 1],
                });
            }
        }

        return success ? 0 : 2;
#endif
    }

    void RemoveOutlierTriangles(const vector<vec3d> &border, const vector<vec3d> &overlap,
                                const vector<pair<int, int>> &border_edges, const vector<vec3i> &border_triangles,
                                int oriN, map<int, int> &vertex_map, vector<vec3d> &final_border,
                                vector<vec3i> &final_triangles)
    {
        deque<pair<int, int>> BFS_edges(border_edges.begin(), border_edges.end());
        map<pair<int, int>, pair<int, int>> edge_map;
        map<pair<int, int>, bool> border_map;
        map<pair<int, int>, bool> same_edge_map;
        map<int, bool> overlap_map;
        const int v_lenth = (int)border.size();
        const int f_lenth = (int)border_triangles.size();
        bool *add_vertex = new bool[v_lenth]();
        bool *remove_map = new bool[f_lenth]();

        for (int i = 0; i < (int)overlap.size(); i++)
            for (int j = 0; j < (int)border.size(); j++)
                if (SamePointDetect(overlap[i], border[j]))
                    overlap_map[j + 1] = true;

        for (int i = 0; i < (int)border_edges.size(); i++)
        {
            int v0 = border_edges[i].first, v1 = border_edges[i].second;
            same_edge_map[std::pair<int, int>(v0, v1)] = true;
        }
        for (int i = 0; i < (int)border_edges.size(); i++)
        {
            int v0 = border_edges[i].first, v1 = border_edges[i].second;
            if (same_edge_map.find(std::pair<int, int>(v1, v0)) == same_edge_map.end())
            {
                border_map[std::pair<int, int>(v0, v1)] = true;
                border_map[std::pair<int, int>(v1, v0)] = true;
            }
        }

        int borderN = border.size();
        for (int i = 0; i < (int)border_triangles.size(); i++)
        {
            int v0, v1, v2;
            v0 = border_triangles[i][0];
            v1 = border_triangles[i][1];
            v2 = border_triangles[i][2];

            if (!(v0 >= 1 && v0 <= borderN && v1 >= 1 && v1 <= borderN && v2 >= 1 && v2 <= borderN))
                continue;

            pair<int, int> edge01 = std::pair<int, int>(v0, v1), edge10 = std::pair<int, int>(v1, v0);
            pair<int, int> edge12 = std::pair<int, int>(v1, v2), edge21 = std::pair<int, int>(v2, v1);
            pair<int, int> edge20 = std::pair<int, int>(v2, v0), edge02 = std::pair<int, int>(v0, v2);

            if (!(same_edge_map.find(edge10) != same_edge_map.end() && same_edge_map.find(edge01) != same_edge_map.end()))
            {
                if (edge_map.find(edge10) == edge_map.end())
                    edge_map[edge01] = std::pair<int, int>(i, -1);
                else
                    edge_map[edge10] = std::pair<int, int>(edge_map[edge10].first, i);
            }

            if (!(same_edge_map.find(edge12) != same_edge_map.end() && same_edge_map.find(edge21) != same_edge_map.end()))
            {
                if (edge_map.find(edge21) == edge_map.end())
                    edge_map[edge12] = std::pair<int, int>(i, -1);
                else
                    edge_map[edge21] = std::pair<int, int>(edge_map[edge21].first, i);
            }

            if (!(same_edge_map.find(edge02) != same_edge_map.end() && same_edge_map.find(edge20) != same_edge_map.end()))
            {
                if (edge_map.find(edge02) == edge_map.end())
                    edge_map[edge20] = std::pair<int, int>(i, -1);
                else
                    edge_map[edge02] = std::pair<int, int>(edge_map[edge02].first, i);
            }
        }

        int i = 0;
        while (!BFS_edges.empty())
        {
            pair<int, int> item = BFS_edges[0];
            BFS_edges.pop_front();
            int v0 = item.first, v1 = item.second;
            int idx;
            pair<int, int> edge01 = std::pair<int, int>(v0, v1), edge10 = std::pair<int, int>(v1, v0);
            if (i < (int)border_edges.size() && edge_map.find(edge10) != edge_map.end())
            {
                idx = edge_map[edge10].second;
                if (idx != -1)
                    remove_map[idx] = true;
                idx = edge_map[edge10].first;
                if (idx != -1 && !remove_map[idx] && !FaceOverlap(overlap_map, border_triangles[idx]))
                {
                    remove_map[idx] = true;
                    final_triangles.push_back(border_triangles[idx]);
                    for (int k = 0; k < 3; k++)
                        add_vertex[border_triangles[idx][k] - 1] = true;

                    int p0 = border_triangles[idx][0], p1 = border_triangles[idx][1], p2 = border_triangles[idx][2];
                    if (p2 != v0 && p2 != v1)
                    {
                        pair<int, int> pt12 = std::pair<int, int>(p1, p2), pt20 = std::pair<int, int>(p2, p0);
                        if (border_map.find(pt12) == border_map.end())
                            BFS_edges.push_back(pt12);
                        if (border_map.find(pt20) == border_map.end())
                            BFS_edges.push_back(pt20);
                    }
                    else if (p1 != v0 && p1 != v1)
                    {
                        pair<int, int> pt12 = std::pair<int, int>(p1, p2), pt01 = std::pair<int, int>(p0, p1);
                        if (border_map.find(pt12) == border_map.end())
                            BFS_edges.push_back(pt12);
                        if (border_map.find(pt01) == border_map.end())
                            BFS_edges.push_back(pt01);
                    }
                    else if (p0 != v0 && p0 != v1)
                    {
                        pair<int, int> pt01 = std::pair<int, int>(p0, p1), pt20 = std::pair<int, int>(p2, p0);
                        if (border_map.find(pt01) == border_map.end())
                            BFS_edges.push_back(pt01);
                        if (border_map.find(pt20) == border_map.end())
                            BFS_edges.push_back(pt20);
                    }
                }
            }
            else if (i < (int)border_edges.size() && edge_map.find(edge01) != edge_map.end())
            {
                idx = edge_map[edge01].first;
                if (idx != -1)
                    remove_map[idx] = true;
                idx = edge_map[edge01].second;
                if (idx != -1 && !remove_map[idx] && !FaceOverlap(overlap_map, border_triangles[idx]))
                {
                    remove_map[idx] = true;
                    final_triangles.push_back(border_triangles[idx]);
                    for (int k = 0; k < 3; k++)
                        add_vertex[border_triangles[idx][k] - 1] = true;

                    int p0 = border_triangles[idx][0], p1 = border_triangles[idx][1], p2 = border_triangles[idx][2];
                    if (p2 != v0 && p2 != v1)
                    {
                        pair<int, int> pt21 = std::pair<int, int>(p2, p1), pt02 = std::pair<int, int>(p0, p2);
                        if (border_map.find(pt21) == border_map.end())
                            BFS_edges.push_back(pt21);
                        if (border_map.find(pt02) == border_map.end())
                            BFS_edges.push_back(pt02);
                    }
                    else if (p1 != v0 && p1 != v1)
                    {
                        pair<int, int> pt21 = std::pair<int, int>(p2, p1), pt10 = std::pair<int, int>(p1, p0);
                        if (border_map.find(pt21) == border_map.end())
                            BFS_edges.push_back(pt21);
                        if (border_map.find(pt10) == border_map.end())
                            BFS_edges.push_back(pt10);
                    }
                    else if (p0 != v0 && p0 != v1)
                    {
                        pair<int, int> pt10 = std::pair<int, int>(p1, p0), pt02 = std::pair<int, int>(p0, p2);
                        if (border_map.find(pt10) == border_map.end())
                            BFS_edges.push_back(pt10);
                        if (border_map.find(pt02) == border_map.end())
                            BFS_edges.push_back(pt02);
                    }
                }
            }
            else if (i >= (int)border_edges.size() && (edge_map.find(edge01) != edge_map.end() ||
                                                       edge_map.find(edge10) != edge_map.end()))
            {
                for (int j = 0; j < 2; j++)
                {
                    if (j == 0)
                        if (edge_map.find(edge01) != edge_map.end())
                            idx = edge_map[edge01].first;
                        else
                            idx = edge_map[edge10].first;
                    else if (edge_map.find(edge01) != edge_map.end())
                        idx = edge_map[edge01].second;
                    else
                        idx = edge_map[edge10].second;
                    if (idx != -1 && !remove_map[idx])
                    {
                        remove_map[idx] = true;
                        final_triangles.push_back(border_triangles[idx]);
                        for (int k = 0; k < 3; k++)
                            add_vertex[border_triangles[idx][k] - 1] = true;

                        int p0 = border_triangles[idx][0], p1 = border_triangles[idx][1], p2 = border_triangles[idx][2];
                        if (p2 != v0 && p2 != v1)
                        {
                            BFS_edges.push_back(std::pair<int, int>(p1, p2));
                            BFS_edges.push_back(std::pair<int, int>(p2, p0));
                        }
                        else if (p1 != v0 && p1 != v1)
                        {
                            BFS_edges.push_back(std::pair<int, int>(p1, p2));
                            BFS_edges.push_back(std::pair<int, int>(p0, p1));
                        }
                        else if (p0 != v0 && p0 != v1)
                        {
                            BFS_edges.push_back(std::pair<int, int>(p0, p1));
                            BFS_edges.push_back(std::pair<int, int>(p2, p0));
                        }
                    }
                }
            }
            i++;
        }

        int index = 0;
        for (int i = 0; i < (int)border.size(); i++)
        {
            if (i < oriN || add_vertex[i] == true)
            {
                final_border.push_back(border[i]);
                vertex_map[i + 1] = ++index;
            }
        }
        delete[] add_vertex;
        delete[] remove_map;
    }
}
