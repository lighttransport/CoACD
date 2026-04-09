#include "clip.h"

namespace coacd
{
    bool Clip(const Model &mesh, Model &pos, Model &neg, Plane &plane, double &cut_area, bool foo)
    {
        (void)foo;
        Model t = mesh;
        vector<vec3d> border;
        vector<vec3d> overlap;
        vector<vec3i> border_triangles, final_triangles;
        vector<pair<int, int>> border_edges;
        map<int, int> border_map;
        vector<vec3d> final_border;

        const int N = (int)mesh.points.size();
        int idx = 0;
        vector<bool> pos_map(N);
        vector<bool> neg_map(N);

        map<pair<int, int>, int> edge_map;
        map<int, int> vertex_map;

        for (int i = 0; i < (int)mesh.triangles.size(); i++)
        {
            int id0, id1, id2;
            id0 = mesh.triangles[i][0];
            id1 = mesh.triangles[i][1];
            id2 = mesh.triangles[i][2];
            vec3d p0, p1, p2;
            p0 = mesh.points[id0];
            p1 = mesh.points[id1];
            p2 = mesh.points[id2];
            short s0 = plane.Side(p0), s1 = plane.Side(p1), s2 = plane.Side(p2);
            short sum = s0 + s1 + s2;
            if (s0 == 0 && s1 == 0 && s2 == 0)
            {
                s0 = s1 = s2 = plane.CutSide(p0, p1, p2, plane);
                sum = s0 + s1 + s2;
                overlap.push_back(p0);
                overlap.push_back(p1);
                overlap.push_back(p2);
            }

            if (sum == 3 || sum == 2 || (sum == 1 && ((s0 == 1 && s1 == 0 && s2 == 0) || (s0 == 0 && s1 == 1 && s2 == 0) || (s0 == 0 && s1 == 0 && s2 == 1))))
            {
                pos_map[id0] = true;
                pos_map[id1] = true;
                pos_map[id2] = true;
                pos.triangles.push_back(mesh.triangles[i]);
                if (sum == 1)
                {
                    if (s0 == 1 && s1 == 0 && s2 == 0)
                    {
                        addPoint(vertex_map, border, p1, id1, idx);
                        addPoint(vertex_map, border, p2, id2, idx);
                        if (vertex_map[id1] != vertex_map[id2])
                            border_edges.push_back(std::pair<int, int>(vertex_map[id1] + 1, vertex_map[id2] + 1));
                    }
                    else if (s0 == 0 && s1 == 1 && s2 == 0)
                    {
                        addPoint(vertex_map, border, p2, id2, idx);
                        addPoint(vertex_map, border, p0, id0, idx);
                        if (vertex_map[id2] != vertex_map[id0])
                            border_edges.push_back(std::pair<int, int>(vertex_map[id2] + 1, vertex_map[id0] + 1));
                    }
                    else if (s0 == 0 && s1 == 0 && s2 == 1)
                    {
                        addPoint(vertex_map, border, p0, id0, idx);
                        addPoint(vertex_map, border, p1, id1, idx);
                        if (vertex_map[id0] != vertex_map[id1])
                            border_edges.push_back(std::pair<int, int>(vertex_map[id0] + 1, vertex_map[id1] + 1));
                    }
                }
            }
            else if (sum == -3 || sum == -2 || (sum == -1 && ((s0 == -1 && s1 == 0 && s2 == 0) || (s0 == 0 && s1 == -1 && s2 == 0) || (s0 == 0 && s1 == 0 && s2 == -1))))
            {
                neg_map[id0] = true;
                neg_map[id1] = true;
                neg_map[id2] = true;
                neg.triangles.push_back(mesh.triangles[i]);
                if (sum == -1)
                {
                    if (s0 == -1 && s1 == 0 && s2 == 0)
                    {
                        addPoint(vertex_map, border, p2, id2, idx);
                        addPoint(vertex_map, border, p1, id1, idx);
                        if (vertex_map[id2] != vertex_map[id1])
                            border_edges.push_back(std::pair<int, int>(vertex_map[id2] + 1, vertex_map[id1] + 1));
                    }
                    else if (s0 == 0 && s1 == -1 && s2 == 0)
                    {
                        addPoint(vertex_map, border, p0, id0, idx);
                        addPoint(vertex_map, border, p2, id2, idx);
                        if (vertex_map[id0] != vertex_map[id2])
                            border_edges.push_back(std::pair<int, int>(vertex_map[id0] + 1, vertex_map[id2] + 1));
                    }
                    else if (s0 == 0 && s1 == 0 && s2 == -1)
                    {
                        addPoint(vertex_map, border, p1, id1, idx);
                        addPoint(vertex_map, border, p0, id0, idx);
                        if (vertex_map[id1] != vertex_map[id0])
                            border_edges.push_back(std::pair<int, int>(vertex_map[id1] + 1, vertex_map[id0] + 1));
                    }
                }
            }
            else
            {
                bool f0, f1, f2;
                vec3d pi0, pi1, pi2;
                f0 = plane.IntersectSegment(p0, p1, pi0);
                f1 = plane.IntersectSegment(p1, p2, pi1);
                f2 = plane.IntersectSegment(p2, p0, pi2);

                if (f0 && f1 && !f2)
                {
                    addEdgePoint(edge_map, border, pi0, id0, id1, idx);
                    addEdgePoint(edge_map, border, pi1, id1, id2, idx);

                    int f0_idx = edge_map[std::pair<int, int>(id0, id1)];
                    int f1_idx = edge_map[std::pair<int, int>(id1, id2)];
                    if (s1 == 1)
                    {
                        if (f1_idx != f0_idx)
                        {
                            border_edges.push_back(std::pair<int, int>(f1_idx + 1, f0_idx + 1));
                            pos_map[id1] = true;
                            neg_map[id0] = true;
                            neg_map[id2] = true;
                            pos.triangles.push_back({id1, -1 * f1_idx - 1, -1 * f0_idx - 1});
                            neg.triangles.push_back({id0, -1 * f0_idx - 1, -1 * f1_idx - 1});
                            neg.triangles.push_back({-1 * f1_idx - 1, id2, id0});
                        }
                        else
                        {
                            neg_map[id0] = true;
                            neg_map[id2] = true;
                            neg.triangles.push_back({-1 * f1_idx - 1, id2, id0});
                        }
                    }
                    else
                    {
                        if (f0_idx != f1_idx)
                        {
                            border_edges.push_back(std::pair<int, int>(f0_idx + 1, f1_idx + 1));
                            neg_map[id1] = true;
                            pos_map[id0] = true;
                            pos_map[id2] = true;
                            neg.triangles.push_back({id1, -1 * f1_idx - 1, -1 * f0_idx - 1});
                            pos.triangles.push_back({id0, -1 * f0_idx - 1, -1 * f1_idx - 1});
                            pos.triangles.push_back({-1 * f1_idx - 1, id2, id0});
                        }
                        else
                        {
                            pos_map[id0] = true;
                            pos_map[id2] = true;
                            pos.triangles.push_back({-1 * f1_idx - 1, id2, id0});
                        }
                    }
                }
                else if (f1 && f2 && !f0)
                {
                    addEdgePoint(edge_map, border, pi1, id1, id2, idx);
                    addEdgePoint(edge_map, border, pi2, id2, id0, idx);

                    int f1_idx = edge_map[std::pair<int, int>(id1, id2)];
                    int f2_idx = edge_map[std::pair<int, int>(id2, id0)];
                    if (s2 == 1)
                    {
                        if (f2_idx != f1_idx)
                        {
                            border_edges.push_back(std::pair<int, int>(f2_idx + 1, f1_idx + 1));
                            pos_map[id2] = true;
                            neg_map[id0] = true;
                            neg_map[id1] = true;
                            pos.triangles.push_back({id2, -1 * f2_idx - 1, -1 * f1_idx - 1});
                            neg.triangles.push_back({id0, -1 * f1_idx - 1, -1 * f2_idx - 1});
                            neg.triangles.push_back({-1 * f1_idx - 1, id0, id1});
                        }
                        else
                        {
                            neg_map[id0] = true;
                            neg_map[id1] = true;
                            neg.triangles.push_back({-1 * f1_idx - 1, id0, id1});
                        }
                    }
                    else
                    {
                        if (f1_idx != f2_idx)
                        {
                            border_edges.push_back(std::pair<int, int>(f1_idx + 1, f2_idx + 1));
                            neg_map[id2] = true;
                            pos_map[id0] = true;
                            pos_map[id1] = true;
                            neg.triangles.push_back({id2, -1 * f2_idx - 1, -1 * f1_idx - 1});
                            pos.triangles.push_back({id0, -1 * f1_idx - 1, -1 * f2_idx - 1});
                            pos.triangles.push_back({-1 * f1_idx - 1, id0, id1});
                        }
                        else
                        {
                            pos_map[id0] = true;
                            pos_map[id1] = true;
                            pos.triangles.push_back({-1 * f1_idx - 1, id0, id1});
                        }
                    }
                }
                else if (f2 && f0 && !f1)
                {
                    addEdgePoint(edge_map, border, pi2, id2, id0, idx);
                    addEdgePoint(edge_map, border, pi0, id0, id1, idx);

                    int f0_idx = edge_map[std::pair<int, int>(id0, id1)];
                    int f2_idx = edge_map[std::pair<int, int>(id2, id0)];
                    if (s0 == 1)
                    {
                        if (f0_idx != f2_idx)
                        {
                            border_edges.push_back(std::pair<int, int>(f0_idx + 1, f2_idx + 1));
                            pos_map[id0] = true;
                            neg_map[id1] = true;
                            neg_map[id2] = true;
                            pos.triangles.push_back({id0, -1 * f0_idx - 1, -1 * f2_idx - 1});
                            neg.triangles.push_back({id1, -1 * f2_idx - 1, -1 * f0_idx - 1});
                            neg.triangles.push_back({-1 * f2_idx - 1, id1, id2});
                        }
                        else
                        {
                            neg_map[id1] = true;
                            neg_map[id2] = true;
                            neg.triangles.push_back({-1 * f2_idx - 1, id1, id2});
                        }
                    }
                    else
                    {
                        if (f2_idx != f0_idx)
                        {
                            border_edges.push_back(std::pair<int, int>(f2_idx + 1, f0_idx + 1));
                            neg_map[id0] = true;
                            pos_map[id1] = true;
                            pos_map[id2] = true;
                            neg.triangles.push_back({id0, -1 * f0_idx - 1, -1 * f2_idx - 1});
                            pos.triangles.push_back({id1, -1 * f2_idx - 1, -1 * f0_idx - 1});
                            pos.triangles.push_back({-1 * f2_idx - 1, id1, id2});
                        }
                        else
                        {
                            pos_map[id1] = true;
                            pos_map[id2] = true;
                            pos.triangles.push_back({-1 * f2_idx - 1, id1, id2});
                        }
                    }
                }
                else if (f0 && f1 && f2)
                {
                    if (s0 == 0 || (s0 != 0 && s1 != 0 && s2 != 0 && SamePointDetect(pi0, pi2)))
                    {
                        addPoint(vertex_map, border, p0, id0, idx);
                        edge_map[std::pair<int, int>(id0, id1)] = vertex_map[id0];
                        edge_map[std::pair<int, int>(id1, id0)] = vertex_map[id0];
                        edge_map[std::pair<int, int>(id2, id0)] = vertex_map[id0];
                        edge_map[std::pair<int, int>(id0, id2)] = vertex_map[id0];

                        addEdgePoint(edge_map, border, pi1, id1, id2, idx);
                        int f1_idx = edge_map[std::pair<int, int>(id1, id2)];
                        int f0_idx = vertex_map[id0];
                        if (s1 == 1)
                        {
                            if (f1_idx != f0_idx)
                            {
                                border_edges.push_back(std::pair<int, int>(f1_idx + 1, f0_idx + 1));
                                pos_map[id1] = true;
                                neg_map[id2] = true;
                                pos.triangles.push_back({id1, -1 * f1_idx - 1, -1 * f0_idx - 1});
                                neg.triangles.push_back({id2, -1 * f0_idx - 1, -1 * f1_idx - 1});
                            }
                        }
                        else
                        {
                            if (f0_idx != f1_idx)
                            {
                                border_edges.push_back(std::pair<int, int>(f0_idx + 1, f1_idx + 1));
                                neg_map[id1] = true;
                                pos_map[id2] = true;
                                neg.triangles.push_back({id1, -1 * f1_idx - 1, -1 * f0_idx - 1});
                                pos.triangles.push_back({id2, -1 * f0_idx - 1, -1 * f1_idx - 1});
                            }
                        }
                    }
                    else if (s1 == 0 || (s0 != 0 && s1 != 0 && s2 != 0 && SamePointDetect(pi0, pi1)))
                    {
                        addPoint(vertex_map, border, p1, id1, idx);
                        edge_map[std::pair<int, int>(id0, id1)] = vertex_map[id1];
                        edge_map[std::pair<int, int>(id1, id0)] = vertex_map[id1];
                        edge_map[std::pair<int, int>(id1, id2)] = vertex_map[id1];
                        edge_map[std::pair<int, int>(id2, id1)] = vertex_map[id1];

                        addEdgePoint(edge_map, border, pi2, id2, id0, idx);
                        int f1_idx = vertex_map[id1];
                        int f2_idx = edge_map[std::pair<int, int>(id2, id0)];
                        if (s0 == 1)
                        {
                            if (f1_idx != f2_idx)
                            {
                                border_edges.push_back(std::pair<int, int>(f1_idx + 1, f2_idx + 1));
                                pos_map[id0] = true;
                                neg_map[id2] = true;
                                pos.triangles.push_back({id0, -1 * f1_idx - 1, -1 * f2_idx - 1});
                                neg.triangles.push_back({id2, -1 * f2_idx - 1, -1 * f1_idx - 1});
                            }
                        }
                        else
                        {
                            if (f2_idx != f1_idx)
                            {
                                border_edges.push_back(std::pair<int, int>(f2_idx + 1, f1_idx + 1));
                                neg_map[id0] = true;
                                pos_map[id2] = true;
                                neg.triangles.push_back({id0, -1 * f1_idx - 1, -1 * f2_idx - 1});
                                pos.triangles.push_back({id2, -1 * f2_idx - 1, -1 * f1_idx - 1});
                            }
                        }
                    }
                    else if (s2 == 0 || (s0 != 0 && s1 != 0 && s2 != 0 && SamePointDetect(pi1, pi2)))
                    {
                        addPoint(vertex_map, border, p2, id2, idx);
                        edge_map[std::pair<int, int>(id1, id2)] = vertex_map[id2];
                        edge_map[std::pair<int, int>(id2, id1)] = vertex_map[id2];
                        edge_map[std::pair<int, int>(id2, id0)] = vertex_map[id2];
                        edge_map[std::pair<int, int>(id0, id2)] = vertex_map[id2];

                        addEdgePoint(edge_map, border, pi0, id0, id1, idx);
                        int f0_idx = edge_map[std::pair<int, int>(id0, id1)];
                        int f1_idx = vertex_map[id2];
                        if (s0 == 1)
                        {
                            if (f0_idx != f1_idx)
                            {
                                border_edges.push_back(std::pair<int, int>(f0_idx + 1, f1_idx + 1));
                                pos_map[id0] = true;
                                neg_map[id1] = true;
                                pos.triangles.push_back({id0, -1 * f0_idx - 1, -1 * f1_idx - 1});
                                neg.triangles.push_back({id1, -1 * f1_idx - 1, -1 * f0_idx - 1});
                            }
                        }
                        else
                        {
                            if (f1_idx != f0_idx)
                            {
                                border_edges.push_back(std::pair<int, int>(f1_idx + 1, f0_idx + 1));
                                neg_map[id0] = true;
                                pos_map[id1] = true;
                                neg.triangles.push_back({id0, -1 * f0_idx - 1, -1 * f1_idx - 1});
                                pos.triangles.push_back({id1, -1 * f1_idx - 1, -1 * f0_idx - 1});
                            }
                        }
                    }
                    else
                    {
                        throw runtime_error("Intersection error. Please report this error to sarahwei0210@gmail.com with your input OBJ and log file.");
                    }
                }
            }
        }

        if (border.size() > 2)
        {
            int oriN = (int)border.size();
            short flag = Triangulation(border, border_edges, border_triangles, plane);
            if (flag == 0)
            {
#if COACD_USE_CDT_TRIANGULATION
                RemoveOutlierTriangles(border, overlap, border_edges, border_triangles, oriN, border_map, final_border, final_triangles);
#else
                (void)overlap;
                (void)oriN;
                final_border.clear();
                border_map.clear();
                for (int i = 0; i < (int)border.size(); ++i)
                {
                    int mapped_index = -1;
                    for (int j = 0; j < (int)final_border.size(); ++j)
                    {
                        if (SamePointDetect(border[i], final_border[j]))
                        {
                            mapped_index = j + 1;
                            break;
                        }
                    }

                    if (mapped_index < 0)
                    {
                        final_border.push_back(border[i]);
                        mapped_index = static_cast<int>(final_border.size());
                    }

                    border_map[i + 1] = mapped_index;
                }

                final_triangles = border_triangles;
                for (vec3i &triangle : final_triangles)
                {
                    triangle = {
                        border_map[triangle[0]],
                        border_map[triangle[1]],
                        border_map[triangle[2]],
                    };
                }
#endif
            }
            else if (flag == 1)
            {
                final_border = border;
            }
            else
            {
                return false;
            }

            cut_area = 0;
        }
        else
        {
            final_border = border;
            cut_area = -10;
        }

        double pos_x_min = INF, pos_x_max = -INF, pos_y_min = INF, pos_y_max = -INF, pos_z_min = INF, pos_z_max = -INF;
        double neg_x_min = INF, neg_x_max = -INF, neg_y_min = INF, neg_y_max = -INF, neg_z_min = INF, neg_z_max = -INF;

        int pos_idx = 0, neg_idx = 0;
        vector<int> pos_proj(N);
        vector<int> neg_proj(N);
        for (int i = 0; i < N; i++)
        {
            if (pos_map[i] == true)
            {
                pos.points.push_back(mesh.points[i]);
                pos_proj[i] = ++pos_idx;

                pos_x_min = min(pos_x_min, mesh.points[i][0]);
                pos_x_max = max(pos_x_max, mesh.points[i][0]);
                pos_y_min = min(pos_y_min, mesh.points[i][1]);
                pos_y_max = max(pos_y_max, mesh.points[i][1]);
                pos_z_min = min(pos_z_min, mesh.points[i][2]);
                pos_z_max = max(pos_z_max, mesh.points[i][2]);
            }
            if (neg_map[i] == true)
            {
                neg.points.push_back(mesh.points[i]);
                neg_proj[i] = ++neg_idx;

                neg_x_min = min(neg_x_min, mesh.points[i][0]);
                neg_x_max = max(neg_x_max, mesh.points[i][0]);
                neg_y_min = min(neg_y_min, mesh.points[i][1]);
                neg_y_max = max(neg_y_max, mesh.points[i][1]);
                neg_z_min = min(neg_z_min, mesh.points[i][2]);
                neg_z_max = max(neg_z_max, mesh.points[i][2]);
            }
        }

        int pos_N = (int)pos.points.size(), neg_N = (int)neg.points.size();

        for (int i = 0; i < (int)final_border.size(); i++)
        {
            pos.points.push_back(final_border[i]);
            neg.points.push_back(final_border[i]);

            pos_x_min = min(pos_x_min, final_border[i][0]);
            pos_x_max = max(pos_x_max, final_border[i][0]);
            pos_y_min = min(pos_y_min, final_border[i][1]);
            pos_y_max = max(pos_y_max, final_border[i][1]);
            pos_z_min = min(pos_z_min, final_border[i][2]);
            pos_z_max = max(pos_z_max, final_border[i][2]);

            neg_x_min = min(neg_x_min, final_border[i][0]);
            neg_x_max = max(neg_x_max, final_border[i][0]);
            neg_y_min = min(neg_y_min, final_border[i][1]);
            neg_y_max = max(neg_y_max, final_border[i][1]);
            neg_z_min = min(neg_z_min, final_border[i][2]);
            neg_z_max = max(neg_z_max, final_border[i][2]);
        }
        pos.bbox[0] = pos_x_min;
        pos.bbox[1] = pos_x_max;
        pos.bbox[2] = pos_y_min;
        pos.bbox[3] = pos_y_max;
        pos.bbox[4] = pos_z_min;
        pos.bbox[5] = pos_z_max;

        neg.bbox[0] = neg_x_min;
        neg.bbox[1] = neg_x_max;
        neg.bbox[2] = neg_y_min;
        neg.bbox[3] = neg_y_max;
        neg.bbox[4] = neg_z_min;
        neg.bbox[5] = neg_z_max;

        for (int i = 0; i < (int)pos.triangles.size(); i++)
        {
            int f0, f1, f2;
            if (pos.triangles[i][0] >= 0)
                f0 = pos_proj[pos.triangles[i][0]] - 1;
            else
                f0 = pos_N + border_map[-1 * pos.triangles[i][0]] - 1;
            if (pos.triangles[i][1] >= 0)
                f1 = pos_proj[pos.triangles[i][1]] - 1;
            else
                f1 = pos_N + border_map[-1 * pos.triangles[i][1]] - 1;
            if (pos.triangles[i][2] >= 0)
                f2 = pos_proj[pos.triangles[i][2]] - 1;
            else
                f2 = pos_N + border_map[-1 * pos.triangles[i][2]] - 1;

            pos.triangles[i] = {f0, f1, f2};
        }
        for (int i = 0; i < (int)neg.triangles.size(); i++)
        {
            int f0, f1, f2;
            if (neg.triangles[i][0] >= 0)
                f0 = neg_proj[neg.triangles[i][0]] - 1;
            else
                f0 = neg_N + border_map[-1 * neg.triangles[i][0]] - 1;
            if (neg.triangles[i][1] >= 0)
                f1 = neg_proj[neg.triangles[i][1]] - 1;
            else
                f1 = neg_N + border_map[-1 * neg.triangles[i][1]] - 1;
            if (neg.triangles[i][2] >= 0)
                f2 = neg_proj[neg.triangles[i][2]] - 1;
            else
                f2 = neg_N + border_map[-1 * neg.triangles[i][2]] - 1;

            neg.triangles[i] = {f0, f1, f2};
        }

        bool cap_aligns_with_plane = true;
        if (!final_triangles.empty())
        {
            const vec3i &triangle = final_triangles[0];
            const vec3d &p0 = final_border[triangle[0] - 1];
            const vec3d &p1 = final_border[triangle[1] - 1];
            const vec3d &p2 = final_border[triangle[2] - 1];
            const vec3d ab = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
            const vec3d ac = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
            const vec3d normal = CrossProduct(ab, ac);
            const double dot = normal[0] * plane.a + normal[1] * plane.b + normal[2] * plane.c;
            cap_aligns_with_plane = dot >= 0.0;
        }

        for (int i = 0; i < (int)final_triangles.size(); i++)
        {
            cut_area += Area(final_border[final_triangles[i][0] - 1], final_border[final_triangles[i][1] - 1], final_border[final_triangles[i][2] - 1]);
            int a, b, c;
#if COACD_USE_CDT_TRIANGULATION
            a = border_map[final_triangles[i][0]];
            b = border_map[final_triangles[i][1]];
            c = border_map[final_triangles[i][2]];
#else
            a = final_triangles[i][0];
            b = final_triangles[i][1];
            c = final_triangles[i][2];
#endif
            if (cap_aligns_with_plane)
            {
                pos.triangles.push_back({pos_N + c - 1, pos_N + b - 1, pos_N + a - 1});
                neg.triangles.push_back({neg_N + a - 1, neg_N + b - 1, neg_N + c - 1});
            }
            else
            {
                pos.triangles.push_back({pos_N + a - 1, pos_N + b - 1, pos_N + c - 1});
                neg.triangles.push_back({neg_N + c - 1, neg_N + b - 1, neg_N + a - 1});
            }
        }

        return true;
    }
}
