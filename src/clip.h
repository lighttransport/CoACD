#pragma once

#include "model_obj.h"
#include <deque>
#include <unordered_map>

namespace coacd {
    struct PairHash {
        std::size_t operator()(const std::pair<int,int>& p) const {
            auto h1 = std::hash<int>{}(p.first);
            auto h2 = std::hash<int>{}(p.second);
            return h1 ^ (h2 * 2654435761u);
        }
    };
    template<typename V>
    using EdgeMap = std::unordered_map<std::pair<int,int>, V, PairHash>;
}

#ifndef COACD_USE_CDT_TRIANGULATION
#define COACD_USE_CDT_TRIANGULATION 0
#endif

// Runtime triangulation method selection.
// 0 = use compile-time default (CDT if COACD_USE_CDT_TRIANGULATION, else built-in)
// 1 = force CDT
// 2 = force built-in ear clipping
#ifndef COACD_CLIP_TRIANGULATION_RUNTIME
#define COACD_CLIP_TRIANGULATION_RUNTIME 1
#endif

#if COACD_CLIP_TRIANGULATION_RUNTIME
#include <atomic>
namespace coacd {
    enum ClipTriangulationMethod {
        kClipTriangulationDefault = 0,
        kClipTriangulationCDT = 1,
        kClipTriangulationBuiltin = 2,
    };
    // Global runtime selection.  Thread-safe to read; set before calling CoACD.
    void SetClipTriangulationMethod(ClipTriangulationMethod method);
    ClipTriangulationMethod GetClipTriangulationMethod();
}
#endif

using std::deque;
using std::endl;

namespace coacd
{
    // Runtime triangulation method query (always available).
    bool UseCDT();
    void SimpleCyclesFromEdges(const vector<pair<int, int>> &edges, vector<vector<int>> &simple_cycles);
    void FindCycleDirection(const vector<vec3d> &border, const vector<vector<int>> &cycles, const Plane &plane,
                            EdgeMap<bool> &cycles_dir);
    void RemoveOutlierTriangles(const vector<vec3d> &border, const vector<vec3d> &overlap,
                                const vector<pair<int, int>> &border_edges, const vector<vec3i> &border_triangles,
                                int oriN, std::unordered_map<int, int> &vertex_map, vector<vec3d> &final_border,
                                vector<vec3i> &final_triangles);
    void RecordIntersectPoint(Model mesh, EdgeMap<int> &edge_map, int i, int ep0, int ep1, int &idx, vector<vec3d> &border, vec3d point);
    bool Clip(const Model &mesh, Model &pos, Model &neg, Plane &plane, double &cut_area, bool foo = false);
    bool CreatePlaneRotationMatrix(const vector<vec3d> &border, const vector<pair<int, int>> &border_edges, vec3d &T,
                                   double R[3][3], Plane &plane);
    short Triangulation(vector<vec3d> &border, const vector<pair<int, int>> &border_edges,
                        vector<vec3i> &border_triangles, Plane &plane);
    void PrintEdgeSet(const vector<pair<int, int>> &edges);

    inline void addPoint(std::unordered_map<int, int> &vertex_map, vector<vec3d> &border, const vec3d &pt, int id, int &idx)
    {
        if (vertex_map.find(id) == vertex_map.end())
        {
            int flag = -1;
            for (int i = 0; i < (int)border.size(); i++)
            {
                if ((fabs(border[i][0] - pt[0])) < 1e-4 && (fabs(border[i][1] - pt[1])) < 1e-4 && (fabs(border[i][2] - pt[2])) < 1e-4)
                {
                    flag = i;
                    break;
                }
            }
            if (flag == -1)
            {
                vertex_map[id] = idx;
                border.push_back(pt);
                idx++;
            }
            else
                vertex_map[id] = flag;
        }
    }

    inline void addEdgePoint(EdgeMap<int> &edge_map, vector<vec3d> &border, const vec3d &pt, int id1, int id2, int &idx)
    {
        pair<int, int> edge1 = make_pair(id1, id2);
        pair<int, int> edge2 = make_pair(id2, id1);
        if (edge_map.find(edge1) == edge_map.end() && edge_map.find(edge2) == edge_map.end())
        {
            int flag = -1;
            for (int i = 0; i < (int)border.size(); i++)
            {
                if ((fabs(border[i][0] - pt[0])) < 1e-4 && (fabs(border[i][1] - pt[1])) < 1e-4 && (fabs(border[i][2] - pt[2])) < 1e-4)
                {
                    flag = i;
                    break;
                }
            }
            if (flag == -1)
            {
                edge_map[edge1] = idx;
                edge_map[edge2] = idx;
                border.push_back(pt);
                idx++;
            }
            else
            {
                edge_map[edge1] = flag;
                edge_map[edge2] = flag;
            }
        }
    }

    inline bool FaceOverlap(const std::unordered_map<int, bool> &overlap_map, const vec3i &triangle)
    {
        int idx0 = triangle[0], idx1 = triangle[1], idx2 = triangle[2];
        if (overlap_map.find(idx0) == overlap_map.end() && overlap_map.find(idx1) == overlap_map.end() &&
            overlap_map.find(idx2) == overlap_map.end())
            return false;
        return true;
    }
}
