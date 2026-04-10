#pragma once
#include <typeinfo>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <thread>
#include <assert.h>

#include <map>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>

#include "nanoflann.hpp"
#include "shape.h"
using namespace std;
using namespace nanoflann;

#define INF numeric_limits<double>::max()

namespace coacd
{

    inline double dist_point2point(vec3d pt, vec3d p)
    {
        double dx = pt[0] - p[0], dy = pt[1] - p[1], dz = pt[2] - p[2];
        return sqrt(dx * dx + dy * dy + dz * dz);
    }

    // Optimized point-to-triangle distance using Voronoi region method.
    // Avoids redundant sqrt/pow calls from the original implementation.
    double dist_point2triangle(vec3d pt, vec3d tri_pt0, vec3d tri_pt1, vec3d tri_pt2, bool flag = false)
    {
        (void)flag;
        // Edge vectors
        double e0[3] = {tri_pt1[0] - tri_pt0[0], tri_pt1[1] - tri_pt0[1], tri_pt1[2] - tri_pt0[2]};
        double e1[3] = {tri_pt2[0] - tri_pt0[0], tri_pt2[1] - tri_pt0[1], tri_pt2[2] - tri_pt0[2]};
        double v[3]  = {tri_pt0[0] - pt[0], tri_pt0[1] - pt[1], tri_pt0[2] - pt[2]};

        double a = e0[0]*e0[0] + e0[1]*e0[1] + e0[2]*e0[2];
        double b = e0[0]*e1[0] + e0[1]*e1[1] + e0[2]*e1[2];
        double c = e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2];
        double d = e0[0]*v[0]  + e0[1]*v[1]  + e0[2]*v[2];
        double e = e1[0]*v[0]  + e1[1]*v[1]  + e1[2]*v[2];

        double det = a * c - b * b;
        double s = b * e - c * d;
        double t = b * d - a * e;

        if (s + t <= det) {
            if (s < 0.0) {
                if (t < 0.0) {
                    // Region 4
                    if (d < 0.0) { s = (-d >= a) ? 1.0 : -d / a; t = 0.0; }
                    else { s = 0.0; t = (e >= 0.0) ? 0.0 : ((-e >= c) ? 1.0 : -e / c); }
                } else {
                    // Region 3
                    s = 0.0; t = (e >= 0.0) ? 0.0 : ((-e >= c) ? 1.0 : -e / c);
                }
            } else if (t < 0.0) {
                // Region 5
                t = 0.0; s = (d >= 0.0) ? 0.0 : ((-d >= a) ? 1.0 : -d / a);
            } else {
                // Region 0 (interior)
                double inv = 1.0 / det;
                s *= inv; t *= inv;
            }
        } else {
            if (s < 0.0) {
                // Region 2
                double tmp0 = b + d, tmp1 = c + e;
                if (tmp1 > tmp0) { double numer = tmp1 - tmp0, denom = a - 2*b + c; s = (numer >= denom) ? 1.0 : numer / denom; t = 1.0 - s; }
                else { s = 0.0; t = (tmp1 <= 0.0) ? 1.0 : ((e >= 0.0) ? 0.0 : -e / c); }
            } else if (t < 0.0) {
                // Region 6
                double tmp0 = b + e, tmp1 = a + d;
                if (tmp1 > tmp0) { double numer = tmp1 - tmp0, denom = a - 2*b + c; t = (numer >= denom) ? 1.0 : numer / denom; s = 1.0 - t; }
                else { t = 0.0; s = (tmp1 <= 0.0) ? 1.0 : ((d >= 0.0) ? 0.0 : -d / a); }
            } else {
                // Region 1
                double numer = (c + e) - (b + d);
                if (numer <= 0.0) { s = 0.0; t = 1.0; }
                else { double denom = a - 2*b + c; s = (numer >= denom) ? 1.0 : numer / denom; t = 1.0 - s; }
            }
        }

        double dx = v[0] + s * e0[0] + t * e1[0];
        double dy = v[1] + s * e0[1] + t * e1[1];
        double dz = v[2] + s * e0[2] + t * e1[2];
        return sqrt(dx * dx + dy * dy + dz * dz);
    }

    // Squared distance variant — avoids sqrt in inner loops.
    inline double dist_point2triangle_sq(const vec3d &pt, const vec3d &tri_pt0, const vec3d &tri_pt1, const vec3d &tri_pt2)
    {
        double e0[3] = {tri_pt1[0] - tri_pt0[0], tri_pt1[1] - tri_pt0[1], tri_pt1[2] - tri_pt0[2]};
        double e1[3] = {tri_pt2[0] - tri_pt0[0], tri_pt2[1] - tri_pt0[1], tri_pt2[2] - tri_pt0[2]};
        double v[3]  = {tri_pt0[0] - pt[0], tri_pt0[1] - pt[1], tri_pt0[2] - pt[2]};

        double a = e0[0]*e0[0] + e0[1]*e0[1] + e0[2]*e0[2];
        double b = e0[0]*e1[0] + e0[1]*e1[1] + e0[2]*e1[2];
        double c = e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2];
        double d = e0[0]*v[0]  + e0[1]*v[1]  + e0[2]*v[2];
        double e = e1[0]*v[0]  + e1[1]*v[1]  + e1[2]*v[2];

        double det = a * c - b * b;
        double s = b * e - c * d;
        double t = b * d - a * e;

        if (s + t <= det) {
            if (s < 0.0) {
                if (t < 0.0) { if (d < 0.0) { s = (-d >= a) ? 1.0 : -d / a; t = 0.0; } else { s = 0.0; t = (e >= 0.0) ? 0.0 : ((-e >= c) ? 1.0 : -e / c); } }
                else { s = 0.0; t = (e >= 0.0) ? 0.0 : ((-e >= c) ? 1.0 : -e / c); }
            } else if (t < 0.0) { t = 0.0; s = (d >= 0.0) ? 0.0 : ((-d >= a) ? 1.0 : -d / a); }
            else { double inv = 1.0 / det; s *= inv; t *= inv; }
        } else {
            if (s < 0.0) { double tmp0 = b + d, tmp1 = c + e; if (tmp1 > tmp0) { double numer = tmp1 - tmp0, denom = a - 2*b + c; s = (numer >= denom) ? 1.0 : numer / denom; t = 1.0 - s; } else { s = 0.0; t = (tmp1 <= 0.0) ? 1.0 : ((e >= 0.0) ? 0.0 : -e / c); } }
            else if (t < 0.0) { double tmp0 = b + e, tmp1 = a + d; if (tmp1 > tmp0) { double numer = tmp1 - tmp0, denom = a - 2*b + c; t = (numer >= denom) ? 1.0 : numer / denom; s = 1.0 - t; } else { t = 0.0; s = (tmp1 <= 0.0) ? 1.0 : ((d >= 0.0) ? 0.0 : -d / a); } }
            else { double numer = (c + e) - (b + d); if (numer <= 0.0) { s = 0.0; t = 1.0; } else { double denom = a - 2*b + c; s = (numer >= denom) ? 1.0 : numer / denom; t = 1.0 - s; } }
        }

        double dx = v[0] + s * e0[0] + t * e1[0];
        double dy = v[1] + s * e0[1] + t * e1[1];
        double dz = v[2] + s * e0[2] + t * e1[2];
        return dx * dx + dy * dy + dz * dz;
    }

    double face_hausdorff_distance(Model &meshA, vector<vec3d> &XA, vector<int> &idA, Model &meshB, vector<vec3d> &XB, vector<int> &idB, bool flag = false)
    {
        (void)flag;
        int nA = static_cast<int>(XA.size());
        int nB = static_cast<int>(XB.size());
        double cmax = 0;

        PointCloud<double> cloudA, cloudB;
        vec2PointCloud(cloudA, XA);
        vec2PointCloud(cloudB, XB);

        typedef KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<double, PointCloud<double>>,
            PointCloud<double>,
            3 /* dim */
            >
            my_kd_tree_t;

        my_kd_tree_t indexA(3 /*dim*/, cloudA, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
        my_kd_tree_t indexB(3 /*dim*/, cloudB, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
        indexA.buildIndex();
        indexB.buildIndex();

        // Stack-allocated KNN result buffers (avoids per-query heap allocation)
        constexpr size_t K = 10;
        size_t ret_index[K];
        double out_dist_sqr[K];

        // Use squared distances in inner loop; sqrt only the final cmin per query.
        for (int i = 0; i < nB; i++)
        {
            double query_pt[3] = {XB[i][0], XB[i][1], XB[i][2]};
            size_t num_results = indexA.knnSearch(&query_pt[0], K, &ret_index[0], &out_dist_sqr[0]);

            double cmin_sq = INF;
            for (size_t j = 0; j < num_results; j++)
            {
                double dsq = dist_point2triangle_sq(XB[i],
                    meshA.points[meshA.triangles[idA[ret_index[j]]][0]],
                    meshA.points[meshA.triangles[idA[ret_index[j]]][1]],
                    meshA.points[meshA.triangles[idA[ret_index[j]]][2]]);
                if (dsq < cmin_sq)
                {
                    cmin_sq = dsq;
                    if (cmin_sq < 1e-28) break;
                }
            }
            double cmin = sqrt(cmin_sq);
            if (cmin > 10)
                cmin = sqrt(out_dist_sqr[0]);
            if (cmin > cmax && INF > cmin)
                cmax = cmin;
        }

        for (int i = 0; i < nA; i++)
        {
            double query_pt[3] = {XA[i][0], XA[i][1], XA[i][2]};
            size_t num_results = indexB.knnSearch(&query_pt[0], K, &ret_index[0], &out_dist_sqr[0]);

            double cmin_sq = INF;
            for (size_t j = 0; j < num_results; j++)
            {
                double dsq = dist_point2triangle_sq(XA[i],
                    meshB.points[meshB.triangles[idB[ret_index[j]]][0]],
                    meshB.points[meshB.triangles[idB[ret_index[j]]][1]],
                    meshB.points[meshB.triangles[idB[ret_index[j]]][2]]);
                if (dsq < cmin_sq)
                {
                    cmin_sq = dsq;
                    if (cmin_sq < 1e-28) break;
                }
            }
            double cmin = sqrt(cmin_sq);
            if (cmin > 10)
                cmin = sqrt(out_dist_sqr[0]);
            if (cmin > cmax && INF > cmin)
                cmax = cmin;
        }

        return cmax;
    }
}