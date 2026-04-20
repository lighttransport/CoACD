// SPDX-License-Identifier: Apache-2.0
// CoACD embind bindings for browser/Node WASM.
//
// JS-facing API:
//
//   const mod = await (await import('./coacd.js')).default();
//   const hulls = mod.decompose({
//     vertices: Float64Array,   // flat xyz
//     indices:  Int32Array,     // flat triangle indices
//     options:  { threshold: 0.05, max_convex_hull: -1, ... }  (optional)
//   });
//   // hulls: Array<{ vertices: Float64Array, indices: Int32Array }>
//
//   mod.set_log_level("off" | "info" | "warn" | "error");
//
// Threading / preprocessing: this build is single-threaded and has
// OpenVDB disabled. The `preprocess` option defaults to "off"; passing
// "on" / "auto" will likely fail at runtime with an error log.

#include <emscripten/bind.h>

#include <cstdint>
#include <string>
#include <vector>

#include "coacd.h"
#include "embind-utils.hpp"

using namespace emscripten;
using embind_utils::float64ArrayFromVector;
using embind_utils::int32ArrayFromVector;
using embind_utils::pick;
using embind_utils::readTypedArray;

namespace {

struct WebOptions {
    double threshold = 0.05;
    int max_convex_hull = -1;
    std::string preprocess = "off";
    int prep_resolution = 50;
    int sample_resolution = 2000;
    int mcts_nodes = 20;
    int mcts_iteration = 150;
    int mcts_max_depth = 3;
    bool pca = false;
    bool merge = true;
    bool decimate = false;
    int max_ch_vertex = 256;
    bool extrude = false;
    double extrude_margin = 0.01;
    std::string apx_mode = "ch";
    unsigned int seed = 0;
    bool real_metric = false;
};

static coacd::Mesh buildInputMesh(const val &input) {
    val v = input["vertices"];
    val f = input["indices"];
    auto verts = readTypedArray<double>(v);
    auto idx = readTypedArray<int32_t>(f);

    coacd::Mesh m;
    if (verts.size() % 3u != 0) {
        return m;
    }
    if (idx.size() % 3u != 0) {
        return m;
    }
    m.vertices.resize(verts.size() / 3u);
    for (size_t i = 0; i < m.vertices.size(); ++i) {
        m.vertices[i] = {verts[3 * i], verts[3 * i + 1], verts[3 * i + 2]};
    }
    m.indices.resize(idx.size() / 3u);
    for (size_t i = 0; i < m.indices.size(); ++i) {
        m.indices[i] = {idx[3 * i], idx[3 * i + 1], idx[3 * i + 2]};
    }
    return m;
}

static val hullToJS(const coacd::Mesh &m) {
    std::vector<double> verts;
    verts.reserve(m.vertices.size() * 3u);
    for (const auto &v : m.vertices) {
        verts.push_back(v[0]);
        verts.push_back(v[1]);
        verts.push_back(v[2]);
    }
    std::vector<int32_t> idx;
    idx.reserve(m.indices.size() * 3u);
    for (const auto &t : m.indices) {
        idx.push_back(t[0]);
        idx.push_back(t[1]);
        idx.push_back(t[2]);
    }

    val obj = val::object();
    obj.set("vertices", float64ArrayFromVector(verts));
    obj.set("indices", int32ArrayFromVector(idx));
    return obj;
}

static WebOptions parseOptions(const val &js) {
    WebOptions o;
    if (js.isUndefined() || js.isNull()) return o;
    pick(js, "threshold", o.threshold);
    pick(js, "max_convex_hull", o.max_convex_hull);
    pick(js, "preprocess", o.preprocess);
    pick(js, "prep_resolution", o.prep_resolution);
    pick(js, "sample_resolution", o.sample_resolution);
    pick(js, "mcts_nodes", o.mcts_nodes);
    pick(js, "mcts_iteration", o.mcts_iteration);
    pick(js, "mcts_max_depth", o.mcts_max_depth);
    pick(js, "pca", o.pca);
    pick(js, "merge", o.merge);
    pick(js, "decimate", o.decimate);
    pick(js, "max_ch_vertex", o.max_ch_vertex);
    pick(js, "extrude", o.extrude);
    pick(js, "extrude_margin", o.extrude_margin);
    pick(js, "apx_mode", o.apx_mode);
    pick(js, "seed", o.seed);
    pick(js, "real_metric", o.real_metric);
    return o;
}

val decompose(val input) {
    coacd::Mesh in_mesh = buildInputMesh(input);

    val opts_js = input["options"];
    WebOptions o = parseOptions(opts_js);

    std::vector<coacd::Mesh> hulls = coacd::CoACD(
        in_mesh, o.threshold, o.max_convex_hull, o.preprocess,
        o.prep_resolution, o.sample_resolution, o.mcts_nodes,
        o.mcts_iteration, o.mcts_max_depth, o.pca, o.merge, o.decimate,
        o.max_ch_vertex, o.extrude, o.extrude_margin, o.apx_mode, o.seed,
        o.real_metric);

    val arr = val::array();
    for (size_t i = 0; i < hulls.size(); ++i) {
        arr.set(i, hullToJS(hulls[i]));
    }
    return arr;
}

void set_log_level(std::string level) {
    coacd::set_log_level(level);
}

}  // namespace

EMSCRIPTEN_BINDINGS(coacd_module) {
    function("decompose", &decompose);
    function("set_log_level", &set_log_level);
}
