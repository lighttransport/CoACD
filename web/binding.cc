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

using namespace emscripten;

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

// Read a JS typed array (Float64Array / Int32Array) into a std::vector.
// We use emscripten::convertJSArrayToNumberVector for portability;
// the caller passes a typed array and we fill a numeric vector.
template <typename T>
static std::vector<T> readTypedArray(const val &arr) {
    const size_t n = arr["length"].as<size_t>();
    std::vector<T> out(n);
    val heap = val::module_property(
        std::is_same<T, double>::value ? "HEAPF64" : "HEAP32");
    val memory = heap["buffer"];
    const size_t bytes_per = sizeof(T);
    val dst = val(typed_memory_view(n, out.data()));
    dst.call<void>("set", arr);
    (void)memory;
    (void)bytes_per;
    return out;
}

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

    val Float64Array = val::global("Float64Array");
    val Int32Array = val::global("Int32Array");

    val v_out = Float64Array.new_(verts.size());
    v_out.call<void>("set", val(typed_memory_view(verts.size(), verts.data())));
    val i_out = Int32Array.new_(idx.size());
    i_out.call<void>("set", val(typed_memory_view(idx.size(), idx.data())));

    val obj = val::object();
    obj.set("vertices", v_out);
    obj.set("indices", i_out);
    return obj;
}

static WebOptions parseOptions(const val &js) {
    WebOptions o;
    if (js.isUndefined() || js.isNull()) return o;
    auto pick = [&](const char *k, auto &dst) {
        val v = js[k];
        if (!v.isUndefined() && !v.isNull()) {
            dst = v.as<std::decay_t<decltype(dst)>>();
        }
    };
    pick("threshold", o.threshold);
    pick("max_convex_hull", o.max_convex_hull);
    pick("preprocess", o.preprocess);
    pick("prep_resolution", o.prep_resolution);
    pick("sample_resolution", o.sample_resolution);
    pick("mcts_nodes", o.mcts_nodes);
    pick("mcts_iteration", o.mcts_iteration);
    pick("mcts_max_depth", o.mcts_max_depth);
    pick("pca", o.pca);
    pick("merge", o.merge);
    pick("decimate", o.decimate);
    pick("max_ch_vertex", o.max_ch_vertex);
    pick("extrude", o.extrude);
    pick("extrude_margin", o.extrude_margin);
    pick("apx_mode", o.apx_mode);
    pick("seed", o.seed);
    pick("real_metric", o.real_metric);
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
