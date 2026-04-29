#include "coacd.h"
#include "clip.h"
#include "model_obj.h"
#include "process.h"

#include <cmath>
#include <cstdlib>
#include <map>
#include <exception>
#include <functional>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using TestFn = std::function<void()>;

struct TestFailure : std::runtime_error {
  using std::runtime_error::runtime_error;
};

void fail(const std::string &message) { throw TestFailure(message); }

void expect_true(bool condition, const std::string &message) {
  if (!condition) {
    fail(message);
  }
}

void expect_eq(size_t actual, size_t expected, const std::string &message) {
  if (actual != expected) {
    fail(message + " expected=" + std::to_string(expected) +
         " actual=" + std::to_string(actual));
  }
}

void expect_near(double actual, double expected, double eps,
                 const std::string &message) {
  if (std::fabs(actual - expected) > eps) {
    fail(message + " expected=" + std::to_string(expected) +
         " actual=" + std::to_string(actual));
  }
}

void expect_throws(const TestFn &fn, const std::string &contains) {
  try {
    fn();
  } catch (const std::exception &e) {
    if (contains.empty() || std::string(e.what()).find(contains) != std::string::npos) {
      return;
    }
    fail("unexpected exception message: " + std::string(e.what()));
  }
  fail("expected exception was not thrown");
}

bool clip_test_debug_enabled() {
  return std::getenv("COACD_DEBUG_CLIP_TEST") != nullptr;
}

void dump_model(const char *label, const coacd::Model &model) {
  if (!clip_test_debug_enabled()) {
    return;
  }

  std::cerr << label << " points=" << model.points.size()
            << " triangles=" << model.triangles.size() << '\n';
  for (size_t i = 0; i < model.points.size(); ++i) {
    const auto &p = model.points[i];
    std::cerr << "  v[" << i << "] = (" << p[0] << ", " << p[1] << ", " << p[2]
              << ")\n";
  }
  for (size_t i = 0; i < model.triangles.size(); ++i) {
    const auto &t = model.triangles[i];
    std::cerr << "  f[" << i << "] = (" << t[0] << ", " << t[1] << ", " << t[2]
              << ")\n";
  }

  std::map<std::pair<int, int>, int> edge_counts;
  for (const auto &t : model.triangles) {
    edge_counts[{t[0], t[1]}]++;
    edge_counts[{t[1], t[2]}]++;
    edge_counts[{t[2], t[0]}]++;
  }

  std::set<std::pair<int, int>> unmatched;
  for (const auto &[edge, count] : edge_counts) {
    (void)count;
    if (!edge_counts.contains({edge.second, edge.first})) {
      unmatched.insert(edge);
    }
  }

  if (!unmatched.empty()) {
    std::cerr << "  unmatched directed edges:\n";
    for (const auto &edge : unmatched) {
      std::cerr << "    (" << edge.first << ", " << edge.second << ")\n";
    }
  }
}

coacd::Mesh make_cube_mesh() {
  coacd::Mesh mesh;
  mesh.vertices = {
      {{-1.0, -1.0, -1.0}},
      {{1.0, -1.0, -1.0}},
      {{1.0, 1.0, -1.0}},
      {{-1.0, 1.0, -1.0}},
      {{-1.0, -1.0, 1.0}},
      {{1.0, -1.0, 1.0}},
      {{1.0, 1.0, 1.0}},
      {{-1.0, 1.0, 1.0}},
  };
  mesh.indices = {
      {{0, 2, 1}}, {{0, 3, 2}}, {{4, 5, 6}}, {{4, 6, 7}},
      {{0, 1, 5}}, {{0, 5, 4}}, {{1, 2, 6}}, {{1, 6, 5}},
      {{2, 3, 7}}, {{2, 7, 6}}, {{3, 0, 4}}, {{3, 4, 7}},
  };
  return mesh;
}

coacd::Mesh make_tetrahedron_mesh() {
  coacd::Mesh mesh;
  mesh.vertices = {
      {{0.0, 0.0, 0.0}},
      {{1.0, 0.0, 0.0}},
      {{0.0, 1.0, 0.0}},
      {{0.0, 0.0, 1.0}},
  };
  mesh.indices = {
      {{0, 2, 1}},
      {{0, 1, 3}},
      {{1, 2, 3}},
      {{2, 0, 3}},
  };
  return mesh;
}

coacd::Model to_model(const coacd::Mesh &mesh) {
  coacd::Model model;
  expect_true(model.Load(mesh.vertices, mesh.indices), "failed to load mesh into Model");
  return model;
}

void test_is_manifold_rejects_empty_mesh() {
  coacd::set_log_level("off");
  coacd::Model model;

  expect_true(!coacd::IsManifold(model), "empty mesh should not be manifold");
}

void test_is_manifold_rejects_out_of_range_triangle() {
  coacd::set_log_level("off");
  coacd::Model model;
  model.points = {
      {{0.0, 0.0, 0.0}},
      {{1.0, 0.0, 0.0}},
      {{0.0, 1.0, 0.0}},
  };
  model.triangles = {{{0, 1, 3}}};

  expect_true(!coacd::IsManifold(model),
              "mesh with out-of-range triangle index should not be manifold");
}

void test_public_api_convex_mesh_returns_single_part() {
  coacd::set_log_level("off");
  const coacd::Mesh mesh = make_tetrahedron_mesh();
  const auto parts = coacd::CoACD(mesh, 0.2, -1, "off", 50, 256, 8, 32, 3,
                                  false, true, false, 64, false, 0.01, "ch",
                                  1234, false);

  expect_eq(parts.size(), 1, "convex tetrahedron should stay as one part");
  expect_eq(parts[0].indices.size(), 4, "tetrahedron hull should keep four faces");
}

void test_public_api_rejects_invalid_threshold() {
  coacd::set_log_level("off");
  const coacd::Mesh mesh = make_tetrahedron_mesh();
  expect_throws(
      [&]() {
        (void)coacd::CoACD(mesh, 1.5, -1, "off", 50, 128, 8, 16, 2, false,
                           true, false, 64, false, 0.01, "ch", 1234, false);
      },
      "threshold > 1");
}

void test_public_api_rejects_empty_mesh() {
  coacd::set_log_level("off");
  const coacd::Mesh mesh;

  expect_throws(
      [&]() {
        (void)coacd::CoACD(mesh, 0.2, -1, "off", 50, 128, 8, 16, 2, false,
                           true, false, 64, false, 0.01, "ch", 1234, false);
      },
      "mesh is empty");
}

void test_public_api_rejects_out_of_range_triangle() {
  coacd::set_log_level("off");
  coacd::Mesh mesh;
  mesh.vertices = {
      {{0.0, 0.0, 0.0}},
      {{1.0, 0.0, 0.0}},
      {{0.0, 1.0, 0.0}},
  };
  mesh.indices = {{{0, 1, 3}}};

  expect_throws(
      [&]() {
        (void)coacd::CoACD(mesh, 0.2, -1, "off", 50, 128, 8, 16, 2, false,
                           true, false, 64, false, 0.01, "ch", 1234, false);
      },
      "triangle index out of range");
}

void test_c_api_convex_mesh_returns_single_part() {
  CoACD_setLogLevel("off");
  const coacd::Mesh mesh = make_tetrahedron_mesh();

  std::vector<double> vertices;
  vertices.reserve(mesh.vertices.size() * 3);
  for (const auto &vertex : mesh.vertices) {
    vertices.push_back(vertex[0]);
    vertices.push_back(vertex[1]);
    vertices.push_back(vertex[2]);
  }

  std::vector<int> triangles;
  triangles.reserve(mesh.indices.size() * 3);
  for (const auto &triangle : mesh.indices) {
    triangles.push_back(triangle[0]);
    triangles.push_back(triangle[1]);
    triangles.push_back(triangle[2]);
  }

  CoACD_Mesh input{};
  input.vertices_ptr = vertices.data();
  input.vertices_count = mesh.vertices.size();
  input.triangles_ptr = triangles.data();
  input.triangles_count = mesh.indices.size();

  const CoACD_MeshArray result =
      CoACD_run(input, 0.2, -1, preprocess_off, 50, 256, 8, 32, 3, false,
                true, false, 64, false, 0.01, apx_ch, 1234, false);

  expect_eq(result.meshes_count, 1, "C API tetrahedron should stay as one part");
  expect_eq(result.meshes_ptr[0].triangles_count, 4,
            "C API tetrahedron hull should keep four faces");
  CoACD_freeMeshArray(result);
}

void test_clip_splits_cube_into_closed_parts() {
  const coacd::Mesh cube = make_cube_mesh();
  coacd::Model input = to_model(cube);

  coacd::Model pos;
  coacd::Model neg;
  coacd::Plane plane(1.0, 0.0, 0.0, 0.0);
  double cut_area = 0.0;

  const bool ok = coacd::Clip(input, pos, neg, plane, cut_area);
  expect_true(ok, "Clip should succeed on a cube");
  expect_true(cut_area > 0.0, "Clip should generate a positive cut area");
  expect_true(!pos.points.empty() && !neg.points.empty(),
              "Clip should populate both output meshes");
  dump_model("pos", pos);
  dump_model("neg", neg);
  expect_true(coacd::IsManifold(pos), "positive clipped mesh should be manifold");
  expect_true(coacd::IsManifold(neg), "negative clipped mesh should be manifold");
  expect_near(cut_area, 4.0, 1e-6, "cube center cut area should be 4");
}

}  // namespace

int main() {
  const std::vector<std::pair<std::string, TestFn>> tests = {
      {"public_api_convex_mesh_returns_single_part",
       test_public_api_convex_mesh_returns_single_part},
      {"public_api_rejects_invalid_threshold",
       test_public_api_rejects_invalid_threshold},
      {"public_api_rejects_empty_mesh", test_public_api_rejects_empty_mesh},
      {"public_api_rejects_out_of_range_triangle",
       test_public_api_rejects_out_of_range_triangle},
      {"c_api_convex_mesh_returns_single_part",
       test_c_api_convex_mesh_returns_single_part},
      {"is_manifold_rejects_empty_mesh", test_is_manifold_rejects_empty_mesh},
      {"is_manifold_rejects_out_of_range_triangle",
       test_is_manifold_rejects_out_of_range_triangle},
      {"clip_splits_cube_into_closed_parts", test_clip_splits_cube_into_closed_parts},
  };

  int failures = 0;
  for (const auto &[name, test] : tests) {
    try {
      test();
      std::cout << "[PASS] " << name << '\n';
    } catch (const std::exception &e) {
      ++failures;
      std::cerr << "[FAIL] " << name << ": " << e.what() << '\n';
    } catch (...) {
      ++failures;
      std::cerr << "[FAIL] " << name << ": unknown exception\n";
    }
  }

  if (failures != 0) {
    std::cerr << failures << " test(s) failed\n";
    return EXIT_FAILURE;
  }

  std::cout << tests.size() << " test(s) passed\n";
  return EXIT_SUCCESS;
}
