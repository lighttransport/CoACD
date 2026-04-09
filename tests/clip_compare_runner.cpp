#include "clip.h"
#include "model_obj.h"
#include "process.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

struct ObjMesh {
  std::vector<std::array<double, 3>> vertices;
  std::vector<std::array<int, 3>> indices;
};

int parse_obj_index(std::string_view token, int vertex_count) {
  const size_t slash = token.find('/');
  token = token.substr(0, slash);
  if (token.empty()) {
    throw std::runtime_error("empty OBJ index");
  }

  const int raw_index = std::stoi(std::string(token));
  if (raw_index > 0) {
    return raw_index - 1;
  }
  if (raw_index < 0) {
    return vertex_count + raw_index;
  }
  throw std::runtime_error("OBJ indices are 1-based");
}

ObjMesh load_obj_mesh(const std::string &path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("failed to open " + path);
  }

  ObjMesh mesh;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }

    std::istringstream stream(line);
    std::string prefix;
    stream >> prefix;
    if (prefix == "v") {
      double x = 0.0;
      double y = 0.0;
      double z = 0.0;
      stream >> x >> y >> z;
      mesh.vertices.push_back({x, y, z});
      continue;
    }

    if (prefix == "f") {
      std::vector<int> face;
      std::string token;
      while (stream >> token) {
        face.push_back(parse_obj_index(token, static_cast<int>(mesh.vertices.size())));
      }

      if (face.size() < 3) {
        throw std::runtime_error("face with fewer than 3 vertices in " + path);
      }

      for (size_t i = 1; i + 1 < face.size(); ++i) {
        mesh.indices.push_back({face[0], face[i], face[i + 1]});
      }
    }
  }

  if (mesh.vertices.empty() || mesh.indices.empty()) {
    throw std::runtime_error("no mesh data found in " + path);
  }

  return mesh;
}

coacd::Model load_example_as_model(const std::string &path) {
  const ObjMesh mesh = load_obj_mesh(path);
  coacd::Model model;
  if (!model.Load(mesh.vertices, mesh.indices)) {
    throw std::runtime_error("Model::Load failed for " + path);
  }
  model.Normalize();
  return model;
}

struct ClipSummary {
  int plane_index = -1;
  bool success = false;
  std::int64_t cut_area_q = 0;
  std::int64_t pos_volume_q = 0;
  std::int64_t neg_volume_q = 0;
  size_t pos_points = 0;
  size_t pos_triangles = 0;
  size_t neg_points = 0;
  size_t neg_triangles = 0;
};

std::int64_t quantize(double value) {
  return static_cast<std::int64_t>(std::llround(value * 1.0e8));
}

std::vector<coacd::Plane> candidate_planes() {
  return {
      coacd::Plane(1.0, 0.0, 0.0, 0.0),
      coacd::Plane(0.0, 1.0, 0.0, 0.0),
      coacd::Plane(0.0, 0.0, 1.0, 0.0),
      coacd::Plane(1.0, 1.0, 1.0, 0.0),
      coacd::Plane(1.0, 2.0, 3.0, 0.0),
      coacd::Plane(2.0, -1.0, 1.0, 0.0),
  };
}

std::optional<int> forced_plane_index() {
  const char *value = std::getenv("COACD_FORCE_PLANE_INDEX");
  if (value == nullptr) {
    return std::nullopt;
  }
  return std::stoi(value);
}

ClipSummary summarize_clip(const coacd::Model &input) {
  const std::vector<coacd::Plane> planes = candidate_planes();
  const std::optional<int> forced_index = forced_plane_index();

  for (size_t plane_index = 0; plane_index < planes.size(); ++plane_index) {
    if (forced_index.has_value() &&
        static_cast<int>(plane_index) != forced_index.value()) {
      continue;
    }

    coacd::Model pos;
    coacd::Model neg;
    coacd::Plane plane = planes[plane_index];
    double cut_area = 0.0;

    if (!coacd::Clip(input, pos, neg, plane, cut_area)) {
      continue;
    }
    if (pos.points.empty() || neg.points.empty() || pos.triangles.empty() ||
        neg.triangles.empty() || cut_area <= 0.0) {
      continue;
    }

    ClipSummary summary;
    summary.plane_index = static_cast<int>(plane_index);
    summary.success = true;
    summary.cut_area_q = quantize(cut_area);
    summary.pos_volume_q = quantize(std::fabs(coacd::MeshVolume(pos)));
    summary.neg_volume_q = quantize(std::fabs(coacd::MeshVolume(neg)));
    summary.pos_points = pos.points.size();
    summary.pos_triangles = pos.triangles.size();
    summary.neg_points = neg.points.size();
    summary.neg_triangles = neg.triangles.size();
    return summary;
  }

  return {};
}

void print_summary_line(const std::string &name, const ClipSummary &summary) {
  std::cout << name
            << "|success=" << (summary.success ? 1 : 0)
            << "|plane=" << summary.plane_index
            << "|cut=" << summary.cut_area_q
            << "|pos_volume=" << summary.pos_volume_q
            << "|neg_volume=" << summary.neg_volume_q
            << "|pos_points=" << summary.pos_points
            << "|pos_tris=" << summary.pos_triangles
            << "|neg_points=" << summary.neg_points
            << "|neg_tris=" << summary.neg_triangles
            << '\n';
}

}  // namespace

int main() {
  const std::vector<std::string> examples = {
      "Bottle.obj",
      "Kettle.obj",
      "KitchenPot.obj",
      "Octocat-v2.obj",
      "SnowFlake.obj",
  };

  try {
    for (const std::string &example : examples) {
      const std::string path = std::string(COACD_SOURCE_DIR) + "/examples/" + example;
      const coacd::Model model = load_example_as_model(path);
      const ClipSummary summary = summarize_clip(model);
      print_summary_line(example, summary);
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }

  return 0;
}
