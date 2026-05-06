#include "progress.h"

#include <mutex>

namespace coacd {

namespace {
CoACD_ProgressCallback g_cb = nullptr;
void* g_userdata = nullptr;
std::mutex g_mutex;  // CoACD's MCTS pool runs under OpenMP; without
                     // this the user callback gets concurrent writes.
} // namespace

void EmitProgress(const char* stage, long long current, long long total) {
  std::lock_guard<std::mutex> lock(g_mutex);
  if (g_cb) {
    g_cb(stage, current, total, g_userdata);
  }
}

} // namespace coacd

extern "C" void CoACD_setProgressCallback(CoACD_ProgressCallback cb,
                                          void* userdata) {
  std::lock_guard<std::mutex> lock(coacd::g_mutex);
  coacd::g_cb = cb;
  coacd::g_userdata = userdata;
}
