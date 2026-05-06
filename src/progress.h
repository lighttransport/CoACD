#pragma once

#include "../public/coacd.h"

namespace coacd {

// Emit a progress event to the registered callback (no-op if none).
// Called synchronously from CoACD's working thread. Cheap when no
// callback is installed, so it's safe to sprinkle freely inside hot
// loops (e.g. the MCTS decomposition pool).
void EmitProgress(const char* stage, long long current, long long total);

} // namespace coacd
