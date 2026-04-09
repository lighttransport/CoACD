#!/usr/bin/env bash
set -euo pipefail

build_dir="build-library"
if [[ $# -gt 0 && "$1" != -* ]]; then
  build_dir="$1"
  shift
fi

cmake -S . -B "${build_dir}" \
  -DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Release}" \
  -DCOACD_LIBRARY_ONLY=ON \
  -DWITH_3RD_PARTY_LIBS=OFF \
  -DCOACD_ENABLE_IO=OFF \
  -DCOACD_BUILD_CLI=OFF \
  -DCOACD_BUILD_SHARED_WRAPPER=OFF \
  -DCOACD_USE_CDT_TRIANGULATION=OFF \
  "$@"
