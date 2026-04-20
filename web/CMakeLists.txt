# CoACD WASM/Emscripten build (single-threaded).
# Invoke via emcmake:
#   emcmake cmake -DCMAKE_BUILD_TYPE=MinSizeRel -Bbuild
#   cmake --build build
cmake_minimum_required(VERSION 3.16)

project(coacd_wasm CXX C)

if(NOT EMSCRIPTEN)
    message(FATAL_ERROR "CoACD web build must be configured with emcmake.")
endif()

set(BUILD_TARGET coacd_js)
set(OUTPUT_NAME coacd)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(COACD_WASM_DEBUG "Emit source maps + assertions" OFF)

# Force the parent CoACD build into library-only mode.
set(COACD_LIBRARY_ONLY ON CACHE BOOL "" FORCE)
set(WITH_3RD_PARTY_LIBS OFF CACHE BOOL "" FORCE)
set(COACD_ENABLE_IO OFF CACHE BOOL "" FORCE)
set(COACD_BUILD_CLI OFF CACHE BOOL "" FORCE)
set(COACD_BUILD_SHARED_WRAPPER OFF CACHE BOOL "" FORCE)
# Enable CDT header so the runtime selector compiles; the built-in
# ear-clipping path remains the default at runtime.
set(COACD_USE_CDT_TRIANGULATION ON CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)

# Silence the parent project's -fsanitize=address in Debug.
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "" FORCE)
endif()

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/.. coacd_core)

add_executable(${BUILD_TARGET} binding.cc)

target_link_libraries(${BUILD_TARGET} PRIVATE coacd)
target_include_directories(${BUILD_TARGET} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../public)

set(COACD_COMPILE_FLAGS "-fno-exceptions=0")

if(COACD_WASM_DEBUG)
    set(COACD_LINK_FLAGS
        "-g -gsource-map -sENVIRONMENT=web,worker,node \
         -sMALLOC=emmalloc -sSTACK_SIZE=1048576 \
         -sALLOW_MEMORY_GROWTH=1 -sFILESYSTEM=0 \
         -sMODULARIZE=1 -sEXPORT_ES6=1 \
         -sINVOKE_RUN=0 --bind \
         -sASSERTIONS=1 \
         -sEXPORTED_RUNTIME_METHODS=['HEAPU8','HEAPF64','HEAP32']")
else()
    set(COACD_LINK_FLAGS
        "-Oz -sENVIRONMENT=web,worker,node \
         -sMALLOC=emmalloc -sSTACK_SIZE=524288 \
         -sALLOW_MEMORY_GROWTH=1 -sFILESYSTEM=0 \
         -sMODULARIZE=1 -sEXPORT_ES6=1 \
         -sINVOKE_RUN=0 --bind \
         -sEXPORTED_RUNTIME_METHODS=['HEAPU8','HEAPF64','HEAP32']")
endif()

set_target_properties(${BUILD_TARGET} PROPERTIES
    OUTPUT_NAME ${OUTPUT_NAME}
    SUFFIX ".js"
    LINK_FLAGS "${COACD_LINK_FLAGS}"
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/js/src)
