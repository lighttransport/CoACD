cmake_minimum_required(VERSION 3.24)
project(CoACD LANGUAGES C CXX)

# ##############################################################################
# General Settings
# ##############################################################################

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(WITH_3RD_PARTY_LIBS "Include 3rd party libraries" ON)
option(COACD_BUILD_CLI "Build the command-line executable" ON)
option(COACD_BUILD_SHARED_WRAPPER
       "Build the shared wrapper used by the Python package" ON)
option(COACD_ENABLE_IO "Enable CLI/file IO code paths" ON)
option(COACD_USE_CDT_TRIANGULATION
       "Use bundled CDT for clipping triangulation instead of the built-in ear clipping path"
       OFF)
option(COACD_LIBRARY_ONLY
       "Configure a minimal library-only build with no 3rd party libs or IO"
       OFF)

if(COACD_LIBRARY_ONLY)
    set(WITH_3RD_PARTY_LIBS OFF CACHE BOOL "Include 3rd party libraries" FORCE)
    set(COACD_ENABLE_IO OFF CACHE BOOL "Enable CLI/file IO code paths" FORCE)
    set(COACD_BUILD_CLI OFF CACHE BOOL "Build the command-line executable" FORCE)
    set(COACD_BUILD_SHARED_WRAPPER OFF CACHE BOOL
        "Build the shared wrapper used by the Python package" FORCE)
endif()

# --- Compiler-specific flags
if(MSVC)
    # For Microsoft Visual Studio (Windows)
    add_compile_options("/bigobj")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0 /D_USE_MATH_DEFINES")
    set(CMAKE_SHARED_LIBRARY_PREFIX "lib")
else()
    # For non-MSVC (Linux, macOS, etc.)
    set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Wall -fpermissive -D_USE_MATH_DEFINES")
    if(APPLE)
        # macOS: add -fPIC
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    elseif(UNIX)
        # Linux/Unix: add -fPIC and disable -Werror=stringop-overflow
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wno-error=stringop-overflow")
    endif()
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
    set(CMAKE_CXX_FLAGS_RELEASE
        "${CMAKE_CXX_FLAGS_RELEASE} -Wno-int-in-bool-context -Wsign-compare")
    set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fsanitize=address")
endif()

if(WITH_3RD_PARTY_LIBS)
    include(zlib)
    include(boost)
    include(openvdb)
    include(spdlog)
endif()

# Set default install postfix
if(UNIX)
    include(GNUInstallDirs)
else()
    set(CMAKE_INSTALL_LIBDIR lib)
    set(CMAKE_INSTALL_INCLUDEDIR include)
    set(CMAKE_INSTALL_BINDIR bin)
endif()

include(CTest)

# ##############################################################################
# BUILD Program
# ##############################################################################

file(GLOB_RECURSE COACD_SRC "src/*.cc" "src/*.cpp")
list(FILTER COACD_SRC EXCLUDE REGEX ".*/io\\.cpp$")
if(NOT WITH_3RD_PARTY_LIBS)
    # Exclude files that match the pattern "preprocess*"
    list(FILTER COACD_SRC EXCLUDE REGEX ".*preprocess.*")
endif()

set(COACD_CDT_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/3rd/cdt/CDT/include")
if(COACD_USE_CDT_TRIANGULATION)
    if(NOT EXISTS "${COACD_CDT_INCLUDE_DIR}/CDT.h")
        message(FATAL_ERROR
                "Missing CDT headers at ${COACD_CDT_INCLUDE_DIR}. "
                "Initialize submodules with `git submodule update --init --recursive`, "
                "or configure with -DCOACD_USE_CDT_TRIANGULATION=OFF.")
    endif()
endif()

find_package(OpenMP)
find_package(Threads)

function(coacd_configure_target target_name)
    set(options USE_CDT)
    cmake_parse_arguments(COACD_TARGET "${options}" "" "" ${ARGN})

    target_include_directories(
      ${target_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/public>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
    if(COACD_TARGET_USE_CDT OR COACD_USE_CDT_TRIANGULATION)
        target_include_directories(${target_name} PRIVATE ${COACD_CDT_INCLUDE_DIR})
        target_compile_definitions(${target_name} PUBLIC COACD_USE_CDT_TRIANGULATION=1)
    else()
        target_compile_definitions(${target_name} PUBLIC COACD_USE_CDT_TRIANGULATION=0)
    endif()
    set_target_properties(${target_name} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

    if(WITH_3RD_PARTY_LIBS)
        target_link_libraries(${target_name} PRIVATE openvdb_static spdlog::spdlog)
        target_compile_definitions(${target_name} PUBLIC WITH_3RD_PARTY_LIBS=1)
    else()
        target_compile_definitions(${target_name} PUBLIC WITH_3RD_PARTY_LIBS=0)
        target_compile_definitions(${target_name} PUBLIC DISABLE_SPDLOG)
    endif()

    if(COACD_ENABLE_IO)
        target_compile_definitions(${target_name} PUBLIC COACD_ENABLE_IO=1)
    else()
        target_compile_definitions(${target_name} PUBLIC COACD_ENABLE_IO=0)
    endif()

    if(OpenMP_CXX_FOUND)
        target_link_libraries(${target_name} PUBLIC OpenMP::OpenMP_CXX)
    endif()

    target_link_libraries(${target_name} PUBLIC Threads::Threads)
endfunction()

set(COACD_LIBRARY_SRC ${COACD_SRC} public/coacd.cpp)

add_library(coacd STATIC ${COACD_LIBRARY_SRC})
coacd_configure_target(coacd)

if(COACD_BUILD_SHARED_WRAPPER)
    add_library(_coacd SHARED ${COACD_LIBRARY_SRC})
    coacd_configure_target(_coacd)
endif()

if(COACD_BUILD_CLI)
    if(NOT COACD_ENABLE_IO)
        message(FATAL_ERROR "COACD_BUILD_CLI requires COACD_ENABLE_IO=ON")
    endif()

    add_executable(main main.cpp src/io.cpp)
    target_link_libraries(main PRIVATE coacd)
endif()

if(BUILD_TESTING)
    add_executable(coacd_test tests/coacd_test.cpp)
    target_include_directories(coacd_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(coacd_test PRIVATE coacd)
    add_test(NAME coacd.unit COMMAND coacd_test)

    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/3rd/cdt/CDT/include/CDT.h")
        add_library(coacd_cdt_testlib STATIC ${COACD_LIBRARY_SRC})
        coacd_configure_target(coacd_cdt_testlib USE_CDT)

        add_executable(clip_compare_builtin tests/clip_compare_runner.cpp)
        target_include_directories(clip_compare_builtin PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
        target_compile_definitions(clip_compare_builtin PRIVATE COACD_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
        target_link_libraries(clip_compare_builtin PRIVATE coacd)

        add_executable(clip_compare_cdt tests/clip_compare_runner.cpp)
        target_include_directories(clip_compare_cdt PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
        target_compile_definitions(clip_compare_cdt PRIVATE COACD_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
        target_link_libraries(clip_compare_cdt PRIVATE coacd_cdt_testlib)

        add_test(
          NAME coacd.clip_compare_examples
          COMMAND ${CMAKE_COMMAND}
                  -DDEFAULT_EXEC=$<TARGET_FILE:clip_compare_builtin>
                  -DCDT_EXEC=$<TARGET_FILE:clip_compare_cdt>
                  -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/compare_clip_outputs.cmake)
    endif()
endif()

# ##############################################################################
# Install target
# ##############################################################################

set(COACD_INSTALL_TARGETS coacd)
if(COACD_BUILD_SHARED_WRAPPER)
    list(APPEND COACD_INSTALL_TARGETS _coacd)
endif()

install(
  TARGETS ${COACD_INSTALL_TARGETS}
  EXPORT CoACDTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(
  EXPORT CoACDTargets
  FILE CoACDTargets.cmake
  NAMESPACE CoACD::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/CoACD)

install(
  FILES public/coacd.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/CoACD)

# ##############################################################################
# Install Cmake Configuration
# ##############################################################################

include(CMakePackageConfigHelpers)

configure_package_config_file(
  ${CMAKE_SOURCE_DIR}/cmake/Config.cmake.in
  "${CMAKE_BINARY_DIR}/CoACDConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/CoACD)

install(FILES "${CMAKE_BINARY_DIR}/CoACDConfig.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/CoACD)
