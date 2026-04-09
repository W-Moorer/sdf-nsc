cmake_minimum_required(VERSION 3.19)
project(baseline_contact_platform LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(SPCC_ENABLE_VDB "Enable OpenVDB + NanoVDB integration" ON)
set(NANOVDB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/third_party/openvdb/nanovdb" CACHE PATH "Path to NanoVDB root (containing nanovdb/NanoVDB.h)")

# Ensure CMake looks in the internal build deps folder FIRST if not externally provided
if(NOT Chrono_DIR)
    set(Chrono_DIR "${CMAKE_CURRENT_SOURCE_DIR}/_deps/chrono-install/lib/cmake/Chrono" CACHE PATH "Path to Project Chrono Config")
endif()

find_package(Chrono CONFIG REQUIRED)

if(NOT Chrono_FOUND)
    message(FATAL_ERROR 
    "Project Chrono not found! \n"
    "Please run: bash scripts/bootstrap_all.sh \n"
    "Or execute bootstrap_chrono.sh and build_chrono_minimal.sh sequentially.")
endif()

message(STATUS "Successfully found Project Chrono at: ${Chrono_DIR}")

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# Common library (Backend implementation)
set(PLATFORM_BACKEND_SOURCES
    src/backend/chrono/ChronoRigidSystemNSC.cpp
    src/backend/chrono/spcc/CompressedContactPipeline.cpp
    src/backend/chrono/spcc/CompressedContactValidation.cpp
    src/backend/chrono/spcc/DenseSurfaceSampler.cpp
    src/models/BaselineDropCase.cpp
    src/models/BaselineHeadOnSphereCase.cpp
    src/models/BaselineRevoluteClearanceCase.cpp
    src/models/BaselineCamCase.cpp
    src/models/BaselineSimpleGearCase.cpp
    src/validation/ReferenceCurve.cpp
    src/validation/ErrorMetrics.cpp
)

if(SPCC_ENABLE_VDB)
    list(APPEND PLATFORM_BACKEND_SOURCES
        src/backend/chrono/spcc/VDBSDFField.cpp
    )
endif()

add_library(platform_backend ${PLATFORM_BACKEND_SOURCES})
target_link_libraries(platform_backend PUBLIC Chrono::Chrono_core)
target_include_directories(platform_backend PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

if(SPCC_ENABLE_VDB)
    if(DEFINED OPENVDB_INCLUDE_DIR)
        unset(OPENVDB_INCLUDE_DIR CACHE)
    endif()
    if(DEFINED OPENVDB_LIBRARY)
        unset(OPENVDB_LIBRARY CACHE)
    endif()
    if(DEFINED TBB_LIBRARY)
        unset(TBB_LIBRARY CACHE)
    endif()

    find_path(OPENVDB_INCLUDE_DIR
        openvdb/openvdb.h
        PATHS /usr/include /usr/local/include
        NO_DEFAULT_PATH)
    if(NOT OPENVDB_INCLUDE_DIR)
        find_path(OPENVDB_INCLUDE_DIR openvdb/openvdb.h)
    endif()

    find_library(OPENVDB_LIBRARY
        openvdb
        PATHS /usr/lib/x86_64-linux-gnu /usr/local/lib
        NO_DEFAULT_PATH)
    if(NOT OPENVDB_LIBRARY)
        find_library(OPENVDB_LIBRARY openvdb)
    endif()

    find_library(TBB_LIBRARY
        tbb
        PATHS /usr/lib/x86_64-linux-gnu /lib/x86_64-linux-gnu /usr/local/lib
        NO_DEFAULT_PATH)
    if(NOT TBB_LIBRARY)
        find_library(TBB_LIBRARY tbb)
    endif()

    if(NOT OPENVDB_INCLUDE_DIR OR NOT OPENVDB_LIBRARY OR NOT TBB_LIBRARY)
        message(FATAL_ERROR
            "SPCC_ENABLE_VDB=ON but OpenVDB/TBB was not found.\n"
            "Install it first (Ubuntu: sudo apt-get install -y libopenvdb-dev libtbb-dev).")
    endif()

    get_filename_component(OPENVDB_LIBRARY_DIR "${OPENVDB_LIBRARY}" DIRECTORY)
    set(SPCC_VDB_RUNTIME_DIR "${OPENVDB_LIBRARY_DIR}")

    if(NOT EXISTS "${NANOVDB_ROOT}/nanovdb/NanoVDB.h")
        message(FATAL_ERROR
            "SPCC_ENABLE_VDB=ON but NanoVDB header not found at: ${NANOVDB_ROOT}/nanovdb/NanoVDB.h\n"
            "Run: bash scripts/bootstrap_vdb.sh")
    endif()

    set(OPENVDB_COORD_HEADER "${OPENVDB_INCLUDE_DIR}/openvdb/math/Coord.h")
    if(NOT EXISTS "${OPENVDB_COORD_HEADER}")
        message(FATAL_ERROR
            "SPCC_ENABLE_VDB=ON but OpenVDB Coord.h was not found at: ${OPENVDB_COORD_HEADER}")
    endif()

    set(OPENVDB_COMPAT_ROOT "${CMAKE_CURRENT_BINARY_DIR}/openvdb_compat")
    set(OPENVDB_COMPAT_COORD_HEADER "${OPENVDB_COMPAT_ROOT}/openvdb/math/Coord.h")
    set(OPENVDB_COMPAT_MESH_TO_VOLUME_HEADER "${OPENVDB_COMPAT_ROOT}/openvdb/tools/MeshToVolume.h")
    set(OPENVDB_COMPAT_POINT_INDEX_GRID_HEADER "${OPENVDB_COMPAT_ROOT}/openvdb/tools/PointIndexGrid.h")
    file(MAKE_DIRECTORY "${OPENVDB_COMPAT_ROOT}/openvdb/math")
    file(MAKE_DIRECTORY "${OPENVDB_COMPAT_ROOT}/openvdb/tools")
    file(MAKE_DIRECTORY "${OPENVDB_COMPAT_ROOT}/openvdb/io")
    file(MAKE_DIRECTORY "${OPENVDB_COMPAT_ROOT}/tbb")

    file(WRITE "${OPENVDB_COMPAT_ROOT}/openvdb/openvdb.h"
"#pragma once\n"
"#include \"${OPENVDB_INCLUDE_DIR}/openvdb/openvdb.h\"\n")

    file(WRITE "${OPENVDB_COMPAT_ROOT}/openvdb/io/File.h"
"#pragma once\n"
"#include \"${OPENVDB_INCLUDE_DIR}/openvdb/io/File.h\"\n")

    file(READ "${OPENVDB_COORD_HEADER}" OPENVDB_COORD_CONTENT)
    string(REPLACE "namespace tbb { class split; } // forward declaration\n" "" OPENVDB_COORD_CONTENT "${OPENVDB_COORD_CONTENT}")
    string(REPLACE "#include <openvdb/Platform.h>" "#include <openvdb/Platform.h>\n#include <tbb/blocked_range.h>" OPENVDB_COORD_CONTENT "${OPENVDB_COORD_CONTENT}")
    string(REPLACE "#include \"Math.h\"" "#include \"${OPENVDB_INCLUDE_DIR}/openvdb/math/Math.h\"" OPENVDB_COORD_CONTENT "${OPENVDB_COORD_CONTENT}")
    string(REPLACE "#include \"Vec3.h\"" "#include \"${OPENVDB_INCLUDE_DIR}/openvdb/math/Vec3.h\"" OPENVDB_COORD_CONTENT "${OPENVDB_COORD_CONTENT}")
    file(WRITE "${OPENVDB_COMPAT_COORD_HEADER}" "${OPENVDB_COORD_CONTENT}")

    file(READ "${OPENVDB_INCLUDE_DIR}/openvdb/tools/MeshToVolume.h" OPENVDB_MESH_TO_VOLUME_CONTENT)
    string(REPLACE "#include \"ChangeBackground.h\"" "#include \"${OPENVDB_INCLUDE_DIR}/openvdb/tools/ChangeBackground.h\"" OPENVDB_MESH_TO_VOLUME_CONTENT "${OPENVDB_MESH_TO_VOLUME_CONTENT}")
    string(REPLACE "#include \"Prune.h\"" "#include \"${OPENVDB_INCLUDE_DIR}/openvdb/tools/Prune.h\"" OPENVDB_MESH_TO_VOLUME_CONTENT "${OPENVDB_MESH_TO_VOLUME_CONTENT}")
    string(REPLACE "#include \"SignedFloodFill.h\"" "#include \"${OPENVDB_INCLUDE_DIR}/openvdb/tools/SignedFloodFill.h\"" OPENVDB_MESH_TO_VOLUME_CONTENT "${OPENVDB_MESH_TO_VOLUME_CONTENT}")
    string(REPLACE "tbb::task::self().cancel_group_execution();" "(void)0;" OPENVDB_MESH_TO_VOLUME_CONTENT "${OPENVDB_MESH_TO_VOLUME_CONTENT}")
    file(WRITE "${OPENVDB_COMPAT_MESH_TO_VOLUME_HEADER}" "${OPENVDB_MESH_TO_VOLUME_CONTENT}")

    file(READ "${OPENVDB_INCLUDE_DIR}/openvdb/tools/PointIndexGrid.h" OPENVDB_POINT_INDEX_GRID_CONTENT)
    string(REPLACE "#include \"PointPartitioner.h\"" "#include \"${OPENVDB_INCLUDE_DIR}/openvdb/tools/PointPartitioner.h\"" OPENVDB_POINT_INDEX_GRID_CONTENT "${OPENVDB_POINT_INDEX_GRID_CONTENT}")
    string(REPLACE "tbb::task::self().cancel_group_execution();" "(void)0;" OPENVDB_POINT_INDEX_GRID_CONTENT "${OPENVDB_POINT_INDEX_GRID_CONTENT}")
    file(WRITE "${OPENVDB_COMPAT_POINT_INDEX_GRID_HEADER}" "${OPENVDB_POINT_INDEX_GRID_CONTENT}")

    file(WRITE "${OPENVDB_COMPAT_ROOT}/tbb/task.h"
"#pragma once\n"
"#include_next <tbb/task.h>\n"
"#include <tbb/task_group.h>\n"
"\n"
"namespace tbb {\n"
"namespace task {\n"
"\n"
"struct openvdb_task_self_compat_proxy {\n"
"    void cancel_group_execution() const {\n"
"        if (auto* ctx = current_context()) {\n"
"            ctx->cancel_group_execution();\n"
"        }\n"
"    }\n"
"};\n"
"\n"
"inline openvdb_task_self_compat_proxy self() {\n"
"    return openvdb_task_self_compat_proxy{};\n"
"}\n"
"\n"
"}  // namespace task\n"
"}  // namespace tbb\n")

    file(WRITE "${OPENVDB_COMPAT_ROOT}/tbb/task_scheduler_init.h"
"#pragma once\n"
"#include <tbb/task_arena.h>\n"
"\n"
"namespace tbb {\n"
"\n"
"class task_scheduler_init {\n"
"  public:\n"
"    static const int automatic = -1;\n"
"\n"
"    explicit task_scheduler_init(int = automatic) {}\n"
"    ~task_scheduler_init() = default;\n"
"\n"
"    void initialize(int = automatic) {}\n"
"    void terminate() {}\n"
"    bool is_active() const { return true; }\n"
"\n"
"    static int default_num_threads() {\n"
"        return task_arena::automatic;\n"
"    }\n"
"};\n"
"\n"
"}  // namespace tbb\n")

    target_include_directories(platform_backend BEFORE PUBLIC ${OPENVDB_COMPAT_ROOT})
    target_include_directories(platform_backend BEFORE PUBLIC ${OPENVDB_INCLUDE_DIR})
    target_include_directories(platform_backend PUBLIC ${NANOVDB_ROOT})
    target_link_libraries(platform_backend PUBLIC ${OPENVDB_LIBRARY} ${TBB_LIBRARY})
    target_compile_definitions(platform_backend PUBLIC SPCC_ENABLE_VDB=1 NANOVDB_USE_OPENVDB=1)

    message(STATUS "SPCC VDB enabled")
    message(STATUS " - OpenVDB include: ${OPENVDB_INCLUDE_DIR}")
    message(STATUS " - OpenVDB library: ${OPENVDB_LIBRARY}")
    message(STATUS " - TBB library    : ${TBB_LIBRARY}")
    message(STATUS " - OpenVDB runtime: ${SPCC_VDB_RUNTIME_DIR}")
    message(STATUS " - OpenVDB compat : ${OPENVDB_COMPAT_COORD_HEADER}")
    message(STATUS " - NanoVDB root   : ${NANOVDB_ROOT}")
endif()

# Minimal target executable
add_executable(baseline_drop_nsc apps/baseline_drop_nsc/main.cpp)
target_link_libraries(baseline_drop_nsc PRIVATE platform_backend Chrono::Chrono_core)
if(SPCC_ENABLE_VDB)
    target_link_options(baseline_drop_nsc PRIVATE "-Wl,-rpath,${SPCC_VDB_RUNTIME_DIR}")
endif()

add_executable(baseline_headon_sphere_nsc apps/baseline_headon_sphere_nsc/main.cpp)
target_link_libraries(baseline_headon_sphere_nsc PRIVATE platform_backend Chrono::Chrono_core)
if(SPCC_ENABLE_VDB)
    target_link_options(baseline_headon_sphere_nsc PRIVATE "-Wl,-rpath,${SPCC_VDB_RUNTIME_DIR}")
endif()

add_executable(baseline_rev_joint_clearance_nsc apps/baseline_rev_joint_clearance_nsc/main.cpp)
target_link_libraries(baseline_rev_joint_clearance_nsc PRIVATE platform_backend Chrono::Chrono_core)
if(SPCC_ENABLE_VDB)
    target_link_options(baseline_rev_joint_clearance_nsc PRIVATE "-Wl,-rpath,${SPCC_VDB_RUNTIME_DIR}")
endif()

# Cam case executable
add_executable(baseline_cam_nsc apps/baseline_cam_nsc/main.cpp)
target_link_libraries(baseline_cam_nsc PRIVATE platform_backend Chrono::Chrono_core)
if(SPCC_ENABLE_VDB)
    target_link_options(baseline_cam_nsc PRIVATE "-Wl,-rpath,${SPCC_VDB_RUNTIME_DIR}")
endif()

# Eccentric disk / roller follower case executable
add_executable(baseline_eccentric_roller_nsc apps/baseline_eccentric_roller_nsc/main.cpp)
target_link_libraries(baseline_eccentric_roller_nsc PRIVATE platform_backend Chrono::Chrono_core)
if(SPCC_ENABLE_VDB)
    target_link_options(baseline_eccentric_roller_nsc PRIVATE "-Wl,-rpath,${SPCC_VDB_RUNTIME_DIR}")
endif()

# Onset stress case executable
add_executable(baseline_onset_stress_nsc apps/baseline_onset_stress_nsc/main.cpp)
target_link_libraries(baseline_onset_stress_nsc PRIVATE platform_backend Chrono::Chrono_core)
if(SPCC_ENABLE_VDB)
    target_link_options(baseline_onset_stress_nsc PRIVATE "-Wl,-rpath,${SPCC_VDB_RUNTIME_DIR}")
endif()

# Onset stress case version B executable
add_executable(baseline_onset_stress_b_nsc apps/baseline_onset_stress_b_nsc/main.cpp)
target_link_libraries(baseline_onset_stress_b_nsc PRIVATE platform_backend Chrono::Chrono_core)
if(SPCC_ENABLE_VDB)
    target_link_options(baseline_onset_stress_b_nsc PRIVATE "-Wl,-rpath,${SPCC_VDB_RUNTIME_DIR}")
endif()

# Simple gear case executable
add_executable(baseline_simple_gear_nsc apps/baseline_simple_gear_nsc/main.cpp)
target_link_libraries(baseline_simple_gear_nsc PRIVATE platform_backend Chrono::Chrono_core)
if(SPCC_ENABLE_VDB)
    target_link_options(baseline_simple_gear_nsc PRIVATE "-Wl,-rpath,${SPCC_VDB_RUNTIME_DIR}")
endif()

add_executable(test_obj2sdf apps/test_obj2sdf/main.cpp)
target_link_libraries(test_obj2sdf PRIVATE platform_backend)
if(SPCC_ENABLE_VDB)
    target_link_options(test_obj2sdf PRIVATE "-Wl,-rpath,${SPCC_VDB_RUNTIME_DIR}")
endif()

add_executable(validate_compressed_contact apps/validate_compressed_contact/main.cpp)
target_link_libraries(validate_compressed_contact PRIVATE platform_backend)
if(SPCC_ENABLE_VDB)
    target_link_options(validate_compressed_contact PRIVATE "-Wl,-rpath,${SPCC_VDB_RUNTIME_DIR}")
endif()

# Offline tool to adaptively convert obj to vdb and extract surface
add_executable(mesh2vdb apps/mesh2vdb/main.cpp)
target_link_libraries(mesh2vdb PRIVATE platform_backend Half)
if(SPCC_ENABLE_VDB)
    target_link_options(mesh2vdb PRIVATE "-Wl,-rpath,${SPCC_VDB_RUNTIME_DIR}")
endif()

include(CTest)

if(BUILD_TESTING)
    set(GTEST_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/third_party/chrono/src/chrono_thirdparty/yaml-cpp/test/googletest-1.13.0/googletest)

    add_library(gtest_local STATIC
        ${GTEST_ROOT}/src/gtest-all.cc
        ${GTEST_ROOT}/src/gtest_main.cc
    )
    target_include_directories(gtest_local
        PUBLIC ${GTEST_ROOT}/include
        PRIVATE ${GTEST_ROOT}
    )
    find_package(Threads REQUIRED)
    target_link_libraries(gtest_local PUBLIC Threads::Threads)

    add_executable(spcc_tests
        tests/spcc/CompressedContactPipeline_test.cpp
    )

    if(SPCC_ENABLE_VDB)
        target_sources(spcc_tests PRIVATE tests/spcc/VDBSDFField_test.cpp)
        target_link_options(spcc_tests PRIVATE "-Wl,-rpath,${SPCC_VDB_RUNTIME_DIR}")
    endif()

    target_link_libraries(spcc_tests PRIVATE platform_backend gtest_local)
    add_test(NAME spcc_tests COMMAND spcc_tests)
endif()
