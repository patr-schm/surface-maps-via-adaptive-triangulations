cmake_minimum_required(VERSION 3.9)

# ===============================================
# Project name

# Set project name if not build within another project
if("${PROJECT_NAME}" STREQUAL "")
    project(SurfaceMaps)
endif()

# ===============================================
# Global settings

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Add -fPIC flag (required to build predicates.c)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Disable some compiler warnings
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-sign-compare")

# Turn some compiler warnings into errors
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Werror=return-type")

# Set glow bin dir
set(GLOW_BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/bin/glow)

# Set output dir for meshes, screenshots, etc...
set(ISM_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/output)

# Disable glfw additionals
option(GLFW_BUILD_EXAMPLES "" OFF)
option(GLFW_BUILD_TESTS "" OFF)
option(GLFW_BUILD_DOCS "" OFF)
option(GLFW_INSTALL "" OFF)

# ===============================================
# Find packages

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(OpenMP REQUIRED)

# ===============================================
# Add submodules

add_subdirectory(libs/OpenMesh)

# Supress an OpenMesh warning
if (${QT_QMAKE_EXECUTABLE})
endif()

if (NOT TARGET glfw)
    add_subdirectory(libs/glfw)
endif()

if (NOT TARGET clean-core)
    add_subdirectory(libs/clean-core)
endif()

if (NOT TARGET typed-geometry)
    add_subdirectory(libs/typed-geometry)
endif()

if (NOT TARGET polymesh)
    add_subdirectory(libs/polymesh)
endif()

if (NOT TARGET glow)
    add_subdirectory(libs/glow)
endif()

if (NOT TARGET imgui)
    add_subdirectory(libs/imgui-lean)
endif()

if (NOT TARGET glow-extras)
    add_subdirectory(libs/glow-extras)
endif()

if (NOT TARGET Eigen3::Eigen)
    add_subdirectory(libs/eigen-lean)
endif()

add_subdirectory(libs/TinyAD)

# ===============================================
# Target SurfaceMaps

file(GLOB_RECURSE ISM_HEADERS src/*.hh libs/ACG_BSP/*.hh libs/predicates/*.h)
file(GLOB_RECURSE ISM_SOURCES src/*.cc libs/ACG_BSP/*.cc libs/predicates/*.c)
add_library(SurfaceMaps ${ISM_SOURCES} ${ISM_HEADERS})
target_compile_features(SurfaceMaps PUBLIC cxx_std_17)
target_compile_definitions(SurfaceMaps PUBLIC SOURCE_PATH_STR="${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_definitions(SurfaceMaps PUBLIC DATA_PATH_STR="${CMAKE_CURRENT_SOURCE_DIR}/data")
target_compile_definitions(SurfaceMaps PUBLIC OUTPUT_PATH_STR="${ISM_OUTPUT_DIR}")
target_include_directories(SurfaceMaps PUBLIC
    libs/libigl/include
    ${OPENMESH_INCLUDE_DIRS}
    libs/predicates
    libs/ACG_BSP
    src
)
target_link_libraries(SurfaceMaps PUBLIC
    eigen
    OpenMP::OpenMP_CXX
    OpenMeshCore
    clean-core
    typed-geometry
    polymesh
    glfw
    glow
    imgui
    glow-extras
    TinyAD
    stdc++fs
)

# ===============================================
# Target apps (Executables)

file(GLOB_RECURSE ISM_APPS apps/*.cc)
foreach(main ${ISM_APPS})
    get_filename_component(name ${main} NAME_WE)
    get_filename_component(dir ${main} DIRECTORY)
    get_filename_component(dir_name ${dir} NAME)
    message(STATUS "Adding app ${name}")

    add_executable(${name} ${main})
    target_compile_definitions(${name} PRIVATE APP_EXECUTABLES)
    add_dependencies(${name} SurfaceMaps)
    target_include_directories(${name} PUBLIC
        "src"
    )
    target_link_libraries(${name} PUBLIC
        SurfaceMaps
    )
endforeach()
