/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Janis Born, Joe Jakobi, Patrick Schmidt
 */
#pragma once

#include <SurfaceMaps/Types.hh>

namespace SurfaceMaps
{

struct ProgressiveMesh
{
    struct Record
    {
        VH vh0; // Not necessary for reconstruction, but used to check consistency of decimation order
        VH vh1;
        VH vhl;
        VH vhr;

        Vec3d p0; // The original rest state position of v0
    };

    // A log of all decimations that were applied to the input mesh. If you want
    // to reconstruct the mesh, you need to iterate over this sequence in
    // reverse order. All indices in this log refer to indices of the
    // progressive mesh.
    std::vector<Record> log;

    // Index maps between vertex indices of the original input mesh and the
    // (reconstructed) progressive mesh.
    std::vector<int> orig_to_prog_idx;
    std::vector<int> prog_to_orig_idx;
};

enum Priority
{
    RANDOM,
    DIHEDRAL_ANGLE_CHANGE,
    UNIFORM_DECIMATION,
    ERROR_QUADRICS,
    ANGLE_DEFECT,
    ANGLE_DEFECT_CHANGE_WITH_MEMORY,
    ANGLE_DEFECT_CHANGE_WITHOUT_MEMORY,
    MAX_MINIMAL_FACE_AREA
};

static std::vector<std::string> PRIORITY_NAMES
{
    "RANDOM",
    "DIHEDRAL_ANGLE_CHANGE",
    "UNIFORM_DECIMATION",
    "ERROR_QUADRICS",
    "ANGLE_DEFECT",
    "ANGLE_DEFECT_CHANGE_WITH_MEMORY",
    "ANGLE_DEFECT_CHANGE_WITHOUT_MEMORY",
    "MAX_MINIMAL_FACE_AREA"
};

struct ProgressiveMeshOptions
{
    Priority priority = ERROR_QUADRICS;
    bool create_independent_sets = true;
    bool weight_with_minimal_face_area = true;
    double weight_coeff = 1e-6;
    double zero_area_delay_thres = 1e-10;
};

//! Decimates the given mesh and returns a decimation history
ProgressiveMesh compute_progressive_mesh(
        TriMesh& _mesh,
        const int _n_target_vertices,
        const ProgressiveMeshOptions& _options,
        std::vector<uint>& _indep_set_sizes,
        std::function<bool(HEH)> _is_legal_collapse = [](auto){return true; },
        std::function<void(HEH)> _pre_collapse = [](auto){},
        std::function<void(VH, VH, VH)> _post_collapse = [](auto, auto, auto){});

ProgressiveMesh compute_progressive_mesh(
        TriMesh& _mesh,
        const int _n_target_vertices,
        const ProgressiveMeshOptions& _options,
        std::function<bool(HEH)> _is_legal_collapse = [](auto){return true; },
        std::function<void(HEH)> _pre_collapse = [](auto){},
        std::function<void(VH, VH, VH)> _post_collapse = [](auto, auto, auto){});

//! Undoes one halfedge collapse (represented by rec). Tries to reconstruct a
//! valid (flip free) parameter domain position (stored in param) for the newly
//! inserted vertex.
//! Undoes one halfedge collapse (represented by rec).
void undo_decimation_step(
    TriMesh& _mesh,
    const ProgressiveMesh::Record& _rec);

void reconstruct_progressive_mesh(
    TriMesh& _mesh,
    const ProgressiveMesh& _pm);

}
