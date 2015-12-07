/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef TEX_SEAMLEVELING_HEADER
#define TEX_SEAMLEVELING_HEADER

#include <vector>

#include <mve/mesh.h>

#include "defines.h"
#include "uni_graph.h"

TEX_NAMESPACE_BEGIN

struct VertexProjectionInfo {
    std::size_t texture_patch_id;
    math::Vec2f projection;
    std::vector<std::size_t> faces;

    bool operator<(VertexProjectionInfo const & other) const {
        return texture_patch_id < other.texture_patch_id;
    }
};

struct EdgeProjectionInfo {
    std::size_t texture_patch_id;
    math::Vec2f p1;
    math::Vec2f p2;

    bool operator<(EdgeProjectionInfo const & other) const {
        return texture_patch_id < other.texture_patch_id;
    }
};

struct MeshEdge {
    std::size_t v1;
    std::size_t v2;
};

void
find_seam_edges(UniGraph const & graph, mve::TriangleMesh::ConstPtr mesh,
    std::vector<MeshEdge> * seam_edges);

void
find_mesh_edge_projections(
    std::vector<std::vector<VertexProjectionInfo> > const & vertex_projection_infos,
    MeshEdge mesh_edge, std::vector<EdgeProjectionInfo> * projected_edge_infos);

TEX_NAMESPACE_END

#endif /* TEX_SEAMLEVELING_HEADER */
