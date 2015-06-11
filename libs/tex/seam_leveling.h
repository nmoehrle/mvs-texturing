#pragma once
#include <vector>

#include <mve/mesh.h>

#include "UniGraph.h"

struct VertexProjectionInfo {
    std::size_t texture_patch_id;
    math::Vec2f projection;
    std::vector<std::size_t> faces;
};

struct ProjectedEdgeInfo {
    std::size_t texture_patch_id;
    math::Vec2f p1;
    math::Vec2f p2;

    bool operator<(ProjectedEdgeInfo const & projected_edge_info) const {
        return texture_patch_id < projected_edge_info.texture_patch_id;
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
    MeshEdge mesh_edge, std::vector<ProjectedEdgeInfo> * projected_edge_infos);
