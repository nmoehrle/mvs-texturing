#pragma once

#include <map>
#include <list>
#include <vector>
#include <iostream>
#include <sstream>

#include "mve/image.h"
#include "mve/scene.h"
#include "mve/mesh.h"
#include "mve/mesh_info.h"
#include "mve/image_io.h"
#include "mve/image_drawing.h"

#include "util/timer.h"
#include "math/accum.h"

#include "mrf/Graph.h"

#include "coldet.h"

#include "util.h"

#include "defines.h"
#include "Settings.h"
#include "RectangularBin.h"
#include "ObjModel.h"
#include "TextureView.h"
#include "UniGraph.h"
#include "TexturePatch.h"
#include "Tri.h"
#include "Timer.h"
#include "ProgressCounter.h"
#include "SparseTable.h"

#include "seam_leveling.h"

typedef SparseTable<std::uint32_t, std::uint16_t, float> ST;

TEX_NAMESPACE_BEGIN

typedef std::vector<TextureView> TextureViews;
typedef std::vector<TexturePatch> TexturePatches;
typedef ObjModel Model;
typedef UniGraph Graph;
typedef ST DataCosts;
typedef std::vector<std::vector<VertexProjectionInfo> > VertexProjectionInfos;
/**
  * prepares the mesh for texturing
  *  -removes duplicated faces
  *  -ensures normals (face and vertex)
  */
void
prepare_mesh(mve::VertexInfoList::Ptr vertex_infos, mve::TriangleMesh::Ptr mesh);

/**
  * Generates TextureViews from the in_scene.
  */
void
generate_texture_views(std::string in_scene, std::vector<TextureView> * texture_views);

/**
  * Builds up the meshes face adjacency graph using the vertex_infos
  */
void
build_adjacency_graph(mve::TriangleMesh::ConstPtr mesh, mve::VertexInfoList::ConstPtr vertex_infos, UniGraph * graph);

/**
 * Calculates the data costs for each face and texture view combination,
 * if the face is visible within the texture view.
 */
void
calculate_data_costs(mve::TriangleMesh::ConstPtr mesh, std::vector<TextureView> * texture_views,
    Settings const & settings, ST * data_costs);

/**
 * Runs the view selection procedure and saves the labeling in the graph
 */
void
view_selection(ST const & data_costs, UniGraph * graph, Settings const & settings);

/**
  * Generates texture patches using the graph to determine adjacent faces with the same label.
  */
void generate_texture_patches(UniGraph const & graph,
    mve::TriangleMesh::ConstPtr mesh,
    mve::VertexInfoList::ConstPtr vertex_infos,
    std::vector<TextureView> * texture_views,
    VertexProjectionInfos * vertex_projection_infos,
    std::vector<TexturePatch> * texture_patches);

/**
  * Runs the seam leveling procedure proposed by Ivanov and Lempitsky
  * [<A HREF="https://www.google.de/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&sqi=2&ved=0CC8QFjAA&url=http%3A%2F%2Fwww.robots.ox.ac.uk%2F~vilem%2FSeamlessMosaicing.pdf&ei=_ZbvUvSZIaPa4ASi7IGAAg&usg=AFQjCNGd4x5HnMMR68Sn2V5dPgmqJWErCA&sig2=4j47bXgovw-uks9LBGl_sA">Seamless mosaicing of image-based texture maps</A>]
  */
void
global_seam_leveling(UniGraph const & graph, mve::TriangleMesh::ConstPtr mesh,
    mve::VertexInfoList::ConstPtr vertex_infos,
    VertexProjectionInfos const & vertex_projection_infos,
    std::vector<TexturePatch> * texture_patches);

void
local_seam_leveling(UniGraph const & graph, mve::TriangleMesh::ConstPtr mesh,
    VertexProjectionInfos const & vertex_projection_infos,
    std::vector<TexturePatch> * texture_patches);

/**
  * Builds up an model for the mesh by constructing materials and
  * texture atlases form the texture_patches
  */
void
build_model(mve::TriangleMesh::ConstPtr mesh,
    std::vector<TexturePatch> const & texture_patches, Model * model);

TEX_NAMESPACE_END
