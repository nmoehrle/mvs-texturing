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

#include "coldet.h"

#include "util.h"

#include "RectangularBin.h"
#include "ObjModel.h"
#include "TextureView.h"
#include "UniGraph.h"
#include "TexturePatch.h"
#include "Tri.h"
#include "Timer.h"
#include "ProgressCounter.h"
#include "MRF.h"
#include "SparseTable.h"

#include "Settings.h"

typedef SparseTable<std::uint32_t, std::uint16_t, int> ST;

struct VertexProjectionInfo {
    std::size_t texture_patch_id;
    math::Vec2f projection;
    std::vector<std::size_t> faces;
};


/**
  * Loads the mesh given by filename,
  * removes duplicated faces and ensures normals (face and vertex).
  */
mve::TriangleMesh::Ptr
load_and_prepare_mesh(const std::string & filename);

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
ST
calculate_data_costs(mve::TriangleMesh::ConstPtr mesh, std::vector<TextureView> & texture_views,
    Settings const & settings);

/**
  * Builds up the mrf graph for the mesh.
  */
void
build_mrf(UniGraph graph, ST const & data_costs, MRF * mrf, Settings const & settings);

/**
  * Runs mrf optimization.
  */
void
run_mrf_optimization(MRF * mrf);

/**
  * Generates texture patches using the graph to determine adjacent faces with the same label.
  */
void generate_texture_patches(UniGraph const & graph,
    std::vector<TextureView> const & texture_views,
    mve::TriangleMesh::ConstPtr mesh,
    std::vector<std::vector<VertexProjectionInfo> > * vertex_projection_infos,
    std::vector<TexturePatch> * texture_patches);

/**
  * Runs the seam leveling procedure proposed by Ivanov and Lempitsky
  * [<A HREF="https://www.google.de/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&sqi=2&ved=0CC8QFjAA&url=http%3A%2F%2Fwww.robots.ox.ac.uk%2F~vilem%2FSeamlessMosaicing.pdf&ei=_ZbvUvSZIaPa4ASi7IGAAg&usg=AFQjCNGd4x5HnMMR68Sn2V5dPgmqJWErCA&sig2=4j47bXgovw-uks9LBGl_sA">Seamless mosaicing of image-based texture maps</A>]
  */
void
global_seam_leveling(UniGraph const & graph, mve::TriangleMesh::ConstPtr mesh,
    mve::VertexInfoList::ConstPtr vertex_infos,
    std::vector<std::vector<VertexProjectionInfo> > const & vertex_projection_infos,
    std::vector<TexturePatch> * texture_patches);

void
local_seam_leveling(UniGraph const & graph, mve::TriangleMesh::ConstPtr mesh,
    std::vector<std::vector<VertexProjectionInfo> > const & vertex_projection_infos,
    std::vector<TexturePatch> * texture_patches);

/**
  * Builds up an obj_model for the mesh by constructing materials and
  * texture atlases form the texture_patches
  */
void
build_obj_model(mve::TriangleMesh::ConstPtr mesh,
    std::vector<TexturePatch> const & texture_patches, ObjModel * obj_model);

