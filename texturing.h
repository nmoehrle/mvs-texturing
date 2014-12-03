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
#include "Arguments.h"

#include "RectangularBin.h"
#include "ObjModel.h"
#include "TextureView.h"
#include "UniGraph.h"
#include "TexturePatch.h"
#include "Tri.h"
#include "Timer.h"
#include "ProgressCounter.h"
#include "MRF.h"

struct VertexProjectionInfo {
    std::size_t texture_patch_id;
    math::Vec2f projection;
    std::vector<std::size_t> faces;
};

/**
  * Loads the mesh given by filename, removes duplicated faces and ensures normals (face and vertex).
  */
mve::TriangleMesh::Ptr
load_and_prepare_mesh(const std::string & filename);

/**
  * Generates TextureViews from the in_scene.
  * For each TextureView a valildity mask and the gradient magnitude image is calculated based on the images.
  */
void
generate_texture_views(Arguments const & conf, std::vector<TextureView> * texture_views);

/**
  * Builds up the meshes face adjacency graph using the vertex_infos
  */
void
build_adjacency_graph(UniGraph * graph, mve::TriangleMesh::ConstPtr mesh,
    mve::VertexInfoList::ConstPtr vertex_infos);

/**
  * Builds up the mrf graph for the mesh.
  * The visibility of faces within texture views is determined for the assignment of data costs.
  */
void
build_mrf(MRF * mrf, mve::TriangleMesh::ConstPtr mesh, std::vector<TextureView> const & texture_views,
    UniGraph const & graph, Arguments const & conf);

/**
  * Runs mrf optimization and save energies to conf.out_prefix + "_mrf_energies.csv".
  */
void
run_mrf_optimization(MRF * mrf, Arguments const & conf);

/**
  * Builds up the obj_model for the mesh constructing materials and texture atlases form the texture_patches
  */
void
build_obj_model(ObjModel * obj_model, std::vector<TexturePatch> const & texture_patches,
    mve::TriangleMesh::ConstPtr mesh);

/**
  * Generates texture patches using the graph to determine adjacent faces with the same label.
  * The created texture atlases are adjusted with the provided adjust_values.
  */
void generate_texture_patches(UniGraph const & graph, std::vector<TextureView> const & texture_views,
    mve::TriangleMesh::ConstPtr mesh, std::vector<std::vector<VertexProjectionInfo> > * vertex_projection_infos,
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
 * Partitions a mesh (represented by a face adjacency graph) into (more or less) convex, equally sized submeshes.
 * The submesh IDs are stored in the UniGraphs vertex labels.
 * This is useful if there are follow-up steps that can be parallelized over submeshes.
 * Uses something similar to k-means clustering with edge hops as distance metric.
 *
 * @param graph the mesh to be partitioned into submeshes
 * @param min_num_partitions the minimum number of partitions into which the mesh will be partitioned. The final number of submeshes may well exceed this if there are small, unconnected regions exist.
 */
void
partition_mesh(UniGraph * graph, std::size_t min_num_partitions);
