#include <iostream>
#include <fstream>
#include <vector>
#include <set>

#include "math/matrix.h"

#include "util/arguments.h"
#include "util/timer.h"
#include "util/file_system.h"

#include "mve/scene.h"

#include "texturing.h"
#include "debug.h"

int main(int argc, char **argv) {

#ifdef RESEARCH
    std::cout << "******************************************************************************" << std::endl
              << " Due to use of the -DRESEARCH=ON compile option, this program is licensed "     << std::endl
              << " for research purposes only. Please pay special attention to the gco license."  << std::endl
              << "******************************************************************************" << std::endl;
#endif

    Timer timer;
    timer.measure("Start");

    util::WallTimer wtimer;

    Arguments conf;
    try {
        conf = parse_args(argc, argv);
    } catch (std::invalid_argument & ia) {
        std::cerr << ia.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (!util::fs::dir_exists(util::fs::dirname(conf.out_prefix).c_str())) {
        std::cerr << "Destination directory does not exist!" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::cout << "Load and prepare mesh: " << std::endl;
    mve::TriangleMesh::ConstPtr mesh = load_and_prepare_mesh(conf.in_mesh);
    mve::VertexInfoList::ConstPtr vertex_infos = mve::VertexInfoList::create(mesh);

    std::size_t const num_faces = mesh->get_faces().size() / 3;

    std::cout << "Generating texture views: " << std::endl;
    std::vector<TextureView> texture_views;
    generate_texture_views(conf, &texture_views);

    write_string_to_file(conf.out_prefix + ".conf", conf.to_string());
    timer.measure("Loading");
    wtimer.reset();

    std::cout << "Building adjacency graph: " << std::endl;
    UniGraph graph(num_faces);
    build_adjacency_graph(&graph, mesh, vertex_infos);

    if (conf.labeling_file.empty()) {
        std::cout << "Building MRF graph:" << std::endl;

        /* Each TextureView is a label and label 0 is undefined */
#ifdef RESEARCH
        GCOWrapper mrf(num_faces, texture_views.size() + 1);
#else
        LBPSolver mrf(num_faces, texture_views.size() + 1);
#endif
        build_mrf(&mrf, mesh, texture_views, graph, conf);
        timer.measure("Calculating data costs");

        std::cout << "Running MRF optimization:" << std::endl;
        run_mrf_optimization(&mrf, conf);
        timer.measure("Running MRF optimization");

        /* Extract resulting labeling from MRF. */
        for (std::size_t i = 0; i < graph.num_nodes(); ++i) {
            int label = mrf.what_label(static_cast<int>(i));
            assert(0 <= label && label < texture_views.size() + 1);
            graph.set_label(i, static_cast<std::size_t>(label));
        }

        /* Fix holes due to geometric inaccuracies. */
        std::vector<std::vector<std::size_t> > subgraphs;
        graph.get_subgraphs(0, &subgraphs);
        #pragma omp parallel for
        for (std::size_t i = 0; i < subgraphs.size(); ++i) {
            /* Only fix small holes, */
            const std::size_t max_hole_size = 10;
            if (subgraphs[i].size() < max_hole_size) {
                /* which are completely contained within another patch. */
                std::set<std::size_t> adj_labels;
                graph.get_adj_labels(subgraphs[i], &adj_labels);
                adj_labels.erase(0);
                if (adj_labels.size() == 1 && graph.get_min_conn(subgraphs[i]) == 3) {
                    const std::size_t new_label = *(adj_labels.begin());
                    for (std::size_t node : subgraphs[i]) {
                        graph.set_label(node, new_label);
                        //TODO check if projection is valid... (very likely)
                    }
                }
            }
        }

        /* Write labeling to file. */
        if (conf.write_intermediate_results) {
            std::vector<std::size_t> labeling(graph.num_nodes());
            for (std::size_t i = 0; i < graph.num_nodes(); ++i) {
                labeling[i] = graph.get_label(i);
            }
            vector_to_file(conf.out_prefix + "_labeling.vec", labeling);
        }

        /* Release superfluous embeddings. */
        for (TextureView & texture_view : texture_views) {
            texture_view.release_validity_mask();
            texture_view.release_gradient_magnitude();
        }

    } else {
        std::cout << "Loading labeling from file... " << std::flush;

        /* Load labeling to file. */
        std::vector<std::size_t> labeling = vector_from_file<std::size_t>(conf.labeling_file);
        if (labeling.size() != graph.num_nodes()) {
            std::cerr << "Wrong labeling file for this mesh/scene combination... aborting!" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        /* Transfer labeling to graph. */
        for (std::size_t i = 0; i < labeling.size(); ++i) {
            const std::size_t label = labeling[i];
            if (label > texture_views.size()){
                std::cerr << "Wrong labeling file for this mesh/scene combination... aborting!" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            graph.set_label(i, label);
        }

        std::cout << "done." << std::endl;
    }

    /* Create texture patches and adjust them. */
    std::vector<TexturePatch> texture_patches;
    {
        std::vector<std::vector<VertexProjectionInfo> > vertex_projection_infos;
        std::cout << "Generating texture patches:" << std::endl;
        generate_texture_patches(graph, texture_views, mesh, &vertex_projection_infos, &texture_patches);
        for (TextureView & texture_view : texture_views) {
            texture_view.release_image();
        }

        if (conf.global_seam_leveling) {
            std::cout << "Running global seam leveling:" << std::endl;
            global_seam_leveling(graph, mesh, vertex_infos, vertex_projection_infos, &texture_patches);
            timer.measure("Running global seam leveling");
        } else {
            ProgressCounter texture_patch_counter("Calculating validity masks for texture patches", texture_patches.size());
            #pragma omp parallel for schedule(dynamic)
            for (std::size_t i = 0; i < texture_patches.size(); ++i) {
                texture_patch_counter.progress<SIMPLE>();
                TexturePatch * texture_patch = &texture_patches[i];
                std::vector<math::Vec3f> patch_adjust_values(texture_patch->get_faces().size() * 3, math::Vec3f(0.0f));
                texture_patch->adjust_colors(patch_adjust_values);
                texture_patch_counter.inc();
            }
            timer.measure("Calculating texture patch validity masks");
        }

        if (conf.local_seam_leveling) {
            std::cout << "Running local seam leveling:" << std::endl;
            local_seam_leveling(graph, mesh, vertex_projection_infos, &texture_patches);
        }
        timer.measure("Running local seam leveling");
    }

    /* Create and write out obj model. */
    {
        ObjModel obj_model;
        std::cout << "Building objmodel:" << std::endl;

        build_obj_model(&obj_model, texture_patches, mesh);
        timer.measure("Building OBJ model");

        std::cout << "\tWriting files... " << std::flush;
        obj_model.save_to_files(conf.out_prefix);
        std::cout << "done." << std::endl;
        timer.measure("Saving");
    }

    std::cout << "Whole texturing procedure took: " << wtimer.get_elapsed_sec() << "s" << std::endl;
    timer.measure("Total");
    if (conf.write_timings) {
        timer.write_to_file(conf.out_prefix + "_timings.csv");
    }

    if (conf.write_view_selection_model) {
        std::cout << "Generating debug texture patches:" << std::endl;
        {
            texture_patches.clear();
            generate_debug_embeddings(&texture_views);
            std::vector<std::vector<VertexProjectionInfo> > vertex_projection_infos; // Will only be written
            generate_texture_patches(graph, texture_views, mesh, &vertex_projection_infos, &texture_patches);
        }

        std::cout << "Building debug objmodel:" << std::endl;
        {
            ObjModel obj_model;
            build_obj_model(&obj_model, texture_patches, mesh);
            std::cout << "\tWriting files... " << std::flush;
            obj_model.save_to_files(conf.out_prefix + "_view_selection");
            std::cout << "done." << std::endl;
        }
    }
}
