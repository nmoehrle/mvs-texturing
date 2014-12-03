#include "Arguments.h"

Arguments parse_args(int argc, char **argv) {
    util::Arguments args;
    args.set_exit_on_error(true);
    args.set_nonopt_maxnum(3);
    args.set_nonopt_minnum(3);
    args.set_helptext_indent(34);
    args.set_description("Textures a mesh given images in form of a 3D scene.");
    args.set_usage("Usage: " + std::string(argv[0]) + " [options] IN_SCENE IN_MESH OUT_MESH_PREFIX"
        "\n\nIN_SCENE := (SCENE_FOLDER | BUNDLE_FILE | MVE_SCENE::EMBEDDING)"
        "\n\nSCENE_FOLDER:"
        "\nWithin a scene folder a .cam file has to be given for each image."
        "\nA .cam file is structured as follows:"
        "\n    tx ty tz R00 R01 R02 R10 R11 R12 R20 R21 R22"
        "\n    f d0 d1 paspect ppx ppy"
        "\nFirst line: Extrinsics - translation vector and rotation matrix"
        "\nSecond line: Intrinsics - focal length, distortion coefficients, pixel aspect ratio and principal point"
        "\nThe focal length is the distance between camera center and image plane normalized by dividing with the larger image dimension."
        "\nFor non zero distortion coefficients the image will be undistorted prior to the texturing process."
        " If only d0 is non zero the Noah Snavely's distortion model is assumed otherwise the distortion model of VSFM is assumed."
        "\nThe pixel aspect ratio is usually 1 or close to 1. If your SfM system doesn't output it, but outputs a different focal length in x and y direction, you have to encode this here."
        "\nThe principal point has to be given in unit dimensions (e.g. 0.5 0.5)."
        "\n\nBUNDLE_FILE:"
        "\nCurrently only NVM bundle files (from VisualSFM, http://ccwu.me/vsfm/) are supported."
        "\nSince the bundle file contains relative paths to the images please make sure you did not move them (relative to the bundle) or rename them after the bundling process."
        "\n\nMVE_SCENE::EMBEDDING:"
        "\nThis is the scene representation we use in our research group: http://www.gris.tu-darmstadt.de/projects/multiview-environment/."
        "\n\nIN_MESH:"
        "\nThe mesh that you want to texture and which needs to be in the same coordinate frame as the camera parameters. You can reconstruct one, e.g. with CMVS: http://www.di.ens.fr/cmvs/"
        "\n\nOUT_MESH_PREFIX:"
        "\nA path and name for the output files, e.g. <path>/<to>/my_textured_mesh"
        "\nDon't append an obj extension. The application does that itself because it outputs multiple files (mesh, material file, texture files)."
        "\n");
    args.add_option('D',"data_cost_file", true,
        "Skip calculation of data costs and use the ones provided in the given file");
    args.add_option('L',"labeling_file", true,
        "Skip view selection and use the labeling provided in the given file");
    args.add_option('d',"data_term", true,
        "Data term: {area, gmi} [gmi]");
    args.add_option('s',"smoothness_term", true,
        "Smoothness term: {potts, edi} [potts]");
    args.add_option('o',"outlier_removal", true,
        "Photometric outlier (pedestrians etc.) removal method: {none, gauss_clamping, gauss_damping} [none]");
    args.add_option('v',"view_selection_model", false,
        "Write out view selection model [false]");
    args.add_option('\0',"skip_global_seam_leveling", false,
        "Skip global seam leveling [false]");
    args.add_option('\0',"skip_geometric_visibility_test", false,
        "Skip geometric visibility test based on ray intersection [false]");
    args.add_option('\0',"skip_local_seam_leveling", false,
        "Skip local seam leveling (Poisson editing) [false]");
    args.parse(argc, argv);

    Arguments conf;
    conf.in_scene = args.get_nth_nonopt(0);
    conf.in_mesh = args.get_nth_nonopt(1);
    conf.out_prefix = args.get_nth_nonopt(2);

    /* Set defaults for optional arguments. */
    conf.data_cost_file = "";
    conf.labeling_file = "";
    conf.data_term = GMI;
    conf.smoothness_term = POTTS;
    conf.outlier_removal = NONE;
    conf.write_view_selection_model = false;
    conf.write_data_term_histograms = false;
    conf.write_mrf_energies = false;
    conf.geometric_visibility_test = true;
    conf.global_seam_leveling = true;
    conf.local_seam_leveling = true;

    /* Handle optional arguments. */
    for (util::ArgResult const* i = args.next_option();
         i != 0; i = args.next_option()) {
        switch (i->opt->sopt) {
        case 'v':
            conf.write_view_selection_model = true;
        break;
        case 'D':
            conf.data_cost_file = i->arg;
        break;
        case 'L':
            conf.labeling_file = i->arg;
        break;
        case 'd':
            conf.data_term = parse_data_term(i->arg);
        break;
        case 's':
            conf.smoothness_term = parse_smoothness_term(i->arg);
        break;
        case 'o':
            conf.outlier_removal = parse_outlier_removal(i->arg);
        break;
        case '\0':
            if (i->opt->lopt == "skip_geometric_visibility_test") {
                conf.geometric_visibility_test = false;
            } else if (i->opt->lopt == "skip_global_seam_leveling") {
                conf.global_seam_leveling = false;
            } else if (i->opt->lopt == "skip_local_seam_leveling") {
                conf.local_seam_leveling = false;
            } else {
                throw std::invalid_argument("Invalid long option");
            }
        break;
        default:
            throw std::invalid_argument("Invalid short option");
        }
    }

    return conf;
}

std::string
bool_to_string(bool b){
    return b ? "True" : "False";
}

std::string
Arguments::to_string(){
    std::stringstream out;
    out << "Input scene: \t" << in_scene << std::endl
        << "Input mesh: \t" << in_mesh << std::endl
        << "Output prefix: \t" << out_prefix << std::endl
        << "Datacost file: \t" << data_cost_file << std::endl
        << "Labeling file: \t" << labeling_file << std::endl
        << "Data term: \t" << DataTermStrings[data_term] << std::endl
        << "Smoothness term: \t" << SmoothnessTermStrings[smoothness_term] << std::endl
        << "Outlier removal method: \t" << OutlierRemovalStrings[outlier_removal] << std::endl
        << "Write view selection model: \t" << bool_to_string(write_view_selection_model) << std::endl
        << "Apply global seam leveling: \t" << bool_to_string(global_seam_leveling) << std::endl
        << "Apply local seam leveling: \t" << bool_to_string(local_seam_leveling) << std::endl;

    return out.str();
}
