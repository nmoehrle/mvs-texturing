# HomeeAI's MVS-Texturing

This repository is a fork of the [official mvs-texturing
repository](https://github.com/nmoehrle/mvs-texturing).

---

## Prerequisites

To build and run `mvs-texturing`, ensure the following dependencies are installed:

- **CMake** >= 3.28  
- **Git**  
- **Make**  
- **GCC/G++** >= 12  
- **Libraries**: `libpng`, `libjpeg`, `libtiff`

Additional dependencies are automatically downloaded and built by the system via
`elibs/CMakeLists.txt`:

- [rayint](https://github.com/nmoehrle/rayint)  
- [Eigen](http://eigen.tuxfamily.org)  
- [Multi-View Environment (MVE)](http://www.gcc.tu-darmstadt.de/home/proj/mve)  
- [mapMAP](http://www.gcc.tu-darmstadt.de/home/proj/mapmap)  

### Install Legacy oneTBB

The current code requires **oneTBB 2019_U9**, which needs to be downloaded and
built manually:

```shell
wget https://github.com/oneapi-src/oneTBB/archive/refs/tags/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cd oneTBB-2019_U9
make -j
```

## Build MVS-Texturing `mvs-texturing`

Use the following command to build the `mvs-texturin`g project:

```
mkdir -p build && cd build && cmake .. && make -j; cd ..
```

**Note:** If there are issues related to oneTBB, pass the correct paths to the
TBB headers and libraries using CMake flags:

```shell
cmake .. -DTBB_INCLUDE_DIRS=../oneTBB-2019_U9/include \
         -DTBB_LIBRARIES=../oneTBB-2019_U9/build/linux_intel64_gcc_*_release
```

## Prepare 3D Model and Camera Parameters

To prepare your input data, follow the preprocessing script
`scripts/mvs_texturing_preprocess_camera_parameters.py` in
the[`inpaint-3d-struct`
repository](https://github.com/homee-ai/inpaint-3d-struct).

Essentially, `mvs-texturing` only generates good results if input data has the
correct format:

- **Mesh Quality:** Ensure the input mesh has appropriate triangle sizes. Large
  triangles may result in incomplete texturing if no single image fully covers a
  triangle.

- **Camera Parameters:** The .cam file must follow this format:

```
tx ty tz r00 r01 r02 r10 r11 r12 r20 r21 r22  # extrinsic
focal_length distort_1 distort_2 pixel_aspect_ratio principal_point_x principal_point_y  # intrinsic
```

Please read [mve/libs/mve/camera.h](https://github.com/simonfuhrmann/mve/blob/master/libs/mve/camera.h) for more information.


## Running MVS-Texturing

After `mvs-texturing` is built and input data is prepared, we can start texturizing meshes. The
following command

```shell
mkdir -p output \
&& ./build/apps/texrecon/texrecon \
--tone_mapping=gamma \
--outlier_removal=gauss_damping \
--view_selection_model \
--keep_unseen_faces \
./input/img_cam/ \
./input/scene.ply \
./output/textured_mesh
```

takes the images, camera parameters, and a `scene.ply` file as input and produces
the results in a destination path prefixed with `./output/textured_mesh`, as
shown below: 

```
textured_mesh.conf
textured_mesh.mtl
textured_mesh.obj
textured_mesh_data_costs.spt
textured_mesh_labeling.vec
textured_mesh_material0000_map_Kd.png
textured_mesh_material0001_map_Kd.png
textured_mesh_view_selection.mtl
textured_mesh_view_selection.obj
textured_mesh_view_selection_material0000_map_Kd.png
textured_mesh_view_selection_material0001_map_Kd.png
```
