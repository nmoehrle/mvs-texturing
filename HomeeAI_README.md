# HomeeAI's MVS-Texturing

This repository is forked from [the official mvs-texturing](https://github.com/nmoehrle/mvs-texturing).

## Dependencies

`mvs-texturing` depends on the following prerequisites:
- `cmake>=3.28`
- git
- make
- `gcc>=12` and `g++>=12`
- `libpng`, `libjpb` and `libtiff`

Furthermore the build system automatically downloads and compiles the following
dependencies in `elibs/CMakeLists.txt`:

- rayint: https://github.com/nmoehrle/rayint
- Eigen: http://eigen.tuxfamily.org
- Multi-View Environment: http://www.gcc.tu-darmstadt.de/home/proj/mve
- mapMAP: http://www.gcc.tu-darmstadt.de/home/proj/mapmap

### Download and Build the Old `oneTBB`

The current code depends on the old version of `oneTBB`. We need to download and
compile `oneTBB` manually. Please follow the instruction below:

```shell
wget https://github.com/oneapi-src/oneTBB/archive/refs/tags/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cd oneTBB-2019_U9
make -j
```

## Build `mvs-texturing`

```
mkdir -p build && cd build && cmake .. && make -j; cd ..
```

<!-- ```
cmake .. -DTBB_INCLUDE_DIRS=../oneTBB-2019_U9/include -DTBB_LIBRARIES=../oneTBB-2019_U9/build/linux_intel64_gcc_cc12.3.0_libc2.35_kernel6.5.0_release

cmake .. -DTBB_LIBRARIES=../oneTBB-2019_U9/build/linux_intel64_gcc_cc12.3.0_libc2.35_kernel6.5.0_release

cmake .. -DTBB_INCLUDE_DIRS=/home/lionlai/sfm_mvs_slam/mvs-texturing/oneTBB-2019_U9/include -DTBB_LIBRARY=/home/lionlai/sfm_mvs_slam/mvs-texturing/oneTBB-2019_U9/build/linux_intel64_gcc_cc12.3.0_libc2.35_kernel6.5.0_release

make -j
``` -->

### Prepare Camera Parameters

```
python3 -m venv venv_preprocess_data \
&& source venv_preprocess_data/bin/activate \
&& python3 -m pip install numpy

python3 preprocess_camera_parameters.py \
--extrinsic-path ./small_meeting_room_7f/colmap/sparse/0/images.txt \
--intrinsic-path ./small_meeting_room_7f/colmap/sparse/0/distort_cameras.txt \
--output-dir-cam ./img_cam
```

### how to run

```
mkdir -p output/small_meeting_room_7f/wall_0 \
&& ./build/apps/texrecon/texrecon \
--tone_mapping=gamma \
--outlier_removal=gauss_damping \
--view_selection_model \
./input_data/small_meeting_room_7f_img_cam/ \
./input_data/small_meeting_room_7f_ply/wall_0.ply \
./output/small_meeting_room_7f/wall_0/textured_mesh

mkdir -p output/victor/wall_2 \
&& ./build/apps/texrecon/texrecon \
--tone_mapping=gamma \
--outlier_removal=gauss_damping \
--view_selection_model \
./input_data/victor_img_cam/ \
./input_data/victor_ply/wall_2.ply \
./output/victor/wall_2/textured_mesh

mkdir -p output/twhg/walls_floors_objects \
&& ./build/apps/texrecon/texrecon \
--tone_mapping=gamma \
--outlier_removal=gauss_damping \
--view_selection_model \
./input_data/twhg_img_cam/ \
./input_data/twhg_ply/walls_floors_objects.ply \
./output/twhg/walls_floors_objects/textured_mesh

mkdir -p output/twhg/walls_floors \
&& ./build/apps/texrecon/texrecon \
--tone_mapping=gamma \
--outlier_removal=gauss_damping \
--view_selection_model \
--keep_unseen_faces \
./input/twhg_img_cam/ \
./input/twhg_ply/walls_floors.ply \
./output/twhg/walls_floors/textured_mesh

rm -r output/twhg/thick \
&& mkdir -p output/twhg/thick \
&& ./build/apps/texrecon/texrecon \
--tone_mapping=gamma \
--outlier_removal=gauss_damping \
--view_selection_model \
--keep_unseen_faces \
./input/twhg_img_cam/ \
./input/twhg_ply/thick.ply \
./output/twhg/thick/textured_mesh

rm -r output/twhg/thick \
&& mkdir -p output/twhg/thick \
&& ./build/apps/texrecon/texrecon \
--tone_mapping=gamma \
--outlier_removal=gauss_damping \
--view_selection_model \
./input/twhg_img_cam/ \
./input/twhg_ply/thick.ply \
./output/twhg/thick/textured_mesh
```

### 
https://gist.github.com/mbinna/c61dbb39bca0e4fb7d1f73b0d66a4fd1

https://cliutils.gitlab.io/modern-cmake/chapters/intro/dodonot.html