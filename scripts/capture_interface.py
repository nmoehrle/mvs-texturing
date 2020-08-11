import json
import os
import argparse
import imageio
import numpy as np
import open3d as o3d
import shutil
import time
import cv2

"""
Reads in the /keyframes directory written by the ondevice keyframer and writes them to the SCENE_FOLDER format expected by texrecon.

SCENE_FOLDER:
Within a scene folder a .cam file has to be given for each image.
A .cam file is structured as follows:
    tx ty tz R00 R01 R02 R10 R11 R12 R20 R21 R22
    f d0 d1 paspect ppx ppy
First line: Extrinsics - translation vector and rotation matrix (the transform from world to camera)
Second line: Intrinsics - focal length, distortion coefficients, pixel aspect ratio and principal point
The focal length is the distance between camera center and image plane normalized by dividing with the larger image dimension.
For non zero distortion coefficients the image will be undistorted prior to the texturing process. If only d0 is non zero the Noah Snavely's distortion model is assumed otherwise the distortion model of VSFM is assumed.
The pixel aspect ratio is usually 1 or close to 1. If your SfM system doesn't output it, but outputs a different focal length in x and y direction, you have to encode this here.
The principal point has to be given in unit dimensions (e.g. 0.5 0.5).

python3 keyframer_interface.py ../data/office/keyframes/ ../data/office/texrecon

"""


def remesh(in_mesh: str, out_mesh: str, points_to_verts: float = 0.5):
    """
    Remeshes an input mesh by sampling it and applying Poisson Surface meshing
    Args:
        in_mesh: Path to input mesh
        out_mesh: Path to write new mesh to
        points_to_verts: The ratio of points in sampled point cloud to verts in original mesh
    """

    print("Remeshing mesh at {}".format(in_mesh))

    start_load = time.time()
    mesh = o3d.io.read_triangle_mesh(in_mesh)
    print("Time to read mesh: {}".format(time.time() - start_load))
    start_normals = time.time()
    mesh.compute_vertex_normals()
    print("Time to compute mesh normals: {}".format(time.time() - start_normals))
    start_sample = time.time()
    num_points = points_to_verts*len(mesh.vertices)
    pcd = mesh.sample_points_uniformly(number_of_points=30000)
    print("Time to sample mesh: {} with {} points".format(
        time.time() - start_sample, num_points))
    start_mesh = time.time()
    with o3d.utility.VerbosityContextManager(o3d.utility.VerbosityLevel.Debug) as cm:
        mesh, densities = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(
            pcd, depth=12)
    print("Time to run Poisson surface reconstruction: {}".format(
        time.time() - start_mesh))
    start_write = time.time()
    o3d.io.write_triangle_mesh(
        out_mesh, mesh)
    print("Time to write mesh: {}".format(time.time() - start_write))


def extract_transform(j):
    transform = np.zeros((4, 4))
    R = np.zeros((3, 3))

    R[0, 0] = j['t_00']
    R[0, 1] = j['t_01']
    R[0, 2] = j['t_02']
    R[1, 0] = j['t_10']
    R[1, 1] = j['t_11']
    R[1, 2] = j['t_12']
    R[2, 0] = j['t_20']
    R[2, 1] = j['t_21']
    R[2, 2] = j['t_22']
    transform[0, 3] = j['t_03']
    transform[1, 3] = j['t_13']
    transform[2, 3] = j['t_23']
    transform[3, 3] = 1
    transform[0:3, 0:3] = R

    # take inverse to get world to cam
    transform = np.linalg.inv(transform)

    # A rotation by 180 degrees around x axis is necessary
    # to remove an onorthodox convention used by ARKit
    Rx = np.zeros((4, 4))
    Rx[0, 0] = 1.0
    Rx[1, 1] = -1.0
    Rx[2, 2] = -1.0
    Rx[3, 3] = 1.0

    transform = np.matmul(Rx, transform)

    return transform


def write_new_cam(new_cam_path, j, img_shape):
    cam_file = open(new_cam_path, 'w+')
    t = extract_transform(j)
    print('transform \n', t)

    first_line = '{} {} {} {} {} {} {} {} {} {} {} {}'.format(
        t[0, 3], t[1, 3], t[2, 3], t[0, 0], t[0, 1], t[0, 2], t[1, 0], t[1, 1], t[1, 2], t[2, 0], t[2, 1], t[2, 2])
    cam_file.write(first_line)
    cam_file.write('\n')
    d0 = 0
    d1 = 0
    # texrecon wants center point normalized by image size
    ppx = j['cx'] / img_shape[1]
    ppy = j['cy'] / img_shape[0]
    f = j['fx'] / max(img_shape)
    paspect = 1.0
    second_line = '{} {} {} {} {} {}'.format(f, d0, d1, paspect, ppx, ppy)
    cam_file.write(second_line)


def main(args):
    os.makedirs(args.output_folder, exist_ok=True)
    in_mesh = os.path.join(args.capture_folder, 'mesh.obj')
    out_mesh = os.path.join(args.capture_folder, 'mesh.ply')
    remesh(in_mesh, out_mesh)

    input_cam_folder = os.path.join(
        args.capture_folder, 'keyframes', 'cameras')
    input_img_folder = os.path.join(
        args.capture_folder, 'keyframes', 'images')
    input_cams = sorted(os.listdir(input_cam_folder))
    for cam_file in input_cams:
        cam_path = os.path.join(input_cam_folder, cam_file)
        data = {}
        with open(cam_path, 'r') as f:
            data = json.load(f)
        new_cam_path = os.path.join(
            args.output_folder, cam_file.replace('.json', '.cam'))
        img_path = os.path.join(
            input_img_folder, cam_file.replace('.json', '.jpg'))
        new_img_path = os.path.join(
            args.output_folder, cam_file.replace('.json', '.jpg'))
        # load image, resize and write to new img path
        img = cv2.imread(img_path, cv2.IMREAD_COLOR)
        shape = img.shape
        print("Original image shape: ", shape)
        scale_factor = args.max_image_dimension / max(shape)
        new_img = cv2.resize(
            img, (int(scale_factor*shape[1]), int(scale_factor*shape[0])))
        cv2.imwrite(new_img_path, new_img)

        new_img_shape = new_img.shape
        print("new image shape: ", new_img_shape)
        write_new_cam(new_cam_path, data, new_img_shape)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('capture_folder')
    parser.add_argument('output_folder')
    parser.add_argument('--max_image_dimension', type=int, default=1920,
                        help='Resize image so it has max_image_dimension')
    args = parser.parse_args()
    main(args)
