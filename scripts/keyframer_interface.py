import json
import os
import argparse
import imageio
import numpy as np
import shutil

""" 
Reads in the camera.json files written by the ondevice keyframer and writes them to the SCENE_FOLDER format expected by texrecon. 

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

python3 keyframer_interface.py ../data/office/keyframes/ ../data/office/keyframes/cameras_texrecon

"""


def write_new_cam(new_cam_path, j, img_shape):
    print(j)
    cam_file = open(new_cam_path, 'w+')
    first_line = '{} {} {} {} {} {} {} {} {} {} {} {}'.format(
        j['t_03'], j['t_13'], j['t_23'], j['t_00'], j['t_01'], j['t_02'], j['t_10'], j['t_11'], j['t_12'], j['t_20'], j['t_21'], j['t_22'])
    cam_file.write(first_line)
    cam_file.write('\n')
    d0 = 0
    d1 = 0
    # texrecon wants center point normalized by image size
    ppx = j['cx'] / img_shape[1]
    ppy = j['cy'] / img_shape[0]
    f = j['fx'] / max(img_shape)
    paspect = 1
    second_line = '{} {} {} {} {} {}'.format(f, d0, d1, paspect, ppx, ppy)
    cam_file.write(second_line)


def main(args):
    os.makedirs(args.output_folder, exist_ok=True)
    input_cam_folder = os.path.join(args.keyframe_folder, 'cameras')
    input_img_folder = os.path.join(args.keyframe_folder, 'images')
    input_cams = os.listdir(input_cam_folder)
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
        shutil.copyfile(img_path, new_img_path)
        img_shape = imageio.imread(img_path).shape
        print(img_shape)
        write_new_cam(new_cam_path, data, img_shape)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('keyframe_folder')
    parser.add_argument('output_folder')
    args = parser.parse_args()
    main(args)
