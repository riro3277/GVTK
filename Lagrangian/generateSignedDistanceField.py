# Standalone python script to get signed distance fields (SDFs) for an anatomical vascular model - SDF is used by classGVTKCollision.py

import sys
import pyvista as pv

import numpy as np
from tqdm import tqdm
import igl

if __name__ == '__main__':
    mesh = pv.read("/mnt/d/lagrangiantest/PoiseuilleFlow/vel_long/u_long.vtu")
    V, F = igl.read_triangle_mesh("/mnt/d/lagrangiantest/PoiseuilleFlow/Surface_long.stl")
    points = mesh.points
    sdf, _, _ = igl.point_mesh_squared_distance(points, V, F)
    sdf =  sdf**0.5
    mesh['sdf'] = sdf
    print(sdf.max(), sdf.min())
    mesh.save("/mnt/d/lagrangiantest/PoiseuilleFlow/SDF_long.vtu")
