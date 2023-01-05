#!/usr/bin/env python3
import numpy as np
import open3d as o3d
output_path = "output/presidio/"
input_path = "output/presidio/features.xyz"
point_cloud = np.loadtxt(input_path, skiprows=1)
pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(point_cloud[:,:3])
#o3d.visualization.draw_geometries([pcd])


pcd.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=0.01, max_nn=30))
pcd.orient_normals_to_align_with_direction()
#bpa_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(pcd,o3d.utility.DoubleVector([radius, radius * 2]))

def poisson():
    poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=15, width=0, scale=1.1, linear_fit=False)[0]
    poisson_mesh.compute_triangle_normals()
    bbox = pcd.get_axis_aligned_bounding_box()
    p_mesh_crop = poisson_mesh.crop(bbox)
    p_mesh_crop.paint_uniform_color([1, 0.706, 0])
    p_mesh_crop.compute_triangle_normals()
    o3d.visualization.draw_geometries([p_mesh_crop])
    o3d.io.write_triangle_mesh('output/presidio/poisson.stl', p_mesh_crop)

def ball_pivot():
    distances = pcd.compute_nearest_neighbor_distance()
    avg_dist = np.mean(distances)
    radius = avg_dist/2
    bpa_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(pcd,o3d.utility.DoubleVector([radius, radius * 2]))    
    bbox = pcd.get_axis_aligned_bounding_box()
    p_mesh_crop = bpa_mesh.crop(bbox)
    p_mesh_crop.paint_uniform_color([1, 0.706, 0])
    p_mesh_crop.compute_triangle_normals()
    o3d.visualization.draw_geometries([p_mesh_crop])
    o3d.io.write_triangle_mesh('output/presidio/ball_pivot.stl', p_mesh_crop)

poisson()
