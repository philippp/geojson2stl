#!/usr/bin/env python3
import re
import pdb
import os
import argparse
import open3d as o3d
import numpy as np
import geopandas
import shapely
import logging
import matplotlib.pyplot as plt
import time

CONV_FT_TO_M = 0.3048  # SFData provides elevation in feet :(

def compute_gdf_pointcloud(gdf: geopandas.GeoDataFrame, z_scale: float) -> o3d.geometry.PointCloud:
    """Compute the Open3D PointCloud of the Geopandas GeoDataFrame.
    """
    min_x, min_y = (99999999,99999999)
    min_elevation = min([e for e in gdf['elevation']])
    
    for f in gdf['geometry']:
        for c in f.coords:
            if c[0] < min_x:
                min_x = c[0]
            if c[1] < min_y:
                min_y = c[1]

    logging.debug(f"min_x={min_x} min_y={min_y} min_elevation={min_elevation}")
    gdf['flat_coordinates'] = gdf['geometry'].combine(
        gdf['elevation'],
        (lambda g, e: [(float(c[0] - min_x),
                        float(c[1] - min_y),
                        float(e) * CONV_FT_TO_M * z_scale) for c in g.coords]))
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector()
    for coordinate_list in gdf['flat_coordinates']:
        pcd.points.extend(o3d.utility.Vector3dVector(coordinate_list))
    pcd.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=0.01, max_nn=30))
    pcd.orient_normals_to_align_with_direction()
    return pcd

def compute_poisson_mesh(pcd: o3d.geometry.PointCloud, depth: int) -> o3d.geometry.TriangleMesh:
    """Compute the mesh of the point cloud.
    depth: The depth of the octree used for the surface reconstruction determines resolution of the resulting triangle mesh.
    """
    poisson_mesh, densities = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(
        pcd, depth=depth, width=0, scale=1.0)
    densities = np.asarray(densities)
    graph_rate = (densities - densities.min()) / (densities.max() - densities.min())
    density_colors = plt.get_cmap('plasma')(graph_rate)
    density_colors = density_colors[:, :3]
    poisson_mesh.vertex_colors = o3d.utility.Vector3dVector(density_colors)
    bbox = pcd.get_axis_aligned_bounding_box()
    p_mesh_crop = poisson_mesh.crop(bbox)
    p_mesh_crop.compute_triangle_normals()
    return p_mesh_crop

def main():
    arg_parser = argparse.ArgumentParser(
        prog=os.path.basename(__file__),
        description="Attempts to convert a GeoJSON file into a shape file.")
    arg_parser.add_argument("-i", "--input_geojson",
                            help="input geojson file.",
                            required=True)
    arg_parser.add_argument("-c", "--clip_coordinates",
                            help="Clip input GeoJSON to these coordinates.",
                            required=False)
    arg_parser.add_argument("-o", "--output_stl",
                            help="Write output stl file.",
                            required=False)
    arg_parser.add_argument("-d", "--detail_level",
                            default=12, type=int,
                            help="Desired level of detail (Octree depth for Poisson surface reconstruction).", required=False)
    arg_parser.add_argument("-z", "--z_scale",
                            default=1, type=float,
                            help="Scaling factor for z-scale",
                            required=False)
    arg_parser.add_argument("-l", "--log",
                            choices=['DEBUG', 'INFO', 'WARNING',
                                     'ERROR', 'CRITICAL'],
                            help="Set loglevel.",
                            default="WARNING")
    arg_parser.add_argument("-g", "--output_clipped_geojson",
                            help="If specified, save the clipped GeoJSON",
                            required=False)

    args = arg_parser.parse_args()
    if getattr(args, 'log'):
        loglevel = args.log
        numeric_loglevel = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_loglevel, int):
            raise ValueError('Invalid log level: %s' % loglevel)
        logging.basicConfig(level=numeric_loglevel)
    logging.info(f"Loading {args.input_geojson}")
    gdf = geopandas.read_file(args.input_geojson)
    if getattr(args, 'clip_coordinates'):
        cc = args.clip_coordinates.strip()
        if (not cc.startswith("(") or not cc.endswith(")") or
            not re.search(r"\)\s*,\s*\(", "), (")):
            logging.critical("Expected (lat_w, long_n),(lat_e, long_s)")
            return
        coords = [float(c.replace(")","").replace("(","").strip()) \
                  for c in cc.split(",")]
        coords = [coords[1], coords[0], coords[3], coords[2]]
        logging.info(f"Clipping to {coords}")
        polygon = shapely.Polygon([(coords[0], coords[1]),
                                   (coords[0], coords[3]),
                                   (coords[2], coords[3]),
                                   (coords[2], coords[1]),
                                   (coords[0], coords[1])])
        gdf = gdf.clip(polygon)
        if getattr(args, 'output_clipped_geojson'):
            gdf.to_file(args.output_clipped_geojson, driver="GeoJSON")
    gdf = gdf.explode(index_parts=True)
    gdf = gdf.to_crs(epsg=3395)  # Mercurial projection
    logging.info(f"{gdf.crs}")
    pcd = compute_gdf_pointcloud(gdf, args.z_scale)
    logging.info(f"Attempting Poisson reconstruction on {len(pcd.points)} points with an octree depth of {args.detail_level}.")
    mesh = compute_poisson_mesh(pcd, args.detail_level)
    target_triangle_count = max(100000, int(len(mesh.triangles)/4))
    logging.info(f"Decimating mesh from {len(mesh.triangles)} to {target_triangle_count}")
    mesh = mesh.simplify_quadric_decimation(target_number_of_triangles=target_triangle_count)
    if args.output_stl:
        o3d.io.write_triangle_mesh(args.output_stl, mesh)
        logging.info(f"Wrote {args.output_stl}")
    o3d.visualization.draw_geometries([mesh])
            
if __name__ == "__main__":
    main()
