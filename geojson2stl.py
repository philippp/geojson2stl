#!/usr/bin/env python3
import pdb
import argparse
import open3d as o3d
import numpy as np
import geopandas
import shapely
import logging
import trimesh

CONV_FT_TO_M = 0.3048

def compute_gdf_pointcloud(gdf: geopandas.GeoDataFrame) -> o3d.geometry.PointCloud:
    """Compute the PointCloud of the GeoDataFrame.
    """
    min_x, min_y = (99999999,99999999)
    min_elevation = min([e for e in gdf['elevation']])
    
    for f in gdf['geometry']:
        for c in f.coords:
            if c[0] < min_x:
                min_x = c[0]
            if c[1] < min_y:
                min_y = c[1]

    logging.info(f"min_x={min_x} min_y={min_y} min_elevation={min_elevation}")
    gdf['flat_coordinates'] = gdf['geometry'].combine(
        gdf['elevation'],
        (lambda g, e: [(float(c[0] - min_x), float(c[1] - min_y), float(e) * CONV_FT_TO_M) for c in g.coords]))
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector()
    for coordinate_list in gdf['flat_coordinates']:
        pcd.points.extend(o3d.utility.Vector3dVector(coordinate_list))
    pcd.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=0.01, max_nn=30))
    pcd.orient_normals_to_align_with_direction()
    return pcd

def compute_poisson_mesh(pcd: o3d.geometry.PointCloud) -> o3d.geometry.TriangleMesh:
    logging.info("Computing Poisson Mesh...")
    poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=12, width=0, scale=1, linear_fit=True)[0]
    bbox = pcd.get_axis_aligned_bounding_box()
    p_mesh_crop = poisson_mesh.crop(bbox)
    p_mesh_crop.paint_uniform_color([1, 0.706, 0])
    p_mesh_crop.compute_triangle_normals()
    o3d.visualization.draw_geometries([p_mesh_crop])
    return p_mesh_crop

def patch_mesh():
    in_trimesh = trimesh.load_mesh()
    in_trimesh.fill_holes()
    return in_trimesh.as_open3d()

def main():
    arg_parser = argparse.ArgumentParser(
        prog="GeoJSON2STL",
        description="Attempts to convert a GeoJSON file into a shape file.")
    arg_parser.add_argument("-i", "--input_geojson",
                            help="input geojson file.",
                            required=True)
    arg_parser.add_argument("-o", "--output_stl",
                            help="output stl file.")
    arg_parser.add_argument("-l", "--log",
                            help="Set loglevel as DEBUG, INFO, WARNING, ERROR, or CRITICAL.")

    args = arg_parser.parse_args()
    if getattr(args, 'log'):
        loglevel = args.log
        numeric_loglevel = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_loglevel, int):
            raise ValueError('Invalid log level: %s' % loglevel)
        logging.basicConfig(level=numeric_loglevel)
    
    gdf = geopandas.read_file(args.input_geojson)
    gdf = gdf.explode(index_parts=True)
    gdf = gdf.to_crs(epsg=3395)  # Mercurial projection
    logging.info(f"{gdf.crs}")
    pcd = compute_gdf_pointcloud(gdf)
    mesh = compute_poisson_mesh(pcd)
    o3d.io.write_triangle_mesh(args.output_stl, mesh)
    logging.info(f"Wrote {args.output_stl}")
if __name__ == "__main__":
    main()
