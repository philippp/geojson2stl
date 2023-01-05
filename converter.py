#!/usr/bin/env python3
import pdb
import json
import math
import pprint
import numpy
import argparse
import os
import datetime
import logging
import geo_xy
import open3d as o3d
import collections
import numpy as np
import math

TOTAL_HEIGHT = 1
scaled_side = 1000
fn_feature_coord = lambda f: f['geometry']['coordinates'][0]
fn_feature_elevation = lambda f: int(f['properties']['elevation'])

class CoordinateConverter:
    def __init__(self, setting="latlon"):
        self._converter = getattr(self, f"convert_{setting}", False)
        self.setting = setting
        assert self._converter

    def convert(self, geojson, lonlat):
        return self._converter(geojson, lonlat)

    def convert_latlon(self, geojson, lonlat):
        return (lonlat[0], lonlat[1])

    def convert_tanlon(self, geojson, lonlat):
        x = (math.tan(lonlat[0]) - math.tan(geojson['metadata']['min_lon']))
        y = (lonlat[1] - geojson['metadata']['min_lat'])
        return (x,y)

    def convert_straight(self, geojson, lonlat):
        x = (lonlat[0] - geojson['metadata']['min_lon'])
        y = (lonlat[1] - geojson['metadata']['min_lat'])
        return (x,y)

# XY[Z] points are specified as (x,y[,z])
# Lines are specified as ((x,y[,z]),(x,y[,z]))

def load_geojson(filename: str, coordinate_converter: CoordinateConverter) -> dict:
    geojson = json.loads(open(filename,'r').read())
    assert geojson['type'] == 'FeatureCollection', geojson.get('type')
    assert 'metadata' in geojson.keys(), "Run geojsonreader.py first"
    return geojson

def render_feature_facets(feature: dict, objectid_to_feature) -> str:
    lines = ""

    elevation_z = feature['properties']['elevation_z']
    height_z = feature['properties']['height_z']
    bottom_z = elevation_z
    top_z = elevation_z + height_z

    for tree_child_id in feature['properties'].get('tree_children',[]):
        child_z = fn_feature_coord(objectid_to_feature[tree_child_id])[0][2]
        top_z = max(top_z, child_z)
    for idx in range(1,len(fn_feature_coord(feature))):
        lines += write_rect_facets(
            fn_feature_coord(feature)[idx-1],
            fn_feature_coord(feature)[idx], top_z, bottom_z)
    lines += write_rect_facets(
        fn_feature_coord(feature)[-1],
        fn_feature_coord(feature)[0], top_z, bottom_z)
    return lines

def write_rect_facets(xy1, xy2, height_top, height_bottom=0):
    f1_p3 = numpy.array([xy2[0],xy2[1],height_top])
    f1_p2 = numpy.array([xy1[0],xy1[1],height_bottom])
    f1_p1 = numpy.array([xy2[0],xy2[1],height_bottom])
    f1_surface = numpy.cross(f1_p2-f1_p1, f1_p3-f1_p1)
    f1_surface_normal = f1_surface / f1_surface.sum()

    f2_p3 = numpy.array([xy2[0],xy2[1],height_top])
    f2_p2 = numpy.array([xy1[0],xy1[1],height_top])
    f2_p1 = numpy.array([xy1[0],xy1[1],height_bottom])
    f2_surface = numpy.cross(f2_p2-f2_p1, f2_p3-f2_p1)
    f2_surface_normal = f2_surface / f2_surface.sum()
    return f"""
facet normal {f1_surface_normal[0]} {f1_surface_normal[1]} {f1_surface_normal[2]}
    outer loop
        vertex {xy2[0]} {xy2[1]} {height_bottom}
        vertex {xy1[0]} {xy1[1]} {height_bottom}
        vertex {xy2[0]} {xy2[1]} {height_top}
    endloop
endfacet
facet normal {f2_surface_normal[0]} {f2_surface_normal[1]} {f2_surface_normal[2]}
    outer loop
        vertex {xy1[0]} {xy1[1]} {height_bottom}
        vertex {xy1[0]} {xy1[1]} {height_top}
        vertex {xy2[0]} {xy2[1]} {height_top}
    endloop
endfacet
"""

def write_stl(feature_list: list, filename: str):
    lines = ""
    objectid_to_feature = dict()
    for f in feature_list:
        objectid_to_feature[f['properties']['objectid']] = f
    for feature in feature_list:
        lines += render_feature_facets(feature, objectid_to_feature)
    with open(filename,'w') as f:
        f.write("solid myshape\n")
        f.write(lines)
        f.write("endsolid myshape\n\n")

def render_feature_svg(feature, min_x, min_y, scaling_factor):
    lines = ""
    poly_coord_string_list = list()
    for cur_coords in fn_feature_coord(feature):
        poly_coord_string_list.append("%f,%f" % (
            (cur_coords[0] - min_x) * scaling_factor,
            (cur_coords[1] - min_y) * scaling_factor))
    poly_coord_string = " ".join(poly_coord_string_list)
    return f'<polyline points="{poly_coord_string}" style="stroke:black;fill:none;stroke-width:1;" />\n'

def write_svg(feature_list: list, filename: str, scale_to=600):
    any_coord = fn_feature_coord(feature_list[0])[0]
    min_x = any_coord[0]
    min_y = any_coord[1]
    max_x = any_coord[0]
    max_y = any_coord[1]
    for feature in feature_list:
        for coord in fn_feature_coord(feature):
            if coord[0] < min_x:
                min_x = coord[0]
            if coord[0] > max_x:
                max_x = coord[0]
            if coord[1] < min_y:
                min_y = coord[1]
            if coord[1] > max_y:
                max_y = coord[1]
    if (max_x - min_x) > (max_y - min_y):
        scaling_factor = scale_to/(max_x - min_x)
    else:
        scaling_factor = scale_to/(max_y - min_y)

    with open(filename,'w') as f:
        f.write(f'<svg viewBox="0 0 {scale_to} {scale_to}" xmlns="http://www.w3.org/2000/svg">\n')
        for feature in feature_list:
            f.write(render_feature_svg(feature, min_x, min_y, scaling_factor))
        f.write('</svg>\n')

def write_json(feature_list: list, filename: str):
    with open(filename,'w') as f:
        f.write(json.dumps(dict(feature_list=feature_list), indent=2))

def scale_features(geojson, converter, scale_to):
    # Scale X and Y. We scale to a fixed width, and allow height to be
    # proportionate.
    geo_meta = geojson['metadata']
    x1, y1 = converter.convert(geojson, (geo_meta['max_lon'], geo_meta['max_lat']))
    x2, y2 = converter.convert(geojson, (geo_meta['min_lon'], geo_meta['min_lat']))
    min_x, max_x = sorted([abs(x1),abs(x2)])
    min_y, max_y = sorted([abs(y1),abs(y2)])
    longest_edge_m = 0
    if (max_y - min_y) > (max_x - min_x):
        # Y is bigger (map is taller)
        unit_size = scale_to / (max_y - min_y)
        longest_edge_m = geo_xy.haversine([[geo_meta['max_lon'],geo_meta['max_lat']],[geo_meta['max_lon'],geo_meta['min_lat']]]) * 1000
        logging.info(f"Y > X, unit_size={unit_size}, distance={longest_edge_m}")
    else:  # X is bigger (map is wider)
        unit_size = scale_to / (max_x - min_x)
        longest_edge_m = geo_xy.haversine([[geo_meta['max_lon'],geo_meta['max_lat']],[geo_meta['min_lon'],geo_meta['max_lat']]]) * 1000
        logging.info(f"X > Y, unit_size={unit_size}, distance={longest_edge_m}")
    if not longest_edge_m:
        longest_edge_m = 100
    geo_meta['scale_unit_size'] = unit_size
    new_z_max = ((geo_meta['max_elevation_m'] - geo_meta['min_elevation_m'])/longest_edge_m)*scale_to
    for f in geojson['features']:
        f_elevation_m = f['properties']['elevation_m']
        z_value = ((f_elevation_m - geo_meta['min_elevation_m']) / longest_edge_m) * scale_to
        f['properties']['elevation_z'] = z_value
        f['properties']['height_z'] = (f['properties']['height_m'] / longest_edge_m) * scale_to
        coordinates = fn_feature_coord(f)
        for idx in range(len(coordinates)):
            new_x, new_y = converter.convert(geojson, coordinates[idx])
            coordinates[idx][0] = new_x * unit_size
            coordinates[idx][1] = new_y * unit_size
            coordinates[idx].append(z_value)
    return geojson

def augment_feature_points(geojson):
    # Get the current average point density and XYZ limits.
    xyz_list = list()
    max_x, max_y, max_z = (0,0,0)
    features = filter(lambda f: f['geometry']['type'].lower() == 'polygon', geojson['features'])
    for coord_list in [fn_feature_coord(f) for f in features]:
        for c in coord_list:
            xyz_list.append([c[0], c[1], c[2]])
            if c[0] > max_x:
                max_x = c[0]
            if c[1] > max_y:
                max_y = c[1]
            if c[2] > max_z:
                max_z = c[2]
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(xyz_list)
    distances = pcd.compute_nearest_neighbor_distance()
    avg_distance = np.mean(distances)
    assert avg_distance > 0, f"Distance is {distance}"
    logging.info(f"max X,Y,Z=({max_x},{max_y},{max_z}) avg_dist={avg_distance}")
    
    # We create a rectangle of points at the highest z-index, then
    # iterate through them and drop them.
    new_points = list()
    cur_x, cur_y = (0,0)
    while cur_x < max_x:
        cur_x += avg_distance
        while cur_y < max_y:
            cur_y += avg_distance
            z_index = z_at_point(features, cur_x, cur_y)
            if z_index >= 0:
                new_points.append((cur_x, cur_y, z_index))
    logging.info(f"Added {len(new_points)} out of a possible {math.floor(max_x/avg_distance)*math.floor(max_y/avg_distance)} points.")
    return new_points

def z_at_point(features, cur_x, cur_y):
    # Sort the polygon layers by z-index, highest to lowest.
    layer_map = collections.defaultdict(list)
    for f in features:
        layer_map[f['properties']['layer_idx']].append(f)
    layer_indices = sorted(layer_map.keys(), reverse=True)
    for layer_index in layer_indices:
        for cur_feature in layer_map[layer_index]:
            if geo_xy.check_inside(fn_feature_coord(cur_feature), (cur_x, cur_y)):
                return fn_feature_coord(cur_feature)[0][2]  # Any Z-index will do.
    return -1

def export_xyz(geojson, filename):
    new_xyz = augment_feature_points(geojson)
    with open(filename,'w') as f:
        f.write("X Y Z\n")
        for feature in geojson['features']:
            for coords in fn_feature_coord(feature):
                f.write(f"{coords[0]} {coords[1]} {coords[2]}\n")
        for coords in new_xyz:
            f.write(f"{coords[0]} {coords[1]} {coords[2]}\n")

def main():
    arg_parser = argparse.ArgumentParser(
        prog="GeoJSON2XYZ",
        description="Attempts to convert a GeoJSON file into a shape file.")
    arg_parser.add_argument("-i", "--input_geojson",
                            help="input geojson file.",
                            required=True)
    arg_parser.add_argument("-f", "--feature_list",
                            help="extract only the given features (CSV by index).")
    arg_parser.add_argument("-o", "--output_dir",
                            help="output directory (will be created if needed).")
    arg_parser.add_argument("-d", "--debug",
                            help="write out each feature as an SVG file.",
                            action="store_true")
    arg_parser.add_argument("-l", "--log",
                            help="Set loglevel as DEBUG, INFO, WARNING, ERROR, or CRITICAL.")
    arg_parser.add_argument("-c", "--coordinate_converter",
                            help="Coordinate converter (latlon, straight, tanlon)",
                            default="latlon")
    arg_parser.add_argument("-s", "--scale",
                            help="Scale the largest dimension to this size",
                            default=1, type=int)

    args = arg_parser.parse_args()
    if getattr(args, 'log'):
        loglevel = args.log
        numeric_loglevel = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_loglevel, int):
            raise ValueError('Invalid log level: %s' % loglevel)
        logging.basicConfig(level=numeric_loglevel)

    coordinate_converter = CoordinateConverter(args.coordinate_converter)
    geojson = load_geojson(args.input_geojson, coordinate_converter)
    scale_features(geojson, coordinate_converter, getattr(args, "scale", 500))

    feature_indices = list()
    if getattr(args, 'feature_list'):
        feature_indices = [int(f.strip()) for f in args.feature_list.split(",")]

    if getattr(args, 'output_dir'):
        output_dir = args.output_dir
    else:
        output_dir = "output-" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    features_to_export = list()
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for idx in range(len(geojson['features'])):
        if feature_indices and not idx in feature_indices:
            continue
        else:
            feature = geojson['features'][idx]
            if getattr(args, 'debug'):
                write_svg([feature], f"{output_dir}/feature_{idx}.svg")
                write_json([feature], f"{output_dir}/feature_{idx}.json")
            else:
                features_to_export.append(feature)
    if features_to_export:
        export_xyz(geojson, f"{output_dir}/features.xyz")        
        write_stl(features_to_export, f"{output_dir}/export.stl")
        write_svg(features_to_export, f"{output_dir}/features.svg")
        write_json(features_to_export, f"{output_dir}/features.json")

if __name__ == "__main__":
    main()
