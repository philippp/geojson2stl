#!/usr/bin/env python3

import pdb
import json
import logging
import geo_xy
import argparse
import math

# Top left: 37.74153079183624, -122.45901985734845
# Botom Right: 37.736389093383096, -122.45053356855081
import inspect


fn_feature_coord = lambda f: f['geometry']['coordinates']
ELEVATION_SCALE_TO_M = 0.3048  # Feet to meters - change depending on fn_feature_elevation source.
fn_feature_elevation = lambda f: int(f['properties']['elevation']) * ELEVATION_SCALE_TO_M

def load(filename: str):
    geojson = json.loads(open(filename,'r').read())
    assert geojson['type'] == 'FeatureCollection', geojson.get('type')
    return update_metadata(geojson)

def update_metadata(geojson):
    """
    Update boundary data in metadata.
    """
    min_lon, min_lat = fn_feature_coord(geojson['features'][0])[0]
    max_lon, max_lat = (min_lon, min_lat)
    max_elevation = -999
    min_elevation = 999
    all_elevations = set()
    geojson['metadata'] = dict()
    for f in geojson['features']:
        for c_lon, c_lat in fn_feature_coord(f):
            if c_lon < min_lon:
                min_lon = c_lon
            if c_lon > max_lon:
                max_lon = c_lon
            if c_lat < min_lat:
                min_lat = c_lat
            if c_lat > max_lat:
                max_lat = c_lat
        cur_feature_elevation = fn_feature_elevation(f)
        f['properties']['elevation_m'] = cur_feature_elevation
        if cur_feature_elevation > max_elevation:
            max_elevation = cur_feature_elevation
        if cur_feature_elevation < min_elevation:
            min_elevation = cur_feature_elevation
        all_elevations.add(cur_feature_elevation)
    geojson['metadata']['max_elevation_m'] = max_elevation
    geojson['metadata']['min_elevation_m'] = min_elevation
    geojson['metadata']['max_lon'] = max_lon
    geojson['metadata']['min_lon'] = min_lon
    geojson['metadata']['max_lat'] = max_lat
    geojson['metadata']['min_lat'] = min_lat
    logging.info(f"Map Bounds: {min_lon} <= lon <= {max_lon}, {min_lat} <= lat <= {max_lat}")

    elevation_list = sorted(all_elevations)
    logging.info(f"Elevations: {elevation_list}")
    step_sizes = set(int(elevation_list[i] - elevation_list[i-1]) for i in range(1, len(elevation_list)))
    if len(step_sizes) > 1:
        logging.info(f"Found {len(step_sizes)} distinct elevation step changes: {list(step_sizes)}")
    logging.info(f"{len(elevation_list)} layers have data.")

    # Longitude gets larger going from west to east (left to right)
    # Latitude gets larger from south to north (bottom to top)
    # We will start our coordinate system with (0,0) at the bottom left.

    for idx in range(len(geojson['features'])):
        f = geojson['features'][idx]
        new_coordinates = list()
        feature_coords = fn_feature_coord(f)
        for c_lon, c_lat in feature_coords:
            new_coordinates.append((c_lon,c_lat))
        if feature_coords[0] == feature_coords[-1]:
            f['geometry']['type'] = "Polygon"
        f['properties']['layer_idx'] = elevation_list.index(
            fn_feature_elevation(f))
    return geojson

def prune(geojson, left_top, right_bottom):
    feature_list = list()
    for feature in geojson['features']:
        feature['geometry']['coordinates'] = _prune_coords(feature['geometry']['coordinates'], left_top, right_bottom)
        if feature['geometry']['coordinates']:
            logging.debug(f"[prune] Keeping feature idx={geojson['features'].index(feature)}")
            feature_list.append(feature)
    logging.info(f"[{inspect.stack()[0][3]}] Kept {len(feature_list)} out of {len(geojson['features'])} features.")
    geojson['features'] = feature_list
    return geojson

def _prune_coords(coords, left_top, right_bottom):
    # Find coords outside of the box
    outside_coords = list()
    for coord_pair in coords:
        if not geo_xy.is_point_in_box(coord_pair, left_top, right_bottom):
            outside_coords.append(coord_pair)
    if len(outside_coords) == len(coords):
        return list() # the entire polygon is outside of the bounding box.
    coords = geo_xy.cut_polygon(coords, left_top[0], 0, 'ge')
    coords = geo_xy.cut_polygon(coords, left_top[1], 1, 'le')
    coords = geo_xy.cut_polygon(coords, right_bottom[0], 0, 'le')
    coords = geo_xy.cut_polygon(coords, right_bottom[1], 1, 'ge')
    return coords

def project_latlon_to_xy(geojson, latlon):
    x = (latlon[0] - geojson['metadata']['min_lat']) * geojson['metadata']['scale_unit_size']
    y = (math.tan(latlon[1]) - math.tan(geojson['metadata']['min_lon'])) * geojson['metadata']['scale_unit_size']
    logging.debug(f"[{inspect.stack()[1][3]}] Translated {latlon} to {x},{y}")
    return (x,y)

def main():
    arg_parser = argparse.ArgumentParser(
        prog="GeoJSONClipper",
        description="Attempts to clip a GeoJSON file to a rectangle, modifying shapes if needed.")
    arg_parser.add_argument("-i", "--input_geojson",
                            help="input geojson file.",
                            required=True)
    arg_parser.add_argument("-b", "--bounding_box",
                            help="Specify bounding box as (lat,long),(lat,long) northwest to southeast")
    arg_parser.add_argument("-o", "--output_geojson",
                            help="output geojson file.")
    arg_parser.add_argument("-l", "--log",
                            help="Set loglevel as DEBUG, INFO, WARNING, ERROR, or CRITICAL.",
                            default="INFO")
    arg_parser.add_argument("-g", "--geometry_type",
                            help="Only export the specified geometry type.")
    args = arg_parser.parse_args()
    if getattr(args, 'log'):
        loglevel = args.log
        numeric_loglevel = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_loglevel, int):
            raise ValueError('Invalid log level: %s' % loglevel)
        logging.basicConfig(level=numeric_loglevel)
    logging.info(f"[main] Loading {args.input_geojson}...")
    geojson = load(args.input_geojson)
    feature_indices = list()
    if getattr(args, 'bounding_box'):
        boundary_start_str, boundary_end_str = [
            s.strip(')').strip('(') for s in args.bounding_box.split("),(")]
        boundary_start_latlon = [float(l) for l in boundary_start_str.split(',')]
        boundary_end_latlon = [float(l) for l in boundary_end_str.split(',')]
        left_top = (boundary_start_latlon[1], boundary_start_latlon[0])
        right_bottom = (boundary_end_latlon[1], boundary_end_latlon[0])
        prune(geojson, left_top, right_bottom)
    if args.geometry_type:
        geojson['features'] = list(
            filter(lambda f: f['geometry']['type'] == args.geometry_type,
                   geojson['features']))
        logging.info(f"[main] {len(geojson['features'])} features have geometry type {args.geometry_type}")
    # After clipping to the boundary, rebuild the metadata.
    update_metadata(geojson)
    open(args.output_geojson,'w').write(json.dumps(geojson, indent=2))

if __name__ == "__main__":
    main()
