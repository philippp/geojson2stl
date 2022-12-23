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

TOTAL_HEIGHT = 10
fn_feature_coord = lambda f: f['geometry']['coordinates']
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
        y = (geojson['metadata']['max_lat'] - lonlat[1])
        return (x,y)

    def convert_straight(self, geojson, lonlat):
        x = (lonlat[0] - geojson['metadata']['min_lon'])
        y = (geojson['metadata']['max_lat'] - lonlat[1])
        return (x,y)



def load(filename: str, converter: CoordinateConverter):
    geojson = json.loads(open(filename,'r').read())
    assert geojson['type'] == 'FeatureCollection', geojson.get('type')
    return preprocess(geojson, converter)    

def preprocess(geojson, converter: CoordinateConverter):
    """
    Convert geojson to x,y coordinates and regularize elevation to Z.
    Note: GeoJSON stores points in (lon,lat)!
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
        if cur_feature_elevation > max_elevation:
            max_elevation = cur_feature_elevation
        if cur_feature_elevation < min_elevation:
            min_elevation = cur_feature_elevation
        all_elevations.add(cur_feature_elevation)
    geojson['metadata']['max_elevation'] = max_elevation
    geojson['metadata']['min_elevation'] = min_elevation
    geojson['metadata']['max_lon'] = max_lon
    geojson['metadata']['min_lon'] = min_lon
    geojson['metadata']['max_lat'] = max_lat
    geojson['metadata']['min_lat'] = min_lat
    logging.info(f"Map Bounds: {min_lon} <= lon <= {max_lon}, {min_lat} <= lat <= {max_lat}")

    # Compute the z-scaling factor, assuming evenly spaced layers.
    elevation_list = sorted(all_elevations)
    logging.info(f"Elevations: {elevation_list}")
    elevation_step_size = elevation_list[1] - elevation_list[0]
    for i in range(1, len(elevation_list)):
        if elevation_list[i] != elevation_list[i-1] + elevation_step_size:
            logging.warning(f"At i={i} we stepped more than {elevation_step_size}")
            logging.warning(f"{elevation_list[i]} != {elevation_list[i-1]} + {elevation_step_size}")
    geojson['metadata']['layer_height'] = TOTAL_HEIGHT / len(elevation_list)
    logging.info(f"{len(elevation_list)} layers, each has height of {geojson['metadata']['layer_height']}")

    # Longitude gets larger going from west to east (left to right)
    # Latitude gets larger from south to north (bottom to top)
    # We will start our coordinate system with (0,0) at the bottom left.
    
    for idx in range(len(geojson['features'])):
        f = geojson['features'][idx]
        new_coordinates = list()
        feature_coords = fn_feature_coord(f)
        for c_lon, c_lat in feature_coords:
            x, y = converter.convert(geojson, (c_lon, c_lat))
            new_coordinates.append(x,y)
        if feature_coords[0] == feature_coords[-1]:
            f['geometry']['type'] = "Polygon"
        f['geometry']['orig_coordinates'] = list(f['geometry']['coordinates'])
        f['geometry']['coordinates'] = new_coordinates
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
    pdb.set_trace()
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

def scale_features(geojson, converter, scale_to):
    for f in geojson['features']:
        f['geometry']['coordinates'] = list(f['geometry']['orig_coordinates'])
    preprocess(geojson, converter)
    # Scale X and Y. We scale to a fixed width, and allow height to be
    # proportionate.
    x1, y1 = converter.convert(geojson, (geojson['metadata']['max_lon'], geojson['metadata']['max_lat']))
    x2, y2 = converter.convert(geojson, (geojson['metadata']['min_lon'], geojson['metadata']['min_lat']))
    min_x, max_x = sorted([abs(x1),abs(x2)])
    min_y, max_y = sorted([abs(y1),abs(y2)])
    
    if (max_y - min_y) > (max_x - min_x):
        unit_size = scale_to / (max_y - min_y)
        logging.info(f"vertical distance (latitude or Y) is bigger, unit_size={unit_size}")
    else:
        unit_size = scale_to / (max_x - min_y)
        logging.info(f"horizontal distance (longitude or X) is bigger, unit_size={unit_size}")        
    geojson['metadata']['scale_unit_size'] = unit_size
    for f in geojson['features']:
        for idx in range(len(f['geometry']['coordinates'])):
            f['geometry']['coordinates'][idx][0] = f['geometry']['coordinates'][idx][0] * unit_size
            f['geometry']['coordinates'][idx][1] = f['geometry']['coordinates'][idx][1] * unit_size            
    return geojson

    
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
    arg_parser.add_argument("-s", "--scale",
                            help="Scale the largest dimension to this size",
                            default=0, type=int)
    arg_parser.add_argument("-g", "--geometry_type",
                            help="Only export the specified geometry type.")
    arg_parser.add_argument("-c", "--coordinate_converter",
                            help="Coordinate converter (latlon, straight, tanlon)",
                            default="latlon")
    
    args = arg_parser.parse_args()
    if getattr(args, 'log'):
        loglevel = args.log
        numeric_loglevel = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_loglevel, int):
            raise ValueError('Invalid log level: %s' % loglevel)
        logging.basicConfig(level=numeric_loglevel)
    logging.info(f"[main] Loading {args.input_geojson}...")
    coordinate_converter = CoordinateConverter(args.coordinate_converter)
    geojson = load(args.input_geojson, coordinate_converter, scale_to=args.scale)

    logging.info(f"[main] Converted lonlat to {args.coordinate_converter}...")    
    feature_indices = list()
    if getattr(args, 'bounding_box'):
        boundary_start_str, boundary_end_str = [
            s.strip(')').strip('(') for s in args.bounding_box.split("),(")]
        boundary_start_latlon = [float(l) for l in boundary_start_str.split(',')]
        boundary_end_latlon = [float(l) for l in boundary_end_str.split(',')]
        left_top = coordinate_converter.convert(
            geojson, (boundary_start_latlon[1], boundary_start_latlon[0]),
            scale_to=geojson['metadata']['scale_unit_size'])
        right_bottom = coordinate_converter.convert(
            geojson, (boundary_end_latlon[1], boundary_end_latlon[0]),
            scale_to=geojson['metadata']['scale_unit_size'])
        prune(geojson, left_top, right_bottom)
    if args.geometry_type:
        geojson['features'] = list(
            filter(lambda f: f['geometry']['type'] == args.geometry_type,
                   geojson['features']))
        logging.info(f"[main] {len(geojson['features']} match geometry type {args.geometry_type}")
    scale_features(geojson, coordinate_converter, 600)
    open(args.output_geojson,'w').write(json.dumps(geojson, indent=2))

if __name__ == "__main__":
    main()

    