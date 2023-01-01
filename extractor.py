#!/usr/bin/env python3
import json
import argparse
import logging
import pdb
import datetime
import collections
import geo_xy
import pprint
import geojsonreader

def load_geojson(filename: str) -> dict:
    geojson = json.loads(open(filename,'r').read())
    assert geojson['type'] == 'FeatureCollection', geojson.get('type')
    assert 'metadata' in geojson.keys(), "Run geojsonreader.py first"
    return geojson

def build_island_trees(geojson):
    """Builds a tree of overlapping polygons."""
    # Step 1: build a map of layers to features for easy iteration.
    layermap = collections.defaultdict(list)
    objectidmap = dict()
    
    for f in geojson['features']:
        layermap[f['properties']['layer_idx']].append(f)
        objectidmap[f['properties']['objectid']] = f
    layer_indices = sorted(layermap.keys(), reverse=True)
    island_roots = list()
    # Step 2: Starting at the top layer, iterate through features in the layer and attempt to map them to the lower layer.
    for layer_idx_idx in range(len(layer_indices)-1):
        layer_idx = layer_indices[layer_idx_idx]
        lower_layer_idx = layer_indices[layer_idx_idx+1]
        for f in layermap[layer_idx]:
            found_base = False
            for f_lower in layermap[lower_layer_idx]:
                # If an arbitrary point on f is inside or on f_lower,
                # f is collected to f_lower.
                f_coords = tuple(f['geometry']['coordinates'][0])
                if geo_xy.check_inside(f_lower['geometry']['coordinates'], f_coords):
                    if 'tree_children' not in f_lower['properties'].keys():
                        f_lower['properties']['tree_children'] = list()
                    f_lower['properties']['tree_children'].append(f['properties']['objectid'])
                    logging.info(f"Feature {f['properties']['objectid']} sits on {f_lower['properties']['objectid']}")
                    found_base = True
            if not found_base:
                logging.info(f"Feature {f['properties']['objectid']} appears to be a floating island starting at {layer_idx}.")
                island_roots.append(f)
    island_roots += layermap[layer_indices[-1]]
    island_trees = list()
    for island_root in island_roots:
        island_trees.append(get_tree_child_nodes(island_root, objectidmap))
    return island_trees

def get_tree_child_nodes(node, objectidmap, nodes=None):
    if not nodes:
        nodes = list()
    nodes.append(node)
    logging.info(f"Entering get_tree_child_nodes for {node['properties']['objectid']}.")    
    if 'tree_children' in node['properties'].keys():
        for objectid in node['properties'].get('tree_children',[]):
            cur_child_node = objectidmap[objectid]
            logging.info(f"Recursing into {cur_child_node['properties']['objectid']} from {node['properties']['objectid']}.")
            nodes += get_tree_child_nodes(cur_child_node, objectidmap)
    else:
        logging.info(f"Feature {node['properties']['objectid']} has no child nodes.")
    return nodes

def build_island_geojson_obj(geojson, island_tree: list):
    geojson_new = dict(geojson)
    geojson_new['features'] = island_tree
    geojson_new['metadata']['tree_root_objectid'] = island_tree[0]['properties']['objectid']
    return geojsonreader.update_metadata(geojson_new)

def main():
    arg_parser = argparse.ArgumentParser(
        prog="GeoJSON2STL",
        description="Attempts to break a GeoJSON file into distinct islands.")
    arg_parser.add_argument("-i", "--input_geojson",
                            help="input geojson file.",
                            required=True)
    arg_parser.add_argument("-o", "--output_dir",
                            help="output directory (will be created if needed).")
    arg_parser.add_argument("-l", "--log",
                            help="Set loglevel as DEBUG, INFO, WARNING, ERROR, or CRITICAL.",
                            default="INFO")

    args = arg_parser.parse_args()
    if getattr(args, 'log'):
        loglevel = args.log
        numeric_loglevel = getattr(logging, loglevel.upper(), None)
        if not isinstance(numeric_loglevel, int):
            raise ValueError('Invalid log level: %s' % loglevel)
        logging.basicConfig(level=numeric_loglevel)
    if getattr(args, 'output_dir'):
        output_dir = args.output_dir
    else:
        output_dir = "output-" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S")

    geojson = load_geojson(args.input_geojson)
    island_trees = build_island_trees(geojson)
    for island_tree in island_trees:
        logging.info(f"Exporting island_{island_tree[0]['properties']['objectid']}")
        island_geojson = build_island_geojson_obj(geojson, island_tree)
        island_id = island_geojson['metadata'].get('tree_root_objectid',island_geojson['features'][0]['properties']['objectid'])
        filename = f"{args.output_dir}/island_{island_id}.json"
        open(filename,'w').write(json.dumps(island_geojson, indent=2))

if __name__ == "__main__":
    main()
