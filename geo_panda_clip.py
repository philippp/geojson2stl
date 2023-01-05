#!/usr/bin/env python3
import matplotlib.pyplot as plt
import geopandas
from shapely.geometry import Polygon

y1=37.790722743099
x1=-122.52012364102214
y2=37.781673725269
x2=-122.49339272007701

gdf = geopandas.read_file('geojson/sf.geojson')
polygon = Polygon([(x1, y1), (x1, y2), (x2, y2), (x2, y1), (x1, y1)])
#poly_gdf = geopandas.GeoDataFrame([1], geometry=[polygon], crs=gdf)
gdf_clipped = gdf.clip(polygon)
gdf_clipped.to_file("geojson/presidio.geojson", driver='GeoJSON')
