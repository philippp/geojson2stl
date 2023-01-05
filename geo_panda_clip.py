#!/usr/bin/env python3
import matplotlib.pyplot as plt
import geopandas
from shapely.geometry import Polygon
#outname="presidio"
#y1=37.790722743099
#x1=-122.52012364102214
#y2=37.781673725269
#x2=-122.49339272007701

#outname = "sutro_peaks"
#y1, x1 = (37.766023988337054, -122.46643656789747)
#y2, x2 = (37.74540047934423, -122.43972063315981)

#outname = "telegraph_hill"
#y1, x1 = (37.80575948592944, -122.41176169834966)
#y2, x2 = (37.79858782981106, -122.40116664264737)

#OUTPUT_NAME=sutro_peaks
outname = "mt_sutro_big"
y1, x1 = (37.76403904569938, -122.46628694839559)
y2, x2 = (37.7483558360443, -122.44232174467042)

gdf = geopandas.read_file('geojson/sf.geojson')
polygon = Polygon([(x1, y1), (x1, y2), (x2, y2), (x2, y1), (x1, y1)])
#poly_gdf = geopandas.GeoDataFrame([1], geometry=[polygon], crs=gdf)
gdf_clipped = gdf.clip(polygon)
gdf_clipped.to_file(f"geojson/{outname}.geojson", driver='GeoJSON')
