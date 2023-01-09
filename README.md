# geojson2stl

This script attempt to reconstruct an STL mesh from GeoJSON polygon data. The GeoJSON data is converted to XYZ coordinates and loaded as an Open3D point cloud. The script then attempts a Poisson reconstruction of a triangle mesh using the Open3d library. 

The state of this software is very much "works on this machine," but I didn't see any good resources or tutorials out there so I decided to share this as well as (a writeup)[https://medium.com/@philippp/converting-geojson-to-stl-for-3d-printing-1ad56e4d3058].
