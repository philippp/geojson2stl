#!/usr/bin/env bash

#COORDINATES="(37.804930911836074,-122.4117456980303),(37.79952149168777, -122.4014226393132)"
#OUTPUT_NAME=telegraph_hill

#COORDINATES="(37.79076129566718, -122.51905187441855),(37.780866070250774,-122.49282844703384)"
#OUTPUT_NAME=lands_end

#COORDINATES="(37.76442037411654, -122.46389229150489),(37.75423124691226, -122.45128635520909)"
#OUTPUT_NAME=mt_sutro

#COORDINATES="(37.766023988337054, -122.46643656789747),(37.74540047934423, -122.43972063315981)"
#OUTPUT_NAME=sutro_peaks

#COORDINATES="(37.83348957189002,-122.3896992971493),(37.80584159740146,-122.35380299136703)"
#OUTPUT_NAME=yerba_buena

COORDINATES="(37.74530140127512, -122.4192342655966),(37.741685578532106, -122.40845878054377)"
OUTPUT_NAME=bernal

mkdir -p geojson/${OUTPUT_NAME}
./geojsonreader.py -i geojson/sf.geojson -b"${COORDINATES}" -o geojson/${OUTPUT_NAME}/polygons.geojson -g Polygon
./extractor.py -i geojson/${OUTPUT_NAME}/polygons.geojson -o geojson/${OUTPUT_NAME} --log info
mkdir -p output/${OUTPUT_NAME}
for i in `ls geojson/${OUTPUT_NAME}/island*.json`; do
    island_id=$(basename $i .json)
    echo "Processing ${island_id}"
    mkdir -p output/${OUTPUT_NAME}/${island_id}
    ./converter.py --input_geojson $i --log INFO -c straight -o output/${OUTPUT_NAME}/${island_id} -s 1
done
./converter.py --input_geojson geojson/${OUTPUT_NAME}/polygons.geojson --log INFO -c straight -o output/${OUTPUT_NAME}
echo "All done!"

