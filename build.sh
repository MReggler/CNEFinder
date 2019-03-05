#!/bin/bash

echo "Building cnefinder image..."
docker build --rm=true --file docker/cnefinder/Dockerfile -t cnefinder .
echo "Built cnefinder image."


echo "Building parse_bed image..."
docker build --rm=true --file docker/parse_bed/Dockerfile -t parse_bed .
echo "Built parse_bed image."

# Add more docker build commands below
