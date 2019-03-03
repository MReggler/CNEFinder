#!/bin/bash

echo "Building cnefinder image..."
docker build --rm=true --file docker/Dockerfile -t cnefinder .
echo "Built cnefinder image."


# Add more docker build commands below
