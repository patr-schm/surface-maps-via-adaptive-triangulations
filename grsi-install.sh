#!/bin/bash

# Install dependencies in a first step.
sudo apt-get install cmake g++
sudo apt install libgl1-mesa-dev mesa-utils libglfw3 libglfw3-dev libxinerama-dev libxcursor-dev libxi-dev

# Build project.
mkdir build
cd build
cmake ..
make -j4
