# Plane Detection from RGBD

### Compile
Use the CMake file from the parent directory for compilation.

### Run test
```
./test_plane_depth [color_file] [depth_file]
```

### Usage
Include Plane.h, and call ExtractPlane, with input and output specified as follow:

Input
1. depth_input, depth image
2. fx fy cx cy, intrinsic parameters

Output, a pair of
1. a mask image with per-pixel label of Plane ID. 255 means not a planar region.
2. a vector of PlaneHelper, in which 3D plane parameters (a,b,c,d) of ax+by+cz+d=0 is included in params.

