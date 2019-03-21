# CENG 789 Assignment \#1

This is the CENG 789 Assignment 1 implemented by Alp Eren SARI.
The project depends on [libigl](http://libigl.github.io/libigl/)  

## Impelemented Stuff

In this project, geodesic distance calculation using dijkstra shortest path
search by using array, minimum heap and fibonacci heap implemented. For the
query point 3 different method is used and calculation duration's printed
on the screen. After that part, FPS runs. Using FPS points later ISO-Curve
Signature and Bilateral Histogram descriptors calculated and printed in a
text file

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `example_bin` binary.

## Run

From within the `build` directory just issue:

    ./example_bin mesh_path source_point destination_point fps_point_number

    where mesh_path is the path of the input mesh source_point is the source
    vertex number, destination_point is the destination vertex number and
    fps_point_number is the fps point numbers to be sampled.

A glfw app should launch displaying a 3D cube.

## Dependencies

The only dependencies are stl, eigen, [libigl](http://libigl.github.io/libigl/) and
the dependencies of the `igl::opengl::glfw::Viewer`.

We recommend you to install libigl using git via:

    git clone https://github.com/libigl/libigl.git
    cd libigl/
    git submodule update --init --recursive
    cd ..

If you have installed libigl at `/path/to/libigl/` then a good place to clone
this library is `/path/to/libigl-example-project/`.
