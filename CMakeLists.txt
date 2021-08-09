cmake_minimum_required(VERSION 3.10)

# set the project name
project(tigre VERSION 1.0)


file(GLOB bfc
    "./src/basic_format_converter.cpp"
    "./src/globals.cpp"
    "./src/math/math.cpp"
    "./src/properties.cpp"
    "./src/atom.cpp"
    "./src/frame.cpp"
    "./src/trajectory.cpp"
    "./src/rw.cpp"
)


set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)


# add the executable
add_executable(bfc  "${bfc}")
target_include_directories(bfc PUBLIC "./src" )

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set (CMAKE_CXX_FLAGS_RELEASE "-Ofast   -pthread")
set (CMAKE_CXX_FLAGS "-Ofast    -pthread")


