cmake_minimum_required (VERSION 2.6)
project (Parallel-Range-Light)

add_definitions(-g -O2 -Wall -std=c++11 -fcilkplus -DCILKP)

add_subdirectory (external)
add_subdirectory (lib)
add_subdirectory (examples)
