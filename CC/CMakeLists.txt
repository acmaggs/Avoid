cmake_minimum_required(VERSION 3.26)
project(CC)

set(CMAKE_CXX_STANDARD 23)

add_executable(CC main.cpp
        lattice.cc
        graphics.cc
        avoid.cc)
