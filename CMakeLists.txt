cmake_minimum_required(VERSION 3.26)
project(equations)

set(CMAKE_CXX_STANDARD 17)

add_executable(equations main.cpp
        Equations.h
        Equations.cpp)
