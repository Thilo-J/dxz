cmake_minimum_required(VERSION 3.4...3.18)
project(dxz)

add_subdirectory(pybind11)

set(Headers
    my_dxz.h
    sparse_matrix.h
    quad_linked_list.h
)

set(Source
    my_dxz.cpp
    sparse_matrix.c
    quad_linked_list.c
)

pybind11_add_module(dxz ${Source})