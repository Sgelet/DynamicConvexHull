cmake_minimum_required(VERSION 3.24)
project(DynamicConvexHull)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3")

find_package(CGAL)

add_executable(DynamicConvexHull include/AvlTree.h include/CHTree.h main.cpp include/util.h include/CQTree.h
        include/inexact.h)

target_link_libraries(DynamicConvexHull CGAL::CGAL)