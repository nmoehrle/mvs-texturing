cmake_minimum_required(VERSION 2.8)

project(GraphCutOptimization)

set(SOURCES
    GCoptimization.cpp
    graph.cpp
    LinkedBlockList.cpp
    maxflow.cpp
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGCO_ENERGYTYPE=float")
set(CMAKE_BUILD_TYPE Release)
add_library(gco STATIC ${SOURCES})
