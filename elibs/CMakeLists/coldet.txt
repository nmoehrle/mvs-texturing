cmake_minimum_required(VERSION 2.8)

project(ColDet)

set(SOURCES
    src/box_bld.cpp
    src/box.cpp
    src/cdmath3d.cpp
    src/coldet_bld.cpp
    src/coldet.cpp
    src/multiobject.cpp
    src/mytritri.cpp
    src/sysdep.cpp
    src/tritri.c
)

if(UNIX)
    add_definitions(-DGCC)
endif()

set(CMAKE_BUILD_TYPE Release)
add_library(coldet STATIC ${SOURCES})
