cmake_minimum_required(VERSION 3.1)

project(Texturing)
include(ExternalProject)

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'RELWITHDEBINFO' as none was specified.")
    set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif()


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED yes)

set(RESEARCH "OFF" CACHE BOOL "Use the gco library for Multi-Label Optimization, which is licensed only for research purposes!")
if(RESEARCH)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DRESEARCH=1")
endif()

FIND_PACKAGE(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wundef -pedantic")
endif()

find_package(PNG REQUIRED)
find_package(JPEG REQUIRED)
find_package(TIFF REQUIRED)

add_subdirectory(elibs)
include_directories(SYSTEM
    ${CMAKE_SOURCE_DIR}/elibs/rayint/libs
    ${CMAKE_SOURCE_DIR}/elibs/mve/libs
    ${CMAKE_SOURCE_DIR}/elibs/eigen
)

include_directories(
    libs
)

link_directories(
    ${CMAKE_SOURCE_DIR}/elibs/mve/libs/mve
    ${CMAKE_SOURCE_DIR}/elibs/mve/libs/util
)

if(RESEARCH)
    include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/elibs/gco)
    link_directories(${CMAKE_BINARY_DIR}/elibs/ext_gco/src/ext_gco-build)
endif()

add_subdirectory(libs)
add_subdirectory(apps)

if(RESEARCH)
    message(
"
******************************************************************************
 Due to use of the -DRESEARCH=ON option, the resulting program is licensed
 for research purposes only. Please pay special attention to the gco license.
******************************************************************************
")
endif()
