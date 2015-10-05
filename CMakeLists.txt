cmake_minimum_required(VERSION 2.8)

project(Texturing)
include(ExternalProject)

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'RELWITHDEBINFO' as none was specified.")
    set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif()

set(RESEARCH "OFF" CACHE BOOL "Use the gco library for Multi-Label Optimization, which is licensed only for research purposes!")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11_FLAG)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X_FLAG)
if(COMPILER_SUPPORTS_CXX11_FLAG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    set(COMPILER_SUPPORTS_CXX11 TRUE)
elseif(COMPILER_SUPPORTS_CXX0X_FLAG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    set(COMPILER_SUPPORTS_CXX11 TRUE)
elseif(MSVC AND NOT (MSVC_VERSION LESS 1800))
    set(COMPILER_SUPPORTS_CXX11 TRUE)
endif()

if(NOT COMPILER_SUPPORTS_CXX11)
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

if(RESEARCH)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DRESEARCH=1")
endif()

FIND_PACKAGE(OpenMP)
if(OPENMP_FOUND)
    foreach(CBT DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
        set(CMAKE_C_FLAGS_${CBT} "${CMAKE_C_FLAGS_${CBT}} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS_${CBT} "${CMAKE_CXX_FLAGS_${CBT}} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS_${CBT} "${CMAKE_EXE_LINKER_FLAGS_${CBT}} ${OpenMP_EXE_LINKER_FLAGS}")
    endforeach(CBT)
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wundef -pedantic")
    set(CMAKE_C_FLAGS_PROFILE "-pg -O3")
    set(CMAKE_CXX_FLAGS_PROFILE "-pg -O3")
    set(CMAKE_EXE_LINKER_FLAGS_PROFILE "-pg")
    set(CMAKE_SHARED_LINKER_FLAGS_PROFILE "-pg")
endif()

find_package(PNG REQUIRED)
find_package(JPEG REQUIRED)
find_package(TIFF REQUIRED)

add_subdirectory(elibs)
include_directories(SYSTEM
    ${CMAKE_SOURCE_DIR}/elibs/coldet/src
    ${CMAKE_SOURCE_DIR}/elibs/mve/libs
    ${CMAKE_SOURCE_DIR}/elibs/eigen
)

include_directories(
    libs
)

link_directories(
    ${CMAKE_BINARY_DIR}/elibs/coldet/src/ext_coldet-build/
    ${CMAKE_SOURCE_DIR}/elibs/mve/libs/mve
    ${CMAKE_SOURCE_DIR}/elibs/mve/libs/util
)

if(RESEARCH)
    include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/elibs/gco)
    link_directories(${CMAKE_BINARY_DIR}/elibs/gco/src/ext_gco-build)
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
