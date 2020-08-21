externalproject_add(ext_mapmap
    PREFIX          ext_mapmap
    GIT_REPOSITORY  https://github.com/PolyCam/mapmap_cpu.git
    GIT_TAG         master
    UPDATE_COMMAND  ""
    SOURCE_DIR      ${CMAKE_SOURCE_DIR}/elibs/mapmap
    CONFIGURE_COMMAND ""
    BUILD_COMMAND   ""
    INSTALL_COMMAND ""
)

externalproject_add(ext_rayint
    PREFIX          ext_rayint
    GIT_REPOSITORY  https://github.com/PolyCam/rayint.git
    GIT_TAG         cuda
    UPDATE_COMMAND  ""
    SOURCE_DIR      ${CMAKE_SOURCE_DIR}/elibs/rayint
    CONFIGURE_COMMAND ""
    BUILD_COMMAND   ""
    INSTALL_COMMAND ""
)

externalproject_add(ext_eigen
    PREFIX          ext_eigen
    GIT_REPOSITORY  https://github.com/PolyCam/eigen
    GIT_TAG         3.3.6
    SOURCE_DIR      ${CMAKE_SOURCE_DIR}/elibs/eigen
    CONFIGURE_COMMAND ""
    BUILD_COMMAND   ""
    INSTALL_COMMAND ""
)

externalproject_add(ext_mve
    PREFIX          mve
    GIT_REPOSITORY  https://github.com/PolyCam/mve.git
    UPDATE_COMMAND  ""
    SOURCE_DIR      ${CMAKE_SOURCE_DIR}/elibs/mve
    CONFIGURE_COMMAND ""
    BUILD_COMMAND   make -C libs/mve && make -C libs/util #not platform independent
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND ""
)
