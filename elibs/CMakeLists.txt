if(RESEARCH)
    externalproject_add(ext_gco
        PREFIX          gco
        URL             http://vision.csd.uwo.ca/code/gco-v3.0.zip
        URL_MD5         10e071892c38f076d30a91ca5351a847
        UPDATE_COMMAND  ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/elibs/CMakeLists/gco.txt CMakeLists.txt
        SOURCE_DIR      ${CMAKE_SOURCE_DIR}/elibs/gco
        INSTALL_COMMAND ""
    )
endif()

externalproject_add(ext_rayint
    PREFIX          rayint
    GIT_REPOSITORY  https://github.com/nmoehrle/rayint.git
    GIT_TAG         cuda
    UPDATE_COMMAND  ""
    SOURCE_DIR      ${CMAKE_SOURCE_DIR}/elibs/rayint
    CONFIGURE_COMMAND ""
    BUILD_COMMAND   ""
    INSTALL_COMMAND ""
)

#find_package eigen
externalproject_add(ext_eigen
    PREFIX          eigen
    URL             https://bitbucket.org/eigen/eigen/get/3.2.8.tar.bz2
    URL_MD5         9e3bfaaab3db18253cfd87ea697b3ab1
    SOURCE_DIR      ${CMAKE_SOURCE_DIR}/elibs/eigen
    CONFIGURE_COMMAND ""
    BUILD_COMMAND   ""
    INSTALL_COMMAND ""
)

externalproject_add(ext_mve
    PREFIX          mve
    GIT_REPOSITORY  https://github.com/nmoehrle/mve.git
    UPDATE_COMMAND  ""
    SOURCE_DIR      ${CMAKE_SOURCE_DIR}/elibs/mve
    CONFIGURE_COMMAND ""
    BUILD_COMMAND   make -C libs/mve && make -C libs/util #not platform independent
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND ""
)
