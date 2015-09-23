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

externalproject_add(ext_coldet
    PREFIX          coldet
    URL             https://downloads.sourceforge.net/project/coldet/coldet/2.0/coldet20.zip
    URL_MD5         37646a7dd046d9c81fca9d55346a108a
    UPDATE_COMMAND  ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/elibs/CMakeLists/coldet.txt CMakeLists.txt
    SOURCE_DIR      ${CMAKE_SOURCE_DIR}/elibs/coldet
    INSTALL_COMMAND ""
)

#find_package eigen
externalproject_add(ext_eigen
    PREFIX          eigen
    URL             https://bitbucket.org/eigen/eigen/get/3.2.5.tar.bz2
    URL_MD5         21a928f6e0f1c7f24b6f63ff823593f5
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
