file (GLOB HEADERS "*.h")
file (GLOB SOURCES "[^_]*.cpp")

set(LIBRARY tex)
add_library(${LIBRARY} STATIC ${SOURCES})
set_property(TARGET ${LIBRARY} PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
add_dependencies(${LIBRARY} ext_mve ext_rayint ext_eigen ext_mapmap)
target_link_libraries(${LIBRARY} -lmve -lmve_util ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${TIFF_LIBRARIES} TBB::tbb)
install(TARGETS ${LIBRARY} ARCHIVE DESTINATION lib)
