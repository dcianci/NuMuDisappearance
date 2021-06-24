add_library (ubLib OscTools.cc XMLTools.cc MicroBooNE_dis.cc)

include_directories ( ${CMAKE_SOURCE_DIR}/inc)
include_directories ( ${CMAKE_SOURCE_DIR}/tinyxml)

target_link_libraries(ubLib ${ROOT_LIBRARIES})
target_link_libraries(ubLib ${GSL_LIBRARIES})
