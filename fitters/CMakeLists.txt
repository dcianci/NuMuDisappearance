add_executable (PL_v1 plotLooper_v1.cc)
add_executable (uBFitter uBfitter.cc)
add_executable (uBFitter_frequensens uBfitter_frequensens.cc)

include_directories (${CMAKE_SOURCE_DIR}/inc)
include_directories (${CMAKE_SOURCE_DIR}/tinyxml)

target_link_libraries (PL_v1 ubLib)
target_link_libraries (PL_v1 TinyXMLLib)

target_link_libraries (uBFitter ubLib)
target_link_libraries (uBFitter TinyXMLLib)

target_link_libraries (uBFitter_frequensens ubLib)
target_link_libraries (uBFitter_frequensens TinyXMLLib)
