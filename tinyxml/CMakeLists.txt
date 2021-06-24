# The extension is already found. Any number of sources could be listed here.
add_library (TinyXMLLib tinyxml2.cpp)

# Make sure the compiler can find include files for our Hello library
# when other libraries or executables link to Hello
#target_include_directories (tinyxmllib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
include_directories ( ${CMAKE_CURRENT_SOURCE_DIR})
