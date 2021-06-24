cmake_minimum_required (VERSION 3.10.2)
project (GlobalFit)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
          set (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
    endif ()
endif ()



FIND_PACKAGE(ROOT)# COMPONENTS Core Cint RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread)
#Get root libs
if(ROOT_FOUND)
	message("Found ROOT using CMAKE FIND_PACKAGE")
	include(${ROOT_USE_FILE})
  set( ROOT_LIBRARIES ${ROOT_LIBRARIES} -lMathMore -lMinuit)
ELSE()
	message("ROOT not compiled using CMAKE, use root-config")
	#Get our ROOT libraries and such
	FIND_PROGRAM(ROOT_CONFIG root-config)
	EXEC_PROGRAM(${ROOT_CONFIG} ARGS --cflags --glibs OUTPUT_VARIABLE ROOT_CXX_FLAGS)
	message("ROOT cxx flags set to: ${ROOT_CXX_FLAGS}")
ENDIF()

#Get our GSL libraries
FIND_PROGRAM(GSL_CONFIG gsl-config)
EXEC_PROGRAM(${GSL_CONFIG} ARGS --cflags OUTPUT_VARIABLE GSL_CXX_FLAGS)
EXEC_PROGRAM(${GSL_CONFIG} ARGS --libs OUTPUT_VARIABLE GSL_LIBS)


# if something is not working, it's probably because you're not running linux and will have to change the line below to reflect the location of libGSL
set (CMAKE_CXX_FLAGS "-m64 -fsigned-char -fPIC -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm -std=c++11")

message("CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

include_directories(${PROJECT_BINARY_DIR})
include_directories(${CMAKE_SOURCE_DIR}/inc/)

add_subdirectory(src)
add_subdirectory(tinyxml)
add_subdirectory(fitters)
