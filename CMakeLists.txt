cmake_minimum_required(VERSION 2.8.12) # Default version on Ubuntu 14.04

project(dymp)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

set(CMAKE_CXX_STANDARD 20)

set(DYMP_BUILD_SAMPLE false CACHE BOOL "build sample programs")

if(WIN32)
  if(NOT CMAKE_DEBUG_POSTFIX)
    set(CMAKE_DEBUG_POSTFIX d)
  endif()
endif()

if(WIN32)
  set(CMAKE_CXX_FLAGS "/openmp /wd4530 /wd4819")
  option(ENABLE_MSVC_AVX  "Enable AVX instructions on VC++." ON)
  option(ENABLE_MSVC_AVX2 "Enable AVX2 instructions on VC++." ON)
  if(ENABLE_MSVC_AVX2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2")
  elseif(ENABLE_MSVC_AVX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX")
  endif()
endif()
if(UNIX)
	set(CMAKE_CXX_FLAGS "-std=c++11 -Wall")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-invalid-offsetof -Wno-sign-compare -Wno-delete-non-virtual-dtor -Wno-reorder")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -fPIC")
	set(CMAKE_CXX_FLAGS_DEBUG "-g")
	set(CMAKE_CXX_FLAGS_RELEASE "-O2")
endif()

#find_package(Boost  REQUIRED)
find_package(Eigen3 REQUIRED)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
find_package(MKL    REQUIRED)

if(DYMP_BUILD_SAMPLE)
	find_package(OpenGL  REQUIRED)
	find_package(glfw3   REQUIRED)
	include_directories(${GLFW_INCLUDE_DIR})
endif()

include_directories(${PROJECT_SOURCE_DIR}/include)
#include_directories(${Boost_INCLUDE_DIR})
include_directories(${MKL_INCLUDE_DIR})
#include_directories(${EIGEN3_INCLUDE_DIR})

set(DYMP_SOURCES
  src/blas.cpp
  src/canvas.cpp
  src/centroid.cpp
  src/constraint.cpp
  src/id.cpp
  src/link.cpp
  src/model.cpp
  src/rollpitchyaw.cpp
  src/solver.cpp
  src/solverddp.cpp
  src/timer.cpp
  src/util.cpp
  src/variable.cpp
  src/wholebody.cpp
  src/world.cpp
  )
set(DYMP_HEADERS
  include/dymp/blas.h
  include/dymp/canvas.h
  include/dymp/centroid.h
  include/dymp/constraint.h
  include/dymp/dymp.h
  include/dymp/id.h
  include/dymp/link.h
  include/dymp/model.h
  include/dymp/rollpitchyaw.h
  include/dymp/solver.h
  include/dymp/timer.h
  include/dymp/types.h
  include/dymp/util.h
  include/dymp/variable.h
  include/dymp/wholebody.h
  include/dymp/world.h
  )

if(DYMP_BUILD_SAMPLE)
	set(DYMP_SOURCES
		${DYMP_SOURCES}
		src/app.cpp
		)
	set(DYMP_HEADERS
		${DYMP_HEADERS}
		include/dymp/app.h
		)	
endif()

add_library(dymp
  ${DYMP_HEADERS}
  ${DYMP_SOURCES}
 )
 
target_link_libraries(dymp
	PRIVATE Eigen3::Eigen
    PRIVATE ${MKL_LIBRARIES}
)
if(DYMP_BUILD_SAMPLE)
	target_link_libraries(dymp
		PRIVATE OpenGL::GL
		PRIVATE glfw      
	    #${GLFW_LIBRARIES}
	)
endif()

install(TARGETS dymp
  EXPORT        dympExport
  RUNTIME       DESTINATION bin
  LIBRARY       DESTINATION lib
  INCLUDES      DESTINATION include/dymp
  PUBLIC_HEADER DESTINATION include/dymp
  )
install(FILES ${DYMP_HEADERS} DESTINATION include/dymp)

install(EXPORT dympExport
  FILE         dympConfig.cmake
  DESTINATION  cmake
  )

if(DYMP_BUILD_SAMPLE)
	add_subdirectory(sample)
endif()

## print all variables
get_cmake_property(_variableNames VARIABLES)
foreach(_variableName ${_variableNames})
	message(STATUS "${_variableName}=${${_variableName}}")
endforeach()

