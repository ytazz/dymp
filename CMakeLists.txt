cmake_minimum_required(VERSION 2.8.12) # Default version on Ubuntu 14.04

project(dymp)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

set(DYMP_BUILD_SAMPLE false CACHE BOOL "build sample programs")

if(WIN32)
  if(NOT CMAKE_DEBUG_POSTFIX)
    set(CMAKE_DEBUG_POSTFIX d)
  endif()
endif()

if(WIN32)
    set(CMAKE_CXX_FLAGS "/openmp /wd4530")
endif()
if(UNIX)
	set(CMAKE_CXX_FLAGS "-std=c++11 -Wall")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-invalid-offsetof -Wno-sign-compare -Wno-delete-non-virtual-dtor -Wno-reorder")
	set(CMAKE_CXX_FLAGS_DEBUG "-g")
	set(CMAKE_CXX_FLAGS_RELEASE "-O2")
endif()

find_package(Boost  REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(MKL    REQUIRED)

if(DYMP_BUILD_SAMPLE)
	find_package(OpenGL  REQUIRED)
	find_package(glfw3   REQUIRED)
endif()

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${Boost_INCLUDE_DIR})
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
  include/blas.h
  include/canvas.h
  include/centroid.h
  include/constraint.h
  include/dymp.h
  include/id.h
  include/link.h
  include/model.h
  include/rollpitchyaw.h
  include/solver.h
  include/timer.h
  include/types.h
  include/util.h
  include/variable.h
  include/wholebody.h
  include/world.h
  )

if(DYMP_BUILD_SAMPLE)
	set(DYMP_SOURCES
		${DYMP_SOURCES}
		src/app.cpp
		)
	set(DYMP_HEADERS
		${DYMP_HEADERS}
		include/app.h
		)	
endif()

add_library(dymp
  ${DYMP_HEADERS}
  ${DYMP_SOURCES}
 )
 
target_link_libraries(dymp
	Eigen3::Eigen
    ${MKL_LIBRARIES}
    ${Boost_LIBRARIES}
)
if(DYMP_BUILD_SAMPLE)
	target_link_libraries(dymp
		OpenGL::GL
	    ${GLFW_LIBRARIES}
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

