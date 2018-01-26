cmake_minimum_required(VERSION 2.8.3)
project(soundloc)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wno-narrowing")

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  message_generation
  json_msgs
)

###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if you package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
## DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES soundloc
   CATKIN_DEPENDS roscpp rospy std_msgs message_runtime json_msgs
#  DEPENDS system_lib
)

###########
## Build ##
###########
INCLUDE(FindPkgConfig)
pkg_check_modules (GLIB2 glib-2.0)
pkg_check_modules (JACK jack)
pkg_check_modules (OPENCV opencv)
pkg_check_modules (PLPLOT plplotd-c++)
pkg_check_modules (FFTW fftw3)
pkg_check_modules (SNDFILE sndfile)
pkg_check_modules (SAMPLERATE samplerate)
pkg_check_modules (JSON jsoncpp)

## Specify additional locations of header files
## Your package locations should be listed before other locations
# include_directories(include)
include_directories(
  src/mcmcda/
  src/multisoundloc/
  ${GLIB2_INCLUDE_DIRS}
  ${JACK_INCLUDE_DIRS}
  ${OPENCV_INCLUDE_DIRS}
  ${PLPLOT_INCLUDE_DIRS}
  ${JSON_INCLUDE_DIRS}
  ${catkin_INCLUDE_DIRS}
)

## Declare a cpp library
add_library(utils STATIC
   src/mcmcda/utils.cpp
)

add_library(munkres STATIC
   src/mcmcda/munkres.cpp
)

add_library(multisoundloc STATIC
   src/multisoundloc/multisoundloc.cpp
)

## Declare a cpp executable
add_executable(soundloc_node src/soundloc_multi.cpp)
add_dependencies(soundloc_node ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

## Add cmake target dependencies of the executable/library
## as an example, message headers may need to be generated before nodes
# add_dependencies(soundloc_node soundloc_generate_messages_cpp)

## Specify libraries to link a library or executable target against
target_link_libraries(munkres
   ${OPENCV_LIBRARIES}
   ${catkin_LIBRARIES}
 )
target_link_libraries(multisoundloc
   munkres
   utils
   ${GLIB2_LIBRARIES}
   ${JACK_LIBRARIES}
   ${OPENCV_LIBRARIES}
   ${PLPLOT_LIBRARIES}
   ${FFTW_LIBRARIES}
   ${SNDFILE_LIBRARIES}
   ${SAMPLERATE_LIBRARIES}
   ${catkin_LIBRARIES}
 )

target_link_libraries(soundloc_node
   multisoundloc
   ${GLIB2_LIBRARIES}
   ${JSON_LIBRARIES}
   ${catkin_LIBRARIES}
 )

