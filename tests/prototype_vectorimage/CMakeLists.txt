cmake_minimum_required(VERSION 2.8)
 
project(prototype_vectorimage)
 
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
 
add_executable(prototype_vectorimage MACOSX_BUNDLE prototype_vectorimage.cxx)
target_link_libraries(prototype_vectorimage ${ITK_LIBRARIES})
