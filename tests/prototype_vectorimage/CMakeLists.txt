project(prototype_vectorimage)

include(find)

add_executable(prototype_vectorimage MACOSX_BUNDLE prototype_vectorimage.cxx)
target_link_libraries(prototype_vectorimage ${ITK_LIBRARIES})
