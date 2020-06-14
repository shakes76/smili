PROJECT(milxImageViewer)
 
#find libraries
include(find)

# Add Headers into List
INCLUDE_DIRECTORIES(${MILXQT_INCLUDE_DIRS})
LINK_DIRECTORIES(${MILXQT_INSTALL_BINARY_DIR})
 
ADD_EXECUTABLE(milxImageViewer milxImageViewer.cpp)
TARGET_LINK_LIBRARIES(milxImageViewer milx-Qt ${QT_LIBRARIES} milx-SMILI ${ITK_LIBRARIES} vtk-ext ${VTK_LIBRARIES})
