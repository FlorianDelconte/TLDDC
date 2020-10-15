
PROJECT(TreeLogSurfaceDefectsSegmentation)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_CXX_FLAGS "-std=c++11 -pthread -O2")

FIND_PACKAGE(Threads)

FIND_PACKAGE(DGtal REQUIRED)

INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})
FIND_PROGRAM(CMAKE_PKGCONFIG_EXECUTABLE pkg-config)

if(CMAKE_PKGCONFIG_EXECUTABLE)
  find_package( PkgConfig )
  pkg_check_modules( EIGEN REQUIRED eigen3 )
ELSE(CMAKE_PKGCONFIG_EXECUTABLE)
        FIND_PACKAGE(EIGEN REQUIRED)
ENDIF(CMAKE_PKGCONFIG_EXECUTABLE)

include_directories( ${EIGEN_INCLUDE_DIRS} )


FIND_PACKAGE(GSL REQUIRED)
INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
LINK_DIRECTORIES(${GSL_LIBRARIES})

find_package(PCL 1.3 REQUIRED COMPONENTS common segmentation)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

FIND_PACKAGE(Boost 1.65.1 COMPONENTS program_options REQUIRED )
if ( Boost_FOUND )
  message(STATUS "Boost and boost_program_options found.")
  include_directories( ${Boost_INCLUDE_DIRS} )
  SET(DGtalToolsLibDependencies ${DGtalToolsLibDependencies} ${Boost_LIBRAIRIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})
  SET(DGtalLibInc ${Boost_INCLUDE_DIRS})
endif( Boost_FOUND)


SET(PCLLib ${PCL_SEGMENTATION_LIBRARIES})


find_package( OpenCV REQUIRED )#core imgproc highgui imgcodecs
if ( OpenCV_FOUND )
  message(STATUS "OpenCV and boost_program_options found.")
endif( OpenCV_FOUND)

include_directories( ${OpenCV_INCLUDE_DIRS} )
#link_directories(${OpenCV_LIBRARY_DIRS})

#find_package(Python2 COMPONENTS Development )
#include_directories(segunroll PRIVATE ${Python2_INCLUDE_DIRS})

#ADD_EXECUTABLE(segmentation Main Statistic IOHelper DefectSegmentation SegmentationAbstract Centerline/Centerline)
#TARGET_LINK_LIBRARIES(segmentation ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies} ${PCLLib} ${GSL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

#ADD_EXECUTABLE(segcyl MainCylinder Statistic IOHelper DefectSegmentationCylinder SegmentationAbstract Centerline/Centerline)
#TARGET_LINK_LIBRARIES(segcyl ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies} ${PCLLib} ${GSL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

ADD_EXECUTABLE(segunroll SegmentationAbstract IOHelper MainUnroll Statistic  DefectSegmentationUnroll UnrolledMap SegmentationAbstract Centerline/Centerline Centerline/MyPoint2DEmbedderIn3D)#ImageAnalyser
TARGET_LINK_LIBRARIES(segunroll ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies} ${PCLLib} ${GSL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} )

ADD_EXECUTABLE(segToMesh segToMesh IOHelper)
TARGET_LINK_LIBRARIES(segToMesh ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies})

#ADD_EXECUTABLE(offToObj off2obj)
#TARGET_LINK_LIBRARIES(offToObj ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies})

#ADD_EXECUTABLE(colorizeMesh colorizeMesh IOHelper)
#TARGET_LINK_LIBRARIES(colorizeMesh ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies})
