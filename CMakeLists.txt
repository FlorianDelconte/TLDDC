
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
  pkg_check_modules( EIGEN3 REQUIRED eigen3 )
ELSE(CMAKE_PKGCONFIG_EXECUTABLE)
        FIND_PACKAGE(EIGEN3 REQUIRED)
ENDIF(CMAKE_PKGCONFIG_EXECUTABLE)

include_directories( ${EIGEN3_INCLUDE_DIRS} )


FIND_PACKAGE(GSL REQUIRED)
INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
LINK_DIRECTORIES(${GSL_LIBRARIES})

find_package(PCL 1.3 REQUIRED COMPONENTS common segmentation)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS program_options)
if ( Boost_FOUND )
  message(STATUS "Boost and boost_program_options found.")
  include_directories( ${Boost_INCLUDE_DIRS} )
      SET(DGtalToolsLibDependencies ${DGtalToolsLibDependencies}
     ${Boost_LIBRAIRIES}
     ${Boost_PROGRAM_OPTIONS_LIBRARY})
   SET(DGtalLibInc ${Boost_INCLUDE_DIRS})
endif( Boost_FOUND)

SET(PCLLib ${PCL_SEGMENTATION_LIBRARIES})

find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )

#ADD_EXECUTABLE(segmentation Main Statistic IOHelper DefectSegmentation SegmentationAbstract Centerline/Centerline)
#TARGET_LINK_LIBRARIES(segmentation ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies} ${PCLLib} ${GSL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

#ADD_EXECUTABLE(segcyl MainCylinder Statistic IOHelper DefectSegmentationCylinder SegmentationAbstract Centerline/Centerline)
#TARGET_LINK_LIBRARIES(segcyl ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies} ${PCLLib} ${GSL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

ADD_EXECUTABLE(segunroll SegmentationAbstract IOHelper MainUnroll Statistic  DefectSegmentationUnroll UnrolledMap SegmentationAbstract  Centerline/Centerline)
TARGET_LINK_LIBRARIES(segunroll  ${OpenCV_LIBS} ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies} ${PCLLib} ${GSL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} )

#ADD_EXECUTABLE(colorizeMesh colorizeMesh IOHelper)
#TARGET_LINK_LIBRARIES(colorizeMesh ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies})
