cmake_minimum_required(VERSION 2.8.3)
project(ekf_tests)

add_definitions(-std=c++11)

find_package(catkin REQUIRED COMPONENTS
  roscpp
)

find_package(cmake_modules REQUIRED)
find_package(Eigen REQUIRED)

catkin_package(
 INCLUDE_DIRS include
 LIBRARIES ekf_tests
#  CATKIN_DEPENDS roscpp
  DEPENDS Eigen
)

include_directories(
include
  ${catkin_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIRS}
)


## EKF Test library
add_library(${PROJECT_NAME}
   include/ekf_tests/tests.cpp
)
 target_link_libraries(${PROJECT_NAME}
   ${catkin_LIBRARIES}
 )
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-O3")

## Sample code
add_executable(ekf_tests_sample src/test_impl.cpp)
target_link_libraries(ekf_tests_sample ${catkin_LIBRARIES} ekf_tests)
