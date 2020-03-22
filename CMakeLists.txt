cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
set(CMAKE_BUILD_TYPE Debug) 
set (CMAKE_CXX_STANDARD 11) 

project(background_subtraction)

find_package(Boost COMPONENTS iostreams serialization chrono REQUIRED) 
find_package(PCL 1.8 REQUIRED) 
find_package(OpenCV REQUIRED ) 

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable (background_subtraction background_subtraction.cpp)
target_link_libraries (background_subtraction ${PCL_LIBRARIES} ${OpenCV_LIBS})
