cmake_minimum_required(VERSION 3.21)
project(GAMES102-hw-7-10 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# partially compile the easy3d components
set(EASY3D_BUILD_VIEWER ON  CACHE BOOL "" FORCE)
set(EASY3D_BUILD_ALGO   ON  CACHE BOOL "" FORCE)
set(EASY3D_BUILD_TESTS  OFF CACHE BOOL "" FORCE)
set(EASY3D_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

add_subdirectory(${CMAKE_SOURCE_DIR}/../third_party/easy3d 
                 ${CMAKE_BINARY_DIR}/easy3d)

add_library(hw_common INTERFACE)
target_link_libraries(hw_common INTERFACE easy3d::viewer easy3d::algo)

add_subdirectory(hw8)
