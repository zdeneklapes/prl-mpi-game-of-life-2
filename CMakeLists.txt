cmake_minimum_required(VERSION 3.6)
project(life)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O0 -g")
#set(CMAKE_PREFIX_PATH "/usr/local/Cellar/open-mpi/5.0.2/lib")
set(CMAKE_PREFIX_PATH "/usr/bin/mpiexec")

#find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

add_executable(${CMAKE_PROJECT_NAME} life.cpp)
#message(STATUS "MPI_LIBRARIES: ${MPI_LIBRARIES}")
#target_link_libraries(${CMAKE_PROJECT_NAME} OpenMP::OpenMP_CXX ${MPI_LIBRARIES})
target_link_libraries(${CMAKE_PROJECT_NAME} ${MPI_LIBRARIES})
