cmake_minimum_required(VERSION 3.14)
project(mpi_test C)

set(CMAKE_C_STANDARD 11)

string(APPEND CMAKE_C_FLAGS "-fopenmp")
add_subdirectory(multi_grid)
add_subdirectory(tSolver)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

add_executable(mpi_test main.c init.c)
target_link_libraries(mpi_test ${MPI_LIBRARIES} m tSolver multi_grid)

if(MPI_COMPILE_FLAGS)
    set_target_properties(mpi_test PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(mpi_test PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

