cmake_minimum_required(VERSION 3.12)
project(StructuredMoyalSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED fftw3)

# Include directories
include_directories(include)
include_directories(${FFTW_INCLUDE_DIRS})

# Source files
set(SOURCES
    src/phase_space.cpp
    src/wigner_df.cpp
    src/propagators.cpp
    src/moyal_solver.cpp
)

# Headers
set(HEADERS
    include/moyal_setup.hpp
    include/phase_space.hpp
    include/wigner_df.hpp
    include/propagators.hpp
    include/moyal_solver.hpp
    include/potentials.hpp
)

# Create library
add_library(moyal_solver ${SOURCES} ${HEADERS})
target_link_libraries(moyal_solver ${FFTW_LIBRARIES})
target_compile_options(moyal_solver PRIVATE ${FFTW_CFLAGS_OTHER})

# Example executable
add_executable(main examples/main.cpp)
target_link_libraries(main moyal_solver)

# Installation
install(TARGETS moyal_solver main
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

install(DIRECTORY include/
        DESTINATION include/moyal_solver
        FILES_MATCHING PATTERN "*.hpp")

# Optional: Add testing
enable_testing()

# Compiler flags
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(moyal_solver PRIVATE -Wall -Wextra -O3)
endif()

# Print configuration info
message(STATUS "FFTW3 found: ${FFTW_FOUND}")
message(STATUS "FFTW3 libraries: ${FFTW_LIBRARIES}")
message(STATUS "FFTW3 include dirs: ${FFTW_INCLUDE_DIRS}")
