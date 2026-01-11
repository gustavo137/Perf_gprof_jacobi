cmake_minimum_required(VERSION 3.16)

project(jacobi_perf LANGUAGES CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optional: warnings
add_compile_options(-Wall -Wextra -Wpedantic)

# Where to put executables: build/bin/
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# MPI
find_package(MPI REQUIRED)

# Executable
add_executable(jacobi.x
    JACOBI/src/jacobi.cpp
)

# Include headers from include/
target_include_directories(jacobi.x
  PRIVATE
  ${PROJECT_SOURCE_DIR}/JACOBI/include
)

# Include MPI headers and link MPI::MPI_CXX
target_link_libraries(jacobi.x 
    PRIVATE
      MPI::MPI_CXX
)

# --- gprof (optional) recomenden with gprof2dot ---
# Activate gprof with: cmake -DENABLE_GPROF=ON ..
option(ENABLE_GPROF "Enable gprof profiling flags (-pg)" OFF)
if(ENABLE_GPROF)
  target_compile_options(jacobi.x PRIVATE -g -pg)
  target_link_options(jacobi.x PRIVATE -g -pg)
endif()

# --- perf (optional) (best option is using -00 or -01) ---
# Activate perf flags with: cmake -DENABLE_PERF=ON ..
option(ENABLE_PERF "Build with flags useful for perf + call stacks" OFF)
if(ENABLE_PERF)
  # target_compile_options(jacobi.x PRIVATE -O3 -g -fno-omit-frame-pointer)
  # If you prefer even more stable stacks for C++:
  target_compile_options(jacobi.x PRIVATE -O2 -g -fno-omit-frame-pointer)
endif()

