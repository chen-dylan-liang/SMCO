cmake_minimum_required(VERSION 3.14)
project(SMCO)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Check if the submodule is actually populated to prevent confusing build errors
if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/eigen-5.0.0/CMakeLists.txt")
    message(FATAL_ERROR "The 'eigen' folder appears empty. Please run: git submodule update --init --recursive")
endif()

# Add the local Eigen folder to the build.
add_subdirectory(eigen-5.0.0)

add_executable(SMCO SMCO.cpp)

# Link against the local Eigen target.
target_link_libraries(SMCO PRIVATE Eigen3::Eigen)

# Include the current directory (optional, for local headers)
target_include_directories(SMCO PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})