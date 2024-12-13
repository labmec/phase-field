cmake_minimum_required(VERSION 3.14)

project (PhaseField LANGUAGES C CXX)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# Finding the pz package
find_package(NeoPZ REQUIRED HINTS ../NeoPZ_install ../neopz_install)

file(GLOB_RECURSE PZ_SRC
     "../neopz/*.h"
     "../neopz/*.cpp"
     "../NeoPzDevelop/*.h"
     "../NeoPzDevelop/*.cpp"
)
add_custom_target( GlobalFiles SOURCES ${PZ_SRC})

# Creating the static library
add_library(phasefield_lib STATIC TPZPhaseField.cpp TPZElasticityPhaseField.cpp TPZPhaseFieldAnalysis.cpp)
target_link_libraries(phasefield_lib PUBLIC NeoPZ::pz)

# creating the projects
add_executable(phasefield phasefield.cpp)
target_link_libraries(phasefield PUBLIC phasefield_lib)

