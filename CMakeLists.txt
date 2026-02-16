cmake_minimum_required(VERSION 3.15)
project(Atoms VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# required packages
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(glfw3 REQUIRED)
find_package(glm REQUIRED)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(atom_2d src/wave_atom_2d.cpp)
target_link_libraries(atom_2d
    OpenGL::GL
    GLEW::GLEW
    glfw
    glm::glm
)

add_executable(atom_realtime src/atom_realtime.cpp)
target_link_libraries(atom_realtime
    OpenGL::GL
    GLEW::GLEW
    glfw
    glm::glm
)

add_executable(atom_raytracer src/atom_raytracer.cpp)
target_link_libraries(atom_raytracer
    OpenGL::GL
    GLEW::GLEW
    glfw
    glm::glm
)

add_executable(atom src/atom.cpp)
target_link_libraries(atom
    OpenGL::GL
    GLEW::GLEW
    glfw
    glm::glm
)

message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Output directory: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
