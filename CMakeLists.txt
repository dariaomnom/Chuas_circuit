cmake_minimum_required(VERSION 3.10)
project(Chuas_circuit)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp)
add_executable(Chuas_circuit ${SOURCE_FILES})
include_directories(/usr/local/include)

find_package(SFML 2.6 COMPONENTS system window graphics network audio REQUIRED)
include_directories(${SFML_INCLUDE_DIRS})
target_link_libraries(Chuas_circuit sfml-system sfml-window sfml-graphics sfml-audio sfml-network)