cmake_minimum_required(VERSION 3.14)
project(StreamingQuantiles)

set(CMAKE_CXX_STANDARD 23)

find_package(fmt)

add_executable(streaming-quantiles main.cpp)
target_link_libraries(streaming-quantiles PRIVATE fmt::fmt)

