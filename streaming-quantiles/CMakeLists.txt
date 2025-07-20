cmake_minimum_required(VERSION 3.14)
project(StreamingQuantiles)

set(CMAKE_CXX_STANDARD 17)

add_executable(streaming-quantiles main.cpp)
target_link_libraries(streaming-quantiles PRIVATE fmt::fmt)
target_include_directories(streaming-quantiles PRIVATE ${cereal_SOURCE_DIR}/include)

# Google Test
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)

FetchContent_Declare(
  fmt
  URL https://github.com/fmtlib/fmt/archive/refs/tags/10.2.1.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
FetchContent_MakeAvailable(fmt)
message(STATUS "cereal_SOURCE_DIR = ${cereal_SOURCE_DIR}")

enable_testing()

add_executable(
  compactor_test
  compactor_test.cpp
)
target_link_libraries(
  compactor_test
  GTest::gtest_main
  fmt::fmt
)

include(GoogleTest)
gtest_discover_tests(compactor_test)

