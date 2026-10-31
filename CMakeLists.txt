cmake_minimum_required(VERSION 3.16)
project(flac_gqh LANGUAGES CXX)

# The portable tests exercise the shared constitutive code, not a FLAC3D DLL.
add_library(gqh_refcurve STATIC src/refcurve.cpp)
target_include_directories(gqh_refcurve PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_compile_features(gqh_refcurve PUBLIC cxx_std_17)

include(CTest)
if(BUILD_TESTING)
  add_executable(gqh_core_tests tests/test_gqh.cpp)
  target_link_libraries(gqh_core_tests PRIVATE gqh_refcurve)
  add_test(NAME gqh_core COMMAND gqh_core_tests)
endif()
