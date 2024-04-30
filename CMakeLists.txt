cmake_minimum_required(VERSION 3.25)
project(convex_hull)

set(CMAKE_VERBOSE_MAKEFILE TRUE)

set(CMAKE_C_STANDARD 99)

add_executable(convex_hull main.c)
target_link_libraries(convex_hull m)

add_executable(convex_hull_stress_test main.c)
target_link_libraries(convex_hull_stress_test m)
target_compile_definitions(convex_hull_stress_test PUBLIC UNITTEST=true)
