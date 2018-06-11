cmake_minimum_required(VERSION 2.6)
project(timedil_harness)

add_executable(timedil_harness main.cpp)

install(TARGETS timedil_harness RUNTIME DESTINATION bin)
