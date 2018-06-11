cmake_minimum_required(VERSION 3.1)
project(timedil_harness)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(timedil_harness main.cpp)
target_link_libraries(timedil_harness Threads::Threads)

install(TARGETS timedil_harness RUNTIME DESTINATION bin)
