cmake_minimum_required(VERSION 3.10)

project(CB-EVO VERSION 1.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc -static-libstdc++")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")

find_package(OpenMP REQUIRED)

# Add other liberaries
add_library(fmt INTERFACE)
target_include_directories(fmt SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/lib/fmt)
target_compile_definitions(fmt INTERFACE FMT_HEADER_ONLY)

add_library(json INTERFACE)
target_include_directories(json SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/lib/json)

add_library(cli11 INTERFACE)
target_include_directories(cli11 SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/lib/cli11)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen-3.4.0)

file(GLOB OTHER_SOURCES "src/*.cpp")
list(REMOVE_ITEM OTHER_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
# Add the executable based on the main source file and the other source files.
add_executable(cb-evo src/main.cpp ${OTHER_SOURCES})
target_link_libraries(cb-evo PRIVATE stdc++fs dl cli11 readline pthread json fmt OpenMP::OpenMP_CXX)
