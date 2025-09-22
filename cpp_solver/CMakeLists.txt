cmake_minimum_required(VERSION 3.15)
project(polc_darcy_solver LANGUAGES CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Helpful for VS Code (clangd, IntelliSense)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build type default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Warnings (no -Werror yet; don't block progress)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  add_compile_options(-Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion)
elseif (MSVC)
  add_compile_options(/W4)
endif()

# Optional sanitizers (enable when you want)
option(POLC_USE_SANITIZERS "Enable Address/UB sanitizers in Debug" ON)
if (POLC_USE_SANITIZERS AND CMAKE_BUILD_TYPE MATCHES "Debug")
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-fsanitize=address,undefined)
    add_link_options(-fsanitize=address,undefined)
  endif()
endif()

# Sources
add_executable(server
  src/server.cpp
)

target_include_directories(server PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Placeholder for gRPC/protobuf later
# add_subdirectory(third_party/grpc ...)  or find_package(gRPC CONFIG REQUIRED)
# target_link_libraries(server PRIVATE gRPC::grpc++ protobuf::libprotobuf)
