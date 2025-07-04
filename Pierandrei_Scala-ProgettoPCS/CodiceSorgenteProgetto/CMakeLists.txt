cmake_minimum_required(VERSION 3.9)

project(Poliedro LANGUAGES CXX C VERSION 1.0.0)

# Set C++ Standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Impose warnings in Debug
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic-errors")

# CXX Flags for Windows
if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif()

# Dependencies
find_package(Eigen3 CONFIG REQUIRED)
find_package(Threads REQUIRED)
find_package(GTest REQUIRED)

# Source subdirectories
add_subdirectory(src)
add_subdirectory(ExportParaview)
add_subdirectory(src_test)

# Collect sources/headers/includes da src e ExportParaview
list(APPEND poliedro_sources ${poliedro_sources_src})
list(APPEND poliedro_headers ${poliedro_headers_src})
list(APPEND poliedro_includes ${poliedro_includes_src})

list(APPEND poliedro_sources ${poliedro_sources_paraview})
list(APPEND poliedro_headers ${poliedro_headers_paraview})
list(APPEND poliedro_includes ${poliedro_includes_paraview})

# Executable: main
add_executable(Poliedro
    main.cpp
    ${poliedro_sources}
    ${poliedro_headers}
)

target_include_directories(Poliedro PRIVATE ${poliedro_includes})
target_link_libraries(Poliedro
    PRIVATE
    Eigen3::Eigen
    Threads::Threads
)

target_compile_options(Poliedro PUBLIC -fPIC)

# Executable: test
add_executable(PoliedroTest
    main_test.cpp
    ${poliedro_sources}
    ${poliedro_headers}
    ${poliedro_sources_test}
    ${poliedro_headers_test}
)

target_include_directories(PoliedroTest
    PRIVATE
    ${poliedro_includes}
    ${poliedro_includes_test}
)

target_link_libraries(PoliedroTest
    PRIVATE
    GTest::gtest
    GTest::gtest_main
    GTest::gmock
    GTest::gmock_main
    Threads::Threads
    Eigen3::Eigen
)

