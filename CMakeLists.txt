cmake_minimum_required(VERSION 4.0)
message(STATUS "Using CMake version ${CMAKE_VERSION}")
cmake_policy(VERSION 4.0)

project(EWF LANGUAGES CXX VERSION 2.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_SHARED_LIBS "Enable compilations of shared libraries" OFF)

option(PYTHON_BINDINGS "Whether to build python bindings" ON)
if (PYTHON_BINDINGS)
    include(FetchContent)
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG        v2.12.1       # or v2.11.1, etc.
    )
    FetchContent_GetProperties(pybind11)
    if (NOT pybind11_POPULATED)
        FetchContent_MakeAvailable(pybind11)
    endif()
endif()

add_subdirectory(src)

#Clang format target
find_program(CLANG_FORMAT_EXE NAMES "clang-format-16" "clang-format-15" "clang-format-14" "clang-format-13" "clang-format-12"
"clang-format-11" "clang-format-10" "clang-format" DOC "Path to clang-format executable"
)
if (NOT CLANG_FORMAT_EXE)
message(STATUS "clang-format not found.")
else ()
message(STATUS "clang-format found: ${CLANG_FORMAT_EXE}")
add_custom_target(
    EWF_clang_format
    COMMAND ${CLANG_FORMAT_EXE} -style=file -i ${CMAKE_SOURCE_DIR}/src/*pp
    COMMAND ${CLANG_FORMAT_EXE} -style=file -i ${CMAKE_SOURCE_DIR}/example/*pp
    COMMAND ${CLANG_FORMAT_EXE} -style=file -i ${CMAKE_SOURCE_DIR}/test/*pp
)
endif ()