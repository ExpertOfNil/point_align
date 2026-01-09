cmake_minimum_required(VERSION 3.16)
project(point_align LANGUAGES C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_C_STANDARD 99)
set(EXTRA_WARNING_FLAGS -Wall -Wextra -Wpedantic)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build" FORCE)
    message(WARNING "Defaulting to Debug build")
endif()

message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    # Use defaults
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    # Use defaults
else()
    if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        message(WARNING "Unknown build type ${CMAKE_BUILD_TYPE}. Defaulting to Debug.")
    endif()
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build" FORCE)
    #list(APPEND EXTRA_COMPILE_FLAGS -DDEBUG -fno-omit-frame-pointer -fsanitize=address,undefined)
    list(APPEND EXTRA_COMPILE_FLAGS -DDEBUG)
endif()
add_compile_options(${EXTRA_COMPILE_FLAGS})
string(JOIN " " EXTRA_COMPILE_FLAGS_STR ${EXTRA_COMPILE_FLAGS})
message(STATUS "Flags: ${CMAKE_C_FLAGS_DEBUG} ${EXTRA_COMPILE_FLAGS_STR}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

set(POINT_ALIGN_SOURCES
    src/point_align.c
)

option(USE_OPENBLAS "Use OpenBLAS for SVD" ON)

add_library(point_align SHARED ${POINT_ALIGN_SOURCES})
set_target_properties(point_align PROPERTIES
    VERSION 0.1.0
    SOVERSION 0
)
target_include_directories(point_align PUBLIC
    ${PROJECT_SOURCE_DIR}/include/
    ${PROJECT_SOURCE_DIR}/vendor/include/
)

add_library(point_align_static STATIC ${POINT_ALIGN_SOURCES})
set_target_properties(point_align_static PROPERTIES
    OUTPUT_NAME point_align
    POSITION_INDEPENDENT_CODE ON
)
target_include_directories(point_align_static PUBLIC
    ${PROJECT_SOURCE_DIR}/include/
    ${PROJECT_SOURCE_DIR}/vendor/include/
)

if(USE_OPENBLAS AND EXISTS "${PROJECT_SOURCE_DIR}/vendor/lib/openblas/libopenblas.a")
    message(STATUS "OpenBLAS: ENABLED")
    add_library(openblas STATIC IMPORTED)
    set_target_properties(openblas PROPERTIES
        IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/vendor/lib/openblas/libopenblas.a"
        INTERFACE_INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}/vendor/include/openblas/"
    )
    target_link_libraries(point_align PRIVATE
        openblas
        pthread
        m
    )
    target_link_libraries(point_align_static PUBLIC
        openblas
        pthread
        m
    )
else()
    if(USE_OPENBLAS)
        message(WARNING "OpenBLAS: DISABLED (libopenblas.a not found)")
    else()
        message(STATUS "OpenBLAS: DISABLED (USE_OPENBLAS=OFF)")
    endif()
    target_link_libraries(point_align PRIVATE m)
    target_link_libraries(point_align_static PUBLIC m)
endif()

# Build executables

file(GLOB EXAMPLE_SOURCES "examples/*.c")
foreach(EXAMPLE_SOURCE ${EXAMPLE_SOURCES})
    get_filename_component(EXAMPLE_NAME ${EXAMPLE_SOURCE} NAME_WE)
    add_executable(${EXAMPLE_NAME} ${EXAMPLE_SOURCE})
    set_target_properties(${EXAMPLE_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples/"
    )
    target_link_libraries(${EXAMPLE_NAME}
        point_align
        m
    )
    target_include_directories(${EXAMPLE_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include/)
endforeach()

file(GLOB TEST_SOURCES "tests/*.c")
foreach(TEST_SOURCE ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_SOURCE})
    set_target_properties(${TEST_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tests/"
    )
    target_link_libraries(${TEST_NAME}
        point_align
        m
    )
    target_include_directories(${TEST_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include/)
endforeach()
