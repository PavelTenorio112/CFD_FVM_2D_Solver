cmake_minimum_required(VERSION 3.16)
project(Geometric_Preprocessing VERSION 0.1.0 LANGUAGES CXX)



set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)



add_executable(Geometric_Preprocessor_Executable
    src/geometric_preprocess/mesh_reader.cpp
    src/geometric_preprocess/mesh_geometric_preprocess.cpp
    src/geometric_preprocess/mesh_data_print.cpp
    src/setup/boundary_conditions.cpp
    src/setup/initial_conditions.cpp
    src/setup/simulation_parameters_reader.cpp
    src/fields_data_print.cpp
    src/inviscid_burgers.cpp
    src/main.cpp
    src/upwind_interpolation.cpp
)

target_include_directories(Geometric_Preprocessor_Executable PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}

)

if(MSVC)
    target_compile_options(Geometric_Preprocessor_Executable PRIVATE /W4)
else()
    target_compile_options(Geometric_Preprocessor_Executable PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()



if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(Geometric_Preprocessor_Executable PRIVATE 
            -fsanitize=address,undefined 
            -fno-omit-frame-pointer 
            -g
        )
        target_link_options(Geometric_Preprocessor_Executable PRIVATE 
            -fsanitize=address,undefined
        )
    endif()
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(MSVC)
        target_compile_options(Geometric_Preprocessor_Executable PRIVATE /O2 /Oi /Ot)
    else()
        target_compile_options(Geometric_Preprocessor_Executable PRIVATE -O3 -march=native)
    endif()
endif()