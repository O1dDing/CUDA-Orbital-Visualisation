cmake_minimum_required(VERSION 3.28)

project(CUDAOrbitalVisualisation VERSION 0.1.0 LANGUAGES CXX)

option(COV_ENABLE_CUDA "Build the CUDA/OpenGL viewer" ON)
option(COV_BUILD_TESTS "Build parser tests" ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(cov_parser
    src/parser/molden_parser.cpp
)
target_include_directories(cov_parser PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_compile_features(cov_parser PUBLIC cxx_std_20)

if(MSVC)
    target_compile_options(cov_parser PRIVATE /W4 /permissive-)
else()
    target_compile_options(cov_parser PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(COV_ENABLE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(NOT CMAKE_CUDA_COMPILER)
        message(FATAL_ERROR "COV_ENABLE_CUDA=ON but no CUDA compiler was found.")
    endif()

    # RTX 5090 / Blackwell first: native sm_120 cubin plus compute_120 PTX.
    # Override at configure time for another CUDA GPU, e.g.
    # -DCMAKE_CUDA_ARCHITECTURES=89
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES "120-real;120-virtual" CACHE STRING
            "CUDA architectures (default: Blackwell sm_120 + PTX)" FORCE)
    endif()

    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 20)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

    find_package(CUDAToolkit 12.8 REQUIRED)
    find_package(OpenGL REQUIRED)

    include(FetchContent)

    set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
        glfw
        GIT_REPOSITORY https://github.com/glfw/glfw.git
        GIT_TAG 3.4
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(glfw)

    FetchContent_Declare(
        imgui
        GIT_REPOSITORY https://github.com/ocornut/imgui.git
        GIT_TAG v1.90.9
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(imgui)

    add_library(cov_imgui STATIC
        ${imgui_SOURCE_DIR}/imgui.cpp
        ${imgui_SOURCE_DIR}/imgui_draw.cpp
        ${imgui_SOURCE_DIR}/imgui_tables.cpp
        ${imgui_SOURCE_DIR}/imgui_widgets.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl2.cpp
    )
    target_include_directories(cov_imgui PUBLIC
        ${imgui_SOURCE_DIR}
        ${imgui_SOURCE_DIR}/backends
    )
    target_link_libraries(cov_imgui PUBLIC glfw OpenGL::GL)

    add_library(cov_gpu
        src/cuda/orbital.cu
        src/render/gl_api.cpp
        src/render/volume_renderer.cpp
    )
    target_include_directories(cov_gpu PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    )
    target_link_libraries(cov_gpu
        PUBLIC
            cov_parser
            CUDA::cudart
            OpenGL::GL
            glfw
    )
    set_target_properties(cov_gpu PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        POSITION_INDEPENDENT_CODE ON
    )
    target_compile_options(cov_gpu PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--use_fast_math>
    )

    add_executable(cov
        src/main.cpp
    )
    target_link_libraries(cov PRIVATE
        cov_gpu
        cov_imgui
    )

    if(WIN32)
        set_property(TARGET cov PROPERTY WIN32_EXECUTABLE FALSE)
    endif()
endif()

if(COV_BUILD_TESTS)
    include(CTest)
    enable_testing()

    add_executable(cov_parser_smoke
        tests/parser_smoke.cpp
    )
    target_link_libraries(cov_parser_smoke PRIVATE cov_parser)

    add_test(
        NAME parser_smoke
        COMMAND cov_parser_smoke ${CMAKE_CURRENT_SOURCE_DIR}/examples/h2.molden
    )
endif()
