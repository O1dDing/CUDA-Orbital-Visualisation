cmake_minimum_required(VERSION 3.28)

project(CUDAOrbitalVisualisation VERSION 0.2.0 LANGUAGES CXX)

option(COV_ENABLE_CUDA "Build the CUDA/OpenGL viewer" ON)
option(COV_BUILD_TESTS "Build parser, orbital, diagram, molecule and UI smoke tests" ON)

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

add_library(cov_orbitals
    src/orbitals/orbital_view.cpp
    src/orbitals/mo_diagram.cpp
)
target_include_directories(cov_orbitals PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_compile_features(cov_orbitals PUBLIC cxx_std_20)

add_library(cov_molecule
    src/render/molecule_style.cpp
)
target_include_directories(cov_molecule PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_compile_features(cov_molecule PUBLIC cxx_std_20)

add_library(cov_platform
    src/platform/file_dialog.cpp
)
target_include_directories(cov_platform PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_compile_features(cov_platform PUBLIC cxx_std_20)
if(WIN32)
    target_link_libraries(cov_platform PRIVATE comdlg32)
endif()

if(MSVC)
    target_compile_options(cov_parser PRIVATE /W4 /permissive- /utf-8)
    target_compile_options(cov_orbitals PRIVATE /W4 /permissive- /utf-8)
    target_compile_options(cov_molecule PRIVATE /W4 /permissive- /utf-8)
    target_compile_options(cov_platform PRIVATE /W4 /permissive- /utf-8)
else()
    target_compile_options(cov_parser PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(cov_orbitals PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(cov_molecule PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(cov_platform PRIVATE -Wall -Wextra -Wpedantic)
endif()

# The localisation/theme layer is also compiled in CPU-only CI, so UI source
# compatibility does not depend on having a CUDA runner available.
if(COV_ENABLE_CUDA OR COV_BUILD_TESTS)
    include(FetchContent)
    FetchContent_Declare(
        imgui
        URL https://codeload.github.com/ocornut/imgui/zip/refs/tags/v1.90.9
        URL_HASH SHA256=4fa247c7abbb00fc255631479004de549b86a41f49fddd7d44cd2cf8d647c5a7
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable(imgui)

    add_library(cov_imgui_core STATIC
        ${imgui_SOURCE_DIR}/imgui.cpp
        ${imgui_SOURCE_DIR}/imgui_draw.cpp
        ${imgui_SOURCE_DIR}/imgui_tables.cpp
        ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    )
    target_include_directories(cov_imgui_core PUBLIC
        ${imgui_SOURCE_DIR}
    )

    add_library(cov_ui STATIC
        src/ui/ui.cpp
        src/ui/orbital_ui.cpp
    )
    target_include_directories(cov_ui PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${imgui_SOURCE_DIR}
    )
    target_link_libraries(cov_ui PUBLIC cov_imgui_core cov_orbitals cov_molecule)
    target_compile_features(cov_ui PUBLIC cxx_std_20)
    target_compile_options(cov_ui PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>
    )
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

    set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
        glfw
        URL https://codeload.github.com/glfw/glfw/zip/refs/tags/3.4
        URL_HASH SHA256=a133ddc3d3c66143eba9035621db8e0bcf34dba1ee9514a9e23e96afd39fd57a
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable(glfw)

    add_library(cov_imgui_backends STATIC
        ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl2.cpp
    )
    target_include_directories(cov_imgui_backends PUBLIC
        ${imgui_SOURCE_DIR}
        ${imgui_SOURCE_DIR}/backends
    )
    target_link_libraries(cov_imgui_backends PUBLIC
        cov_imgui_core
        glfw
        OpenGL::GL
    )

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
            cov_molecule
            CUDA::cudart
            OpenGL::GL
            glfw
    )
    set_target_properties(cov_gpu PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON
        POSITION_INDEPENDENT_CODE ON
    )
    target_compile_options(cov_gpu PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--use_fast_math>
    )

    add_executable(cov
        src/main.cpp
    )
    target_include_directories(cov PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    )
    target_link_libraries(cov PRIVATE
        cov_gpu
        cov_orbitals
        cov_molecule
        cov_platform
        cov_ui
        cov_imgui_backends
    )

    target_compile_options(cov PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>
    )

    if(WIN32)
        set_property(TARGET cov PROPERTY WIN32_EXECUTABLE FALSE)
    endif()
endif()

if(COV_BUILD_TESTS)
    include(CTest)
    enable_testing()

    add_executable(cov_parser_smoke tests/parser_smoke.cpp)
    target_link_libraries(cov_parser_smoke PRIVATE cov_parser)
    add_test(NAME parser_smoke
        COMMAND cov_parser_smoke ${CMAKE_CURRENT_SOURCE_DIR}/examples/h2.molden)

    add_executable(cov_orbital_view_smoke tests/orbital_view_smoke.cpp)
    target_link_libraries(cov_orbital_view_smoke PRIVATE cov_orbitals)
    add_test(NAME orbital_view_smoke COMMAND cov_orbital_view_smoke)

    add_executable(cov_mo_diagram_smoke tests/mo_diagram_smoke.cpp)
    target_link_libraries(cov_mo_diagram_smoke PRIVATE cov_orbitals)
    add_test(NAME mo_diagram_smoke COMMAND cov_mo_diagram_smoke)

    add_executable(cov_molecule_style_smoke tests/molecule_style_smoke.cpp)
    target_link_libraries(cov_molecule_style_smoke PRIVATE cov_molecule)
    add_test(NAME molecule_style_smoke COMMAND cov_molecule_style_smoke)

    add_executable(cov_ui_smoke tests/ui_smoke.cpp)
    target_link_libraries(cov_ui_smoke PRIVATE cov_ui)
    target_compile_options(cov_ui_smoke PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>
    )
    add_test(NAME ui_smoke COMMAND cov_ui_smoke)
endif()
