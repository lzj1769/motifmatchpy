cmake_minimum_required(VERSION 3.15...3.31)
project(moods_reference LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(Python 3.12
  REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)

# Locate nanobind's CMake package via the Python module installed in the build env.
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

# Vendored MOODS 1.9.4.1 core (github.com/jhkorhonen/MOODS), see core/VENDORED.md
set(MOODS_CORE_SOURCES
  core/match_types.cpp
  core/moods_misc.cpp
  core/moods_parsers.cpp
  core/moods_scan.cpp
  core/moods_tools.cpp
  core/motif_0.cpp
  core/motif_h.cpp
  core/scanner.cpp)

nanobind_add_module(_core
  NB_STATIC
  STABLE_ABI
  LTO
  bindings.cpp
  ${MOODS_CORE_SOURCES})

target_include_directories(_core PRIVATE core)

# Upstream MOODS builds with -march=native, which produces wheels that crash on
# older CPUs. Opt in explicitly with -DMOODS_REFERENCE_NATIVE=ON for local builds.
option(MOODS_REFERENCE_NATIVE "Optimise for the building machine's CPU" OFF)
if(MOODS_REFERENCE_NATIVE AND NOT MSVC)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-march=native" HAVE_MARCH_NATIVE)
  if(HAVE_MARCH_NATIVE)
    target_compile_options(_core PRIVATE -march=native)
  endif()
endif()

# Compile with strict IEEE semantics. By default the compiler may contract
# `a + b * c` into a fused multiply-add, which is more accurate but makes the
# result depend on the compiler and target -- and puts it a ULP away from what
# Python computes. Turning contraction off is what lets the test suite assert
# bitwise equality between motifmatchpy and this reference rather than
# "close enough".
if(MSVC)
  target_compile_options(_core PRIVATE /O2 /fp:precise)
else()
  target_compile_options(_core PRIVATE -O3 -ffp-contract=off)
endif()

install(TARGETS _core LIBRARY DESTINATION moods_reference)
