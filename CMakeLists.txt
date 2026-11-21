cmake_minimum_required(VERSION 3.10)
# Keep VERSION in sync with xfoil_c/__init__.py (tests/test_skeleton.py checks this).
project(xfoil_c VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
# Export only symbols marked XFOIL_API. Use double precision everywhere,
# matching the Fortran build (-fdefault-real-8).
set(CMAKE_C_VISIBILITY_PRESET hidden)
# =====================================================================================================
#  BUILD VARIANTS: how to choose the optimizations
# =====================================================================================================
#  Two independent switches, set at configure time. Leaving both unset gives the EXACT build.
#
#  XFOIL_ARCH   (string, default "")   CPU target
#      ""          portable x86-64 code, no FMA (-ffp-contract=off): results are bit-identical to the
#                  Fortran reference XFOIL. Runs on any x86-64 CPU.            -> part of the EXACT build
#      native      -march=native: all instructions of THIS machine's CPU (AVX2/AVX-512, FMA) and
#                  a*b+c fused into FMA (-ffp-contract=fast). ~10-15 % faster; results differ from the
#                  reference in the last digits (measured: <= 1.2e-5 relative). The library may crash
#                  (illegal instruction) on an older CPU: rebuild it on each machine.
#      <cpu>       any gcc -march value, e.g. x86-64-v3 (AVX2 + FMA, Intel Haswell / AMD Zen and newer),
#                  znver5, ...: as `native`, but for a named CPU family (portable to all CPUs of it).
#
#  XFOIL_FAST   (ON/OFF, default OFF)   source-level fast paths (code in #if XFOIL_FAST ... #else exact)
#      O2  blsolv: the BL Newton-system blocks VM stored row-contiguous (src/xlayout.h), so the
#          elimination vectorises (blsolv 0.52 -> 0.15 ms per call, ~3.5x; a viscous polar ~16 %
#          faster). Same arithmetic: on its own bit-identical to the exact build.
#      Further fast paths are listed in MIGRATION_PLAN.md (optimization phase); each is documented
#      where it is implemented.
#
#  Combined (measured, viscous NACA 0012 aseq, 80 points):   exact 1.57 s | FAST 1.30 s |
#      ARCH=native 1.38 s | ARCH=native + FAST 1.09 s (-30 %)
#
#  HOW TO BUILD (from xfoil-c/; setup.py always passes both switches, so an old CMake cache never
#  keeps a previous choice; the -D form works as well, e.g. -DXFOIL_ARCH=native):
#      python setup.py build_ext --inplace                                    # exact (default)
#      XFOIL_FAST=ON python setup.py build_ext --inplace                      # fast paths, portable
#      XFOIL_ARCH=native XFOIL_FAST=ON python setup.py build_ext --inplace    # fastest, this machine
#      XFOIL_ARCH=native XFOIL_FAST=ON pip install .                          # same, installed
#  Plain CMake:  cmake -S . -B build -DXFOIL_ARCH=native -DXFOIL_FAST=ON && cmake --build build
#
#  CHECK WHICH VARIANT IS LOADED:   python -c "import xfoil_c; print(xfoil_c.build_info())"
#      -> "exact"  or  "arch=native fast=ON" (compiled in as XFOIL_BUILD_INFO below)
#
#  TESTS:  python -m pytest tests -q
#      exact build: the full suite, bit-identical to the Fortran library.
#      other builds: the bit-identical tests are skipped automatically; tests/test_reference.py checks
#      the variant against the reference results (deviation <= 1e-3, no more unconverged points).
#      Rebuild the exact variant before committing code changes and run the full suite on it.
#
#  ALWAYS ON (all variants): -O3 -fno-math-errno (bit-identical, ~5 % faster than -O2; clean-up C7).
# =====================================================================================================
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fno-math-errno -Wall -Wextra")

set(XFOIL_ARCH "" CACHE STRING "CPU target (-march); empty = portable and exact. See the header above.")
option(XFOIL_FAST "Source-level fast paths (#if XFOIL_FAST). See the header above." OFF)
if(XFOIL_ARCH STREQUAL "")
    # never fuse a*b+c into an FMA: results stay bit-identical to the Fortran build on every CPU
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffp-contract=off")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=${XFOIL_ARCH} -ffp-contract=fast")
endif()
set(XFOIL_BUILD_INFO "exact")
if(NOT XFOIL_ARCH STREQUAL "" OR XFOIL_FAST)
    set(XFOIL_BUILD_INFO "arch=${XFOIL_ARCH} fast=${XFOIL_FAST}")
endif()
message(STATUS "xfoil_c build variant: ${XFOIL_BUILD_INFO}")

add_library(xfoil_c SHARED
        src/xstate.c
        src/xglobals.c
        src/xerror.c
        src/xspline.c
        src/xutils.c
        src/xuserio.c
        src/xnaca.c
        src/xgeom.c
        src/xsolve.c
        src/sxbl.c
        src/xpanel.c
        src/xblsys.c
        src/xbl.c
        src/sxfoil.c
        src/sxoper.c
        src/xoper.c
        src/xfoil.c
        src/xgdes.c
        src/xfoil_ext.c
        src/xfoil_api.c)

target_include_directories(xfoil_c PRIVATE src)
target_compile_definitions(xfoil_c PRIVATE XFOIL_C_VERSION="${PROJECT_VERSION}"
                           XFOIL_BUILD_INFO="${XFOIL_BUILD_INFO}")
if(XFOIL_FAST)
    target_compile_definitions(xfoil_c PRIVATE XFOIL_FAST=1)
endif()
target_link_libraries(xfoil_c PRIVATE m)
if(MINGW)
    # Windows (MinGW-w64 gcc): link libgcc (and winpthread, if pulled in) statically, so the DLL
    # needs only system DLLs (KERNEL32 + the UCRT) and loads from a python.org Python as is.
    target_link_options(xfoil_c PRIVATE -static-libgcc -static)
    if(NOT XFOIL_ARCH STREQUAL "")
        # With AVX/AVX-512 gcc may spill vectors to the stack with aligned moves, but the Windows ABI
        # only aligns the stack to 16 bytes (gcc bug 54412): have the assembler emit unaligned moves
        # instead (same speed on AVX CPUs, same results).
        target_compile_options(xfoil_c PRIVATE -Wa,-muse-unaligned-vector-move)
    endif()
endif()
