# The top-level CMake file is there to bring all modules into scope. That
# means, adding the subdirectories for all CMake projects in this tree, and
# finding external libraries and turning them into imported targets.

cmake_minimum_required(VERSION 3.15)

# set preference for clang compiler and intel compiler over gcc and other compilers
include(Platform/${CMAKE_SYSTEM_NAME}-Determine-C OPTIONAL)
include(Platform/${CMAKE_SYSTEM_NAME}-C OPTIONAL)
set(CMAKE_C_COMPILER_NAMES clang icc cc ${CMAKE_C_COMPILER_NAMES})

include(Platform/${CMAKE_SYSTEM_NAME}-Determine-CXX OPTIONAL)
include(Platform/${CMAKE_SYSTEM_NAME}-CXX OPTIONAL)
set(CMAKE_CXX_COMPILER_NAMES clang++ icpc c++ ${CMAKE_CXX_COMPILER_NAMES})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(utils-highs)
set_version(VERSION)

project(HIGHS VERSION ${VERSION} LANGUAGES CXX C)

set(PROJECT_NAMESPACE highs)
message(STATUS "${PROJECT_NAME} version: ${PROJECT_VERSION}")

# use C++11 standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default Build Type to be Release
get_property(isMultiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(isMultiConfig)
  if(NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_CONFIGURATION_TYPES "Release;Debug" CACHE STRING
    "Choose the type of builds, options are: Debug Release RelWithDebInfo MinSizeRel. (default: Release;Debug)"
    FORCE)
  endif()
  message(STATUS "Configuration types: ${CMAKE_CONFIGURATION_TYPES}")
else()
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel. (default: Release)"
    FORCE)
  endif()
  message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
endif()

if(CMAKE_BUILD_TYPE STREQUAL Release)
    set(HiGHSRELEASE ON)
endif()

# Layout build dir like install dir
include(GNUInstallDirs)

if(UNIX)
  option(BUILD_SHARED_LIBS "Build shared libraries (.so or .dyld)." ON)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
  # for multi-config build system (e.g. xcode)
  foreach(OUTPUTCONFIG IN LISTS CMAKE_CONFIGURATION_TYPES)
    string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/${CMAKE_INSTALL_LIBDIR})
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/${CMAKE_INSTALL_LIBDIR})
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/${CMAKE_INSTALL_BINDIR})
  endforeach()
else()
  # Currently Only support static build for windows
  option(BUILD_SHARED_LIBS "Build shared libraries (.dll)." OFF)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
  # for multi-config builds (e.g. msvc)
  foreach(OUTPUTCONFIG IN LISTS CMAKE_CONFIGURATION_TYPES)
    string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/${CMAKE_INSTALL_BINDIR})
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/${CMAKE_INSTALL_BINDIR})
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}/${CMAKE_INSTALL_BINDIR})
  endforeach()
endif()

# if(BUILD_SHARED_LIBS AND MSVC)
#   set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
# endif()

# option(BUILD_SHARED_LIBS "Build shared libraries." ON)

# todo ? 
# set(CMAKE_PLATFORM_USES_PATH_WHEN_NO_SONAME FALSE)


if(MSVC)
	# This option is only available when building with MSVC. By default, highs
	# is build using the cdecl calling convention, which is useful if you're
	# writing C. However, the CLR and Win32 API both expect stdcall.
	#
	# If you are writing a CLR program and want to link to highs, you'll want
	# to turn this on by invoking CMake with the "-DSTDCALL=ON" argument.
	option( STDCALL			"Build highs with the __stdcall convention"   OFF  )
endif()

# By default only build the C++ library.
option(BUILD_CXX "Build C++ library" ON)
message(STATUS "Build C++ library: ${BUILD_CXX}")

option(PYTHON "Build Python interface" OFF)
message(STATUS "Build Python: ${PYTHON}")
option(FORTRAN "Build Fortran interface" OFF)
message(STATUS "Build Fortran: ${FORTRAN}")
option(CSHARP "Build CSharp interface" OFF)
message(STATUS "Build CSharp: ${CSHARP}")

# ZLIB can be switched off, for building interfaces
option(ZLIB "Fast build: " ON)

# If wrapper are built, we need to have the install rpath in BINARY_DIR to package
if(PYTHON OR FORTRAN OR CSHARP)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
endif()

# For Python interface
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

set(CMAKE_MACOSX_RPATH ON)

# todo rpaht again
# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")


# Basic type
include(CMakePushCheckState)
cmake_push_check_state(RESET)
set(CMAKE_EXTRA_INCLUDE_FILES "cstdint")
include(CheckTypeSize)
check_type_size("long" SIZEOF_LONG LANGUAGE CXX)
message(STATUS "Found long size: ${SIZEOF_LONG}")
check_type_size("long long" SIZEOF_LONG_LONG LANGUAGE CXX)
message(STATUS "Found long long size: ${SIZEOF_LONG_LONG}")
check_type_size("int64_t" SIZEOF_INT64_T LANGUAGE CXX)
message(STATUS "Found int64_t size: ${SIZEOF_INT64_T}")

check_type_size("unsigned long" SIZEOF_ULONG LANGUAGE CXX)
message(STATUS "Found unsigned long size: ${SIZEOF_ULONG}")
check_type_size("unsigned long long" SIZEOF_ULONG_LONG LANGUAGE CXX)
message(STATUS "Found unsigned long long size: ${SIZEOF_ULONG_LONG}")
check_type_size("uint64_t" SIZEOF_UINT64_T LANGUAGE CXX)
message(STATUS "Found uint64_t size: ${SIZEOF_UINT64_T}")

check_type_size("int *" SIZEOF_INT_P LANGUAGE CXX)
message(STATUS "Found int * size: ${SIZEOF_INT_P}")
cmake_pop_check_state()


### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file).
    Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

option(BUILD_TESTING "Build Tests" ON)

if(DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION)
    message(STATUS "IPO / LTO as requested by user: ${CMAKE_INTERPROCEDURAL_OPTIMIZATION}")
elseif(LINUX AND (NOT MSVC) AND (NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
    include(CheckIPOSupported)
    check_ipo_supported(RESULT ipo_supported OUTPUT error)

    if(ipo_supported)
        message(STATUS "IPO / LTO enabled")
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    else()
        message(STATUS "IPO / LTO not supported: <${error}>")
    endif()
endif()


# emscripten
option(EMSCRIPTEN_HTML "Emscripten HTML output" OFF)


include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
    "#include <immintrin.h>
    int main () {
        _mm_pause();
        return 0;
    }"
    HIGHS_HAVE_MM_PAUSE)

if(MSVC)
    check_cxx_source_compiles(
        "#include <intrin.h>
        #pragma intrinsic(_BitScanReverse)
        #pragma intrinsic(_BitScanReverse64)
        int main () {
            unsigned long x = 5;
            unsigned long y;
            _BitScanReverse(&y, x);
            _BitScanReverse64(&x, y);
            return 0;
        }"
        HIGHS_HAVE_BITSCAN_REVERSE)
else()
    check_cxx_source_compiles(
        "#include <cstdint>
         int main () {
            unsigned int x = 5;
            unsigned long long y = __builtin_clz(x);
            x = __builtin_clzll(y);
            return 0;
        }"
        HIGHS_HAVE_BUILTIN_CLZ)
endif()

include(CheckCXXCompilerFlag)

if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86\_64|i686)")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpopcnt")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(ppc64|powerpc64)" AND NOT APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpopcntd")
else()
  message("FLAG_MPOPCNT_SUPPORTED is not available on this architecture")
endif()

option(DEBUGSOL "check the debug solution" OFF)

if(DEBUGSOL)
  add_definitions("-DHIGHS_DEBUGSOL")
endif()

option(HIGHSINT64 "Use 64 bit integers indexing" OFF)
if (NOT (${HIGHSINT64} STREQUAL  "OFF"))
    message(STATUS "HIGHSINT64: " ${HIGHSINT64})
endif()

# If Visual Studio targets are being built.
if(MSVC)
    add_definitions(/W4)
    add_definitions(/wd4018 /wd4061 /wd4100 /wd4101 /wd4127 /wd4189 /wd4244 /wd4245 /wd4267 /wd4324 /wd4365 /wd4389 /wd4456 /wd4457 /wd4458 /wd4459 /wd4514 /wd4701 /wd4820)
    add_definitions(/MP)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    if (STDCALL)
	# /Gz - stdcall calling convention
        add_definitions(/Gz)
    endif()

#
# This fouls up building HiGHS as part of a larger CMake project: see #1129.
# Setting it will speed up run-time in debug (and compile-time?)
#
#    add_definitions(-D_ITERATOR_DEBUG_LEVEL=0)
endif()

if (FORTRAN)
  include(CheckLanguage)
  if(NOT MSVC)
      check_language("Fortran")
  endif()
  if(CMAKE_Fortran_COMPILER)
      enable_language(Fortran)
      set(FORTRAN_FOUND ON)
  else()
      set(FORTRAN_FOUND OFF)
  endif(CMAKE_Fortran_COMPILER)
endif()

if (CSHARP)
  check_language("CSharp")
  if(CMAKE_CSharp_COMPILER)
      enable_language(CSharp)
      set(CSHARP_FOUND ON)
  else()
      set(CSHARP_FOUND OFF)
  endif(CMAKE_CSharp_COMPILER)
endif()

check_cxx_compiler_flag("-fno-omit-frame-pointer" NO_OMIT_FRAME_POINTER_FLAG_SUPPORTED)
if(NO_OMIT_FRAME_POINTER_FLAG_SUPPORTED)
    set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS_RELEASE} -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE} -fno-omit-frame-pointer")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer")
else()
    set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS_RELEASE}")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

# uncomment for memory debugging
# set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined")
# set (CMAKE_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_LINKER_FLAGS_RELWITHDEBINFO} -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined")
# set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined")
# set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined")

# if zlib is found, then we can enable reading zlib-compressed input
if (ZLIB)
  find_package(ZLIB 1.2.3)
endif()

include(CPack)
set(CPACK_RESOURCE_FILE_LICENSE "${HIGHS_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION_MAJOR "${HIGHS_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${HIGHS_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${HIGHS_VERSION_PATCH}")
set(CPACK_PACKAGE_VENDOR "University of Edinburgh")

find_program(GIT git)
if((GIT) AND (EXISTS ${HIGHS_SOURCE_DIR}/.git))
    execute_process(
        COMMAND ${GIT} describe --always --dirty
        WORKING_DIRECTORY ${HIGHS_SOURCE_DIR}
        OUTPUT_VARIABLE GITHASH OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(REGEX REPLACE "^.*-g" "" GITHASH ${GITHASH})
else()
    set(GITHASH "n/a")
endif()
message(STATUS "Git hash: " ${GITHASH})

string(TIMESTAMP TODAY "%Y-%m-%d")
message(STATUS "Compilation date: " ${TODAY})

option(ALL_TESTS "Enable full ctest set" ON)
# todo: make dependent on build_testing

include(CTest)


# Targets
include(cpp-highs)
include(c-highs)

add_subdirectory(app)
add_subdirectory(examples)

if (ALL_TESTS)
    enable_testing()
    # set BUILD_TESTING ON 
elseif (BUILD_TESTING)
    enable_testing()
endif()

if(BUILD_TESTING AND ALL_TESTS)
    enable_testing()
    add_subdirectory(check)
endif()

option(EXP "Experimental mode: run unit tests with doctest." OFF)

include(coverage)

include(CMakeDependentOption)

option(BUILD_EXAMPLES "Build examples" ON)
message(STATUS "Build examples: ${BUILD_EXAMPLES}")

CMAKE_DEPENDENT_OPTION(BUILD_CXX_EXAMPLE "Build cxx example" ON "BUILD_EXAMPLES;BUILD_CXX" OFF)
message(STATUS "Build C++: ${BUILD_CXX_EX}")
CMAKE_DEPENDENT_OPTION(BUILD_PYTHON_EXAMPLE "Build Python example" ON "BUILD_EXAMPLES;PYTHON" OFF)
message(STATUS "Build Python: ${BUILD_PYTHON_EXAMPLE}")
CMAKE_DEPENDENT_OPTION(BUILD_CSHARP_EXAMPLE "Build CSharp example" ON "BUILD_EXAMPLES;CSHARP" OFF)
message(STATUS "Build CSharp: ${BUILD_CSHARP_EXAMPLE}")

# IF BUILD_DEPS=ON THEN Force all BUILD_*=ON
CMAKE_DEPENDENT_OPTION(BUILD_ZLIB "Build the ZLIB dependency Library" OFF
  "NOT BUILD_DEPS" ON)
message(STATUS "Build ZLIB: ${BUILD_ZLIB}")

if(PYTHON)
  CMAKE_DEPENDENT_OPTION(BUILD_pybind11 "Build the pybind11 dependency Library" OFF
    "NOT BUILD_DEPS" ON)
  message(STATUS "Python: Build pybind11: ${BUILD_pybind11}")

  CMAKE_DEPENDENT_OPTION(BUILD_VENV "Create Python venv in BINARY_DIR/python/venv" OFF
    "NOT BUILD_TESTING" ON)
  message(STATUS "Python: Create venv: ${BUILD_VENV}")

  option(VENV_USE_SYSTEM_SITE_PACKAGES "Python venv can use system site packages" OFF)
  message(STATUS "Python: Allow venv to use system site packages: ${VENV_USE_SYSTEM_SITE_PACKAGES}")

  option(FETCH_PYTHON_DEPS "Install python required modules if not available" ${BUILD_DEPS})
  message(STATUS "Python: Fetch dependencies: ${FETCH_PYTHON_DEPS}")
endif()

option(JULIA "Build library and executable for Julia" OFF)

if (PYTHON)
  include(python-highs)
endif()

if (EXP)
    add_executable(doctest)
    #target_sources(doctest PRIVATE check/doctest/TestPresolveColumnSingletons.cpp)
    target_sources(doctest PRIVATE check/doctest/TestPresolveIssue.cpp)

    if (NOT APPLE)
        # triggers hanging on macOS
        target_sources(doctest PRIVATE check/doctest/TestGas11.cpp)
    endif()

    target_include_directories(doctest PRIVATE extern)
    target_link_libraries(doctest highs)
endif()
