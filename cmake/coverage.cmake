
# Coverage part
# 'make coverage' to start the coverage process
option(HIGHS_COVERAGE "Activate the code coverage compilation" OFF)
if (HIGHS_COVERAGE)
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_C_FLAGS_DEBUG    "${CMAKE_C_FLAGS_DEBUG}    -O0 --coverage")
    set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG}  -O0 --coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -O0 --coverage")
  endif ()
  if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -fprofile-arcs -ftest-coverage -Xclang -coverage-cfg-checksum -Xclang -coverage-no-function-names-in-data -Xclang -coverage-version='408*'")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage -Xclang -coverage-cfg-checksum -Xclang -coverage-no-function-names-in-data -Xclang -coverage-version='408*'")
  endif ()
endif ()

if (HIGHS_COVERAGE)
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(FATAL_ERROR "Warning: to enable coverage, you must compile in DEBUG mode")
  endif ()
endif ()

if (HIGHS_COVERAGE)
  if (WIN32)
    message(FATAL_ERROR "Error: code coverage analysis is only available under Linux for now.")
  endif ()

  find_program(GCOV_PATH gcov)
  find_program(LCOV_PATH lcov)
  find_program(GENHTML_PATH genhtml)

  if (NOT GCOV_PATH)
    message(FATAL_ERROR "gcov not found! Please install lcov and gcov. Aborting...")
  endif ()

  if (NOT LCOV_PATH)
    message(FATAL_ERROR "lcov not found! Please install lcov and gcov. Aborting...")
  endif ()

  if (NOT GENHTML_PATH)
    message(FATAL_ERROR "genhtml not found! Please install lcov and gcov. Aborting...")
  endif ()

  # Capturing lcov counters and generating report
  if (NOT CI)
    add_custom_target(coverage
                    COMMAND ${LCOV_PATH} --directory ${CMAKE_BINARY_DIR} --zerocounters
                    COMMAND ${LCOV_PATH} --capture --initial --directory ${CMAKE_BINARY_DIR}/bin --output-file ${CMAKE_BINARY_DIR}/coverage.info
                    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR} ${CMAKE_CTEST_COMMAND} -LE "(LONG|FAIL)" || true
                    COMMAND ${LCOV_PATH} --capture --directory ${CMAKE_BINARY_DIR}/bin --directory ${CMAKE_BINARY_DIR}/highs --directory ${CMAKE_BINARY_DIR}/app --directory ${CMAKE_BINARY_DIR}/check --output-file ${CMAKE_BINARY_DIR}/coverage.info
                    COMMAND ${LCOV_PATH} --remove "*/usr/include/*" --output-file ${CMAKE_BINARY_DIR}/coverage.info.cleaned
                    COMMAND ${GENHTML_PATH} -o ${CMAKE_BINARY_DIR}/coverage ${CMAKE_BINARY_DIR}/coverage.info.cleaned
                    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/coverage.info.cleaned
            VERBATIM
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                    COMMENT "Resetting code coverage counters to zero.
    Processing code coverage counters and generating report.
    You can zip the directory ${CMAKE_BINARY_DIR}/coverage and upload the content to a web server.")
  else()
    add_custom_target(ci_cov
                    COMMAND ${LCOV_PATH} --directory ${CMAKE_BINARY_DIR} --zerocounters
                    COMMAND ${LCOV_PATH} --capture --initial --directory ${CMAKE_BINARY_DIR}/bin --output-file ${CMAKE_BINARY_DIR}/coverage.info
                    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR} ${CMAKE_CTEST_COMMAND} -LE "(LONG|FAIL)" || true
                    COMMAND ${LCOV_PATH} --capture --directory ${CMAKE_BINARY_DIR}/bin --directory ${CMAKE_BINARY_DIR}/highs --directory ${CMAKE_BINARY_DIR}/app --directory ${CMAKE_BINARY_DIR}/check --output-file ${CMAKE_BINARY_DIR}/coverage.info
		    VERBATIM
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                    COMMENT "Resetting code coverage counters to zero.")
  endif()
endif ()