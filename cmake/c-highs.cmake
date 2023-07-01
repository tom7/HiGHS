# add_c_example()
function(add_c_example FILE_NAME)
  message(STATUS "Configuring test ${FILE_NAME}: ...")
  get_filename_component(TEST_NAME ${FILE_NAME} NAME_WE)
  add_executable(${TEST_NAME} ${FILE_NAME})
  
  target_include_directories(${TEST_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  target_include_directories(${TEST_NAME} PUBLIC ${PROJECT_BINARY_DIR})

  target_link_libraries(${TEST_NAME} PRIVATE ${PROJECT_NAMESPACE}::highs)
  
  message(STATUS "Configuring test ${FILE_NAME}: ...DONE")
endfunction()