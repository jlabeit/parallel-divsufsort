## Add libprange
set(libprange_cmake_file ${CMAKE_CURRENT_LIST_DIR}/libprange/CMakeLists.txt)
if(NOT EXISTS ${libprange_cmake_file})
    EXECUTE_PROCESS(COMMAND git submodule update --init external/libprange
                    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/..
                    OUTPUT_QUIET
            )
endif(NOT EXISTS ${libprange_cmake_file})
add_subdirectory(libprange)
