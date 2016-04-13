configure_file("${CMAKE_CURRENT_SOURCE_DIR}/divsufsort.h" "${CMAKE_CURRENT_BINARY_DIR}/divsufsort.h" COPYONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/divsufsort.h" DESTINATION include)
