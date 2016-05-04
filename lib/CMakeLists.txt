include_directories("${CMAKE_CURRENT_BINARY_DIR}/../include"
	"${CMAKE_CURRENT_BINARY_DIR}/../external/libprange/include"
	"${CMAKE_CURRENT_BINARY_DIR}"
)
configure_file(
	"${CMAKE_CURRENT_SOURCE_DIR}/parallel.hpp.cmake"
	"${CMAKE_CURRENT_BINARY_DIR}/parallel.hpp"
	)
file(GLOB libFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
add_library(divsufsort ${libFiles})
install(TARGETS divsufsort
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)
