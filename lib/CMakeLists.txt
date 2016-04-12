include_directories("${CMAKE_CURRENT_BINARY_DIR}/../include"
	"${CMAKE_CURRENT_BINARY_DIR}/../external/libprange/include"
)
file(GLOB libFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

add_library(divsufsort ${libFiles})

#target_include_directories (libprange PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

