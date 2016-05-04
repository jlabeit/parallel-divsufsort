add_executable(divsufsortDemo main.cpp)
target_link_libraries (divsufsortDemo cilkrts divsufsort libprange)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../include")
