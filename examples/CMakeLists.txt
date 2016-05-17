add_executable(timing timing.cpp)
target_link_libraries (timing cilkrts divsufsort libprange)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../include")
