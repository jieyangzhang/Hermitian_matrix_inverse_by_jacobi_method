cmake_minimum_required(VERSION 3.8)
project( jacobi_inverse )

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
LINK_LIBRARIES(m)

add_executable(jacobi_inverse_test test/jacobi_inverse_test.c src/jacobi_inverse.c)
target_include_directories(jacobi_inverse_test PUBLIC ${pub_includes})
target_compile_options(jacobi_inverse_test PUBLIC ${DEFINES_OPTIONS})
target_link_libraries(jacobi_inverse_test PUBLIC ${jacobi_inverse_test_DEPS} ${ld_options})

INSTALL(TARGETS jacobi_inverse_test RUNTIME DESTINATION bin)