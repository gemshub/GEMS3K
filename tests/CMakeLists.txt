add_executable(pstall_default_test pstall_default_test.cpp)

if(TARGET GEMS3K-static)
    target_link_libraries(pstall_default_test PRIVATE GEMS3K-static)
else()
    target_link_libraries(pstall_default_test PRIVATE GEMS3K)
endif()

target_compile_definitions(
    pstall_default_test
    PRIVATE GEMS3K_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/nodearray-gem-build/TestRW"
)

add_test(NAME pstall_default_reset COMMAND pstall_default_test)
