list(APPEND poliedro_sources_test
    ${CMAKE_CURRENT_SOURCE_DIR}/test.hpp
)

list(APPEND poliedro_includes_test
    ${CMAKE_CURRENT_SOURCE_DIR}
)

set(poliedro_sources_test ${poliedro_sources_test} PARENT_SCOPE)
set(poliedro_headers_test ${poliedro_headers_test} PARENT_SCOPE)
set(poliedro_includes_test ${poliedro_includes_test} PARENT_SCOPE)