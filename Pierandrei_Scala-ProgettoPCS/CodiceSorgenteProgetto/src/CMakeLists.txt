list(APPEND poliedro_sources_src
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils.cpp
)

list(APPEND poliedro_headers_src
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Poliedro.hpp
)

list(APPEND poliedro_includes_src
    ${CMAKE_CURRENT_SOURCE_DIR}
)

set(poliedro_sources_src ${poliedro_sources_src} PARENT_SCOPE)
set(poliedro_headers_src ${poliedro_headers_src} PARENT_SCOPE)
set(poliedro_includes_src ${poliedro_includes_src} PARENT_SCOPE)
