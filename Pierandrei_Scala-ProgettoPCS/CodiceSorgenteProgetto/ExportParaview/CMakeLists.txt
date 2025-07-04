list(APPEND poliedro_sources_paraview
    ${CMAKE_CURRENT_SOURCE_DIR}/UCDUtilities.cpp
)

list(APPEND poliedro_headers_paraview
    ${CMAKE_CURRENT_SOURCE_DIR}/UCDUtilities.hpp
)

list(APPEND poliedro_includes_paraview
    ${CMAKE_CURRENT_SOURCE_DIR}
)

set(poliedro_sources_paraview ${poliedro_sources_paraview} PARENT_SCOPE)
set(poliedro_headers_paraview ${poliedro_headers_paraview} PARENT_SCOPE)
set(poliedro_includes_paraview ${poliedro_includes_paraview} PARENT_SCOPE)

