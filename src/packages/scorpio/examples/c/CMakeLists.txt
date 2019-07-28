add_executable(simplex simplex.c)
target_link_libraries(simplex scorpio)

add_executable(attrex attrex.c)
target_link_libraries(attrex scorpio)

add_executable(strex strex.c)
target_link_libraries(strex scorpio)
