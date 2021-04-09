add_executable (dmk_p1p0_from_file${EXECUTABLE_EXTENSION} srcs/main.f90)
target_link_libraries( dmk_p1p0_from_file${EXECUTABLE_EXTENSION} dmk linalg agmg arpack)
