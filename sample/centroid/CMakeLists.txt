add_executable(dymp_centroid
  main.cpp
  )

target_link_libraries(dymp_centroid
  PRIVATE dymp
  PRIVATE Eigen3::Eigen
  PRIVATE glfw
)
