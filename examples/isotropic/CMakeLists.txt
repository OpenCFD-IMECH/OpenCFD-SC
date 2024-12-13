install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/isotropic)

# If testing active add test for isotropic case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/isotropic")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY opencfd2.in DESTINATION ${case_dir})
  file(COPY opencfd.dat DESTINATION ${case_dir})
  file(COPY OCFD-grid.dat DESTINATION ${case_dir})
endif()
