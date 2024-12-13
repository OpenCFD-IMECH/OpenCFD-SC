install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/blunt-cylinder)

# If testing active add test for blunt-cylinder case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/blunt-cylinder")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY opencfd2.in DESTINATION ${case_dir})
  file(COPY opencfd.dat DESTINATION ${case_dir})
  file(COPY OCFD-grid.dat DESTINATION ${case_dir})
endif()
