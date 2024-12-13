install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Taylor-Green-Vortex)

# If testing active add test for isotropic case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/TGV")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY opencfd2.in DESTINATION ${case_dir})
endif()
