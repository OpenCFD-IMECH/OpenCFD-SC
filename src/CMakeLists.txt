
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

add_executable(opencfd 
                    OCFD_parameters.f90
                    opencfd.f90
                    OCFD_readpara.f90
                    OCFD_BC_Def.f90
                    OCFD_BC.f90
                    OCFD_dxyz_viscous.f90
                    OCFD_flux_split.f90
                    OCFD_flowinit.f90
                    OCFD_handle_NegativeT.f90
                    OCFD_init.f90
                    OCFD_inviscous.f90
                    OCFD_inviscous_character.f90
                    OCFD_IO.f90
                    OCFD_IO_mpi.f90
                    OCFD_mesh.f90
                    OCFD_mpi.f90
                    OCFD_NS_Solver.f90
                    OCFD_User_Def.f90
                    OCFD_Schemes.f90
                    OCFD_Schemes_1.f90
                    OCFD_Time_Adv.f90
                    OCFD_viscous.f90
                    OCFD_ANA.f90
                    OCFD_BC_Boundarylayers.f90
                    OCFD_ANA_savedata.f90
                    OCFD_filtering.f90
                    OCFD_SetHybrid.f90
                    OCFD_GhostCell.f90
                    OCFD_GhostCell_UserDef.f90)

target_link_libraries(opencfd)
if (MPI_FOUND)
  target_link_libraries(opencfd PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

install(TARGETS opencfd
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
