f77=mpif90
#f77 = ifort
opt= -O3   
srs= OCFD_parameters.f90 opencfd.f90 OCFD_readpara.f90 OCFD_BC_Def.f90  \
OCFD_BC.f90 OCFD_dxyz_viscous.f90 OCFD_flux_split.f90 OCFD_handle_NegativeT.f90 \
OCFD_init.f90 OCFD_inviscous.f90 OCFD_inviscous_character.f90  \
OCFD_IO.f90 OCFD_IO_mpi.f90 OCFD_mesh.f90 OCFD_mpi.f90  \
OCFD_NS_Solver.f90 OCFD_User_Def.f90  \
OCFD_Schemes.f90  OCFD_Schemes_1.f90 \
OCFD_Time_Adv.f90 OCFD_viscous.f90 OCFD_ANA.f90 OCFD_BC_Boundarylayers.f90 \
OCFD_ANA_savedata.f90 OCFD_filtering.f90 OCFD_SetHybrid.f90 OCFD_GhostCell.f90 OCFD_GhostCell_UserDef.f90
	 
OBJS=$(srs:.f90=.o)

%.o:%.f90
	$(f77) $(opt) -c $<

default: $(OBJS)
	$(f77) -O3  -o opencfd-2.2b.out $(OBJS)
test0: $(OBJS)
	$(f77) -O3  -o test-dxyz.out OCFD_parameters.f90 opencfd.f90 OCFD_readpara.f90 OCFD_mpi.f90 OCFD_dxyz_viscous.f90 test-dxyz.f90
test1: $(OBJS)
	$(f77) -O3  -o test-dx-invis.out test-dx-invis.f90 OCFD_parameters.o opencfd.f90 OCFD_readpara.o OCFD_mpi.o  OCFD_inviscous.o OCFD_Schemes.o  OCFD_Schemes_1.o \
	OCFD_flux_split.o OCFD_User_Def.o OCFD_SetHybrid.o  OCFD_dxyz_viscous.o
clean:
	rm -f *.out *.o *.mod
