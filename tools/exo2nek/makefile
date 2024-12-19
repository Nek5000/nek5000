prefix = $(bin_nek_tools)
INC = -I./3rd_party/seacas-exodus/packages/seacas/libraries/exodus_for/include
LIBS += ./3rd_party/seacas-exodus/build/packages/seacas/libraries/exoIIv2for32/libexoIIv2for32.a
LIBS += ./3rd_party/seacas-exodus/build/packages/seacas/libraries/exodus/libexodus.a
LIBS += ./3rd_party/netcdf/install/lib/libnetcdf.a

OBJS = mod_SIZE.o exo2nek.o tet2hex.o periodicity.o quadratic.o non_right_hand_check.o  byte.o speclib.o mxm.o mod_OCTUPLESIZE.o octuple.o shell.o

all: lib exo2nek

exo2nek: $(OBJS)
	$(FC) $(FFLAGS) -o $(prefix)/exo2nek $^ $(LIBS) $(LDFLAGS)

clean:
	@rm -f *.o *.mod 
	@cd ./3rd_party ; rm -rf seacas-exodus netcdf

lib:
	@cd 3rd_party; env CC="$(CC)" CFLAGS="$(CFLAGS)" FC="$(FC)" FFLAGS="$(FFLAGS)" ./install

mod_SIZE.o      : mod_SIZE.f90			;  $(FC) -c $(FFLAGS) mod_SIZE.f90
exo2nek.o	: exo2nek.f90			;  $(FC) -c $(FFLAGS) $(INC) exo2nek.f90
tet2hex.o	: tet2hex.f90			;  $(FC) -c $(FFLAGS) $(INC) tet2hex.f90
periodicity.o	: periodicity.f90   ;  $(FC) -c $(FFLAGS) $(INC) periodicity.f90
quadratic.o     : quadratic.f90     ;  $(FC) -c $(FFLAGS) $(INC) quadratic.f90
non_right_hand_check.o     : non_right_hand_check.f90     ;  $(FC) -c $(FFLAGS) $(INC) non_right_hand_check.f90
byte.o		: ../../core/byte.c		;  $(CC) -c $(CFLAGS) ../../core/byte.c
speclib.o	: ../../core/speclib.f		;  $(FC) -c $(FFLAGS) ../../core/speclib.f
mxm.o		: mxm.f				;  $(FC) -c $(FFLAGS) mxm.f
mod_OCTUPLESIZE.o      : mod_OCTUPLESIZE.f90			;  $(FC) -c $(FFLAGS) mod_OCTUPLESIZE.f90
octuple.o	: octuple.f90			;  $(FC) -c $(FFLAGS) $(INC) octuple.f90
shell.o	: shell.f90			;  $(FC) -c $(FFLAGS) $(INC) shell.f90