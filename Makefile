FC=           ifort
FC_PAR=       mpif08

FFLAGS=       -O3

EXE_PAR=      percolate_par
EXE_SER=      percolate_ser

FSRC_COMMON=  src/uni.f90 src/cli_info.f90 src/io.f90
FSRC_PAR=			src/cart_comm.f90

OBJ_COMMON=   $(FSRC_COMMON:src/%.f90=%.o)
OBJ_PAR=      $(FSRC_PAR:src/%.f90=%.o)


.PHONY: all
all: $(EXE_PAR) $(EXE_SER) clean

.PHONY: clean
clean:
	@rm *.mod *.o
	@echo Done with cleanup.


$(EXE_PAR): $(OBJ_COMMON) $(OBJ_PAR)
	@$(FC_PAR) $^ src/$(EXE_PAR).f90 -o $(EXE_PAR) $(FFLAGS)
	@echo Done making $(EXE_PAR).

$(EXE_SER): $(OBJ_COMMON)
	@$(FC) $^ src/$(EXE_SER).f90 -o $(EXE_SER) $(FFLAGS)
	@echo Done making $(EXE_SER).


$(OBJ_COMMON): $(FSRC_COMMON)
	@$(FC) -c $^ $(FFLAGS)

$(OBJ_PAR): $(FSRC_PAR)
	@$(FC_PAR) -c $^ $(FFLAGS)
