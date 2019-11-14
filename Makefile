FC=           mpifort
FFLAGS=       -O3
EXE_PAR=      percolate_par
EXE_SER=      percolate_ser
FSRC=         src/uni.f90 src/cli_info.f90 src/io.f90
OBJ=          $(FSRC:src/%.f90=%.o)

.PHONY: all
all: $(EXE_PAR) $(EXE_SER) clean

.PHONY: clean
clean:
	@rm *.mod *.o
	@echo Done with cleanup.

$(EXE_PAR): $(OBJ)
	@$(FC) $^ src/$(EXE_PAR).f90 -o $(EXE_PAR) $(FFLAGS)
	@echo Done making $(EXE_PAR).

$(EXE_SER): $(OBJ)
	@$(FC) $^ src/$(EXE_SER).f90 -o $(EXE_SER) $(FFLAGS)
	@echo Done making $(EXE_SER).

$(OBJ): $(FSRC)
	@$(FC) -c $^ $(FFLAGS)
