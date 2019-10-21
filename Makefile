FC=           gfortran
FFLAGS=       -O3
EXE=          percolate
INSTALL_PATH= /usr/bin/
FSRC=         src/uni.f90 src/map_class.f90 src/sorted_clusters_class.f90 \
			        src/color_map_class.f90 src/cli_info.f90 src/io.f90 src/pgm_out.f90 \
							src/percolate.f90
OBJ=          $(FSRC:src/%.f90=%.o)

##all:     compile $(EXE). This is the default rule that
##         gets executed when using make.
##
.PHONY: all
all: $(EXE) clean

##install: installs $(EXE) to a PATH directory. Default is
##         /usr/bin/. Can be changed with the
##         variable INSTALL_PATH.
##
.PHONY: install
install:
	@cp $(EXE) $(INSTALL_PATH)

##clean:   remove all object and module files generated
##         during the build process.
##
.PHONY: clean
clean:
	@rm *.mod *.o
	@echo Done with cleanup.

##help:    displays this help message.
##
.PHONY: help
help: Makefile
	@sed -n 's/^##//p' $<

$(EXE): $(OBJ)
	@$(FC) $^ -o $(EXE) $(FFLAGS)
	@echo Done making $(EXE).

$(OBJ): $(FSRC)
	@$(FC) -c $^ $(FFLAGS)
