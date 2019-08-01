# ============================================================================
# module swap PrgEnv-cray PrgEnv-gnu
# module load zlib cray-tpsl cray-trilinos cray-hdf5-parallel fox
# ============================================================================
export CRAYPE_LINK_TYPE=dynamic

# Link FOX
FOX=/ivec/cle52/magnus/apps/PrgEnv-gnu/5.2.25/fox/4.1.2
FOXMOD=$(FOX)/finclude
FOXLIB=$(FOX)/lib
FOXLIBS= -lFoX_wkml -lFoX_dom -lFoX_sax -lFoX_wcml -lFoX_wxml -lFoX_common -lFoX_utils -lFoX_fsys
FOXFLAGS= -I$(FOXMOD) 
LDFOXFLAGS= -L$(FOXLIB) $(FOXLIBS)

# ============================================================================

# Link HDF5
H5LDFLAGS = -dynamic 
H5LIBS = -lhdf5_fortran_parallel -lhdf5_parallel  -lhdf5hl_fortran_parallel -lhdf5_hl_parallel -lz

# ============================================================================

# Link Zoltan
ZOLTANLDFLAGS = -dynamic  
ZOLTANFLAGS =
ZOLTANLIBS = -lzoltan 

# ============================================================================

# Link METIS
METIS=/group/partner985/SPM/metis-gnu
METISLDFLAGS = -dynamic -L${METIS}/lib 
METISLIBS =  -lmetis

# ============================================================================

# Executable name
EXE=badlands-gnu

# C compiler
BADLANDS_C = cc 
BADLANDS_F = ftn

# Rules to make library
AR = ar -rcs

# Fortran optimisation flags
#FCFLAGS = -O3 -cpp -fpic
FCFLAGS = -O3 -cpp -fpic -ffree-line-length-none
#FCFLAGS= -O0 -g  -Wall -fbacktrace -lstdc++ -cpp -fcheck=bounds -finit-real=nan\
	-ffpe-trap=zero,overflow,invalid -ffree-form -fno-common\
	-Wtabs -Wunused-parameter -Wuninitialized  -ffree-line-length-none \
	-fdump-fortran-optimized -fdump-tree-original 

# C optimisation flags
CFLAGS = -O3 -fpic 







