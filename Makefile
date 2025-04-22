SHELL = sh

# ATTENTION: path to SUNDIALS installation directory
prefix     = /data/project/general/ans/COD/SUNDIALS/INSTDIR
includedir = ${prefix}/fortran
libdir     = ${prefix}/lib64

# ATTENTION: path to FORTRAN compiler
F90      = /opt/psi/Programming/gcc/13.1.0/bin/gfortran
F90FLAGS = -O3 
F90LIBS  = -lm

# ------------------------------------------------------------------------------

INCLUDES  = -I${includedir}
LIBRARIES = -lsundials_fida_mod -lsundials_ida -lsundials_fcore_mod ${F90LIBS}
LINKFLAGS = -Wl,-rpath,$(libdir)

# ------------------------------------------------------------------------------

FILES =  main
FILES_DEPENDENCIES = globals FRED

OBJECTS = ${FILES:=.o}
OBJECTS_DEPENDENCIES = ${FILES_DEPENDENCIES:=.o}

EXE = FRED1.4.0.x

# ------------------------------------------------------------------------------

.SUFFIXES : .o .f90

.f90.o :
	${F90} ${F90FLAGS} ${INCLUDES} -c $<

# ------------------------------------------------------------------------------

all: ${OBJECTS_DEPENDENCIES} ${OBJECTS}
	@for i in ${FILES} ; do \
	  echo "${F90} -o ${EXE} $${i}.o ${OBJECTS_DEPENDENCIES} ${F90FLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS}" ; \
	  ${F90} -o ${EXE} $${i}.o ${OBJECTS_DEPENDENCIES} ${F90FLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS} ; \
	done
