#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

include ../../../Makefile.def

SEARCHDIR = -I${SHAREDIR} -I${MAGNETOGRAMDIR}

include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libEMPIRICALEE.a

-include Makefile.DEPEND
-include Makefile.RULES

OBJECTS = \
	EEE_ModCommonVariables.o \
	EEE_ModMain.o \
	EEE_ModTD99.o \
	EEE_ModGL98.o \
	EEE_ModArch.o \
	EEE_ModGetB0.o \
	EEE_ModShearFlow.o \
	EEE_ModCms.o

OBJECTS_EXE= flux_rope_magnetogram.o
LIBSHARE  = ${LIBDIR}/libSHARE.a
LIBMAGNETOGRAM = ${LIBDIR}/libMAGNETOGRAM.a

DEPEND: ${LIBSHARE} ${LIBMAGNETOGRAM}
	@perl ${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${OBJECTS} ${OBJECTS_EXE}

LIB:
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}:  DEPEND ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}


${LIBSHARE}:
	cd ${SHAREDIR}; make LIB

${LIBMAGNETOGRAM}:
	cd ${MAGNETOGRAMDIR}; make LIB

flux_rope_magnetogram.o : ${LIBSHARE} ${LIBMAGNETOGRAM} ${MY_LIB}

FRM:	DEPEND
	cd ${SHAREDIR}; make LIB
	cd ${MAGNETOGRAMDIR}; make LIB
	make ${BINDIR}/FRMAGNETOGRAM.exe

${BINDIR}/FRMAGNETOGRAM.exe: flux_rope_magnetogram.o
	${LINK.f90} -o ${BINDIR}/FRMAGNETOGRAM.exe \
		flux_rope_magnetogram.o \
		-L${LIBDIR} -lSHARE -lEMPIRICALEE -lMAGNETOGRAM ${Lflag}

#######################################################################
# Setting BLESS=YES or Y will copy the solution into the reference solution
BLESS=NO

DIFFNUM = ${SCRIPTDIR}/DiffNum.pl -BLESS=${BLESS}

test_td22.o: test_td22.f90 ${LIBSHARE} ${MY_LIB}
	${COMPILE.f90} ${Cflag0} test_td22.f90

test_td22_exe: test_td22.o
	${LINK.f90} -o test_td22.exe test_td22.o \
		-L${LIBDIR} -lSHARE -lEMPIRICALEE ${Lflag}

test_td22_check:
	-@(${DIFFNUM} -b -r=1e-7 \
		test_fields.out test_fields.ref \
                        >> test_td22.diff)

test_td22:
	@echo test_td22 starting     >  test_td22.diff
	@make test_td22_exe	       >> test_td22.diff
	@echo ./test_td22.exe        >> test_td22.diff
	-@(rm -f test_fields.out)
	-@(./test_td22.exe > test_td22.diff)
	-@(make test_td22_check)
	@ls -l test_td22.diff
#######################################################################

clean: cleanfiles

distclean: clean

