#  Copyright (C) 2002 Regents of the University of Michigan,
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

include ../../../Makefile.def

SEARCHDIR = -I${SHAREDIR} -I${MAGNETOGRAMDIR}

include ../../../Makefile.conf

MY_LIB = ${LIBDIR}/libEMPIRICALSC.a

-include Makefile.DEPEND
-include Makefile.RULES

OBJECTS = \
	ModWsa.o

OBJECTS_EXE= stand_alone_wsa.o

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${OBJECTS} ${OBJECTS_EXE}

LIBSHARE:
	cd ${SHAREDIR}; make LIB

LIBMAGNETOGRAM: LIBSHARE
	cd ${MAGNETOGRAMDIR}; make LIB

LIB:	DEPEND LIBMAGNETOGRAM
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}: ${OBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${OBJECTS}

################## WSA   ###################

WSA: 	LIB
	make ${WSAEXE}
	@echo ${WSAEXE} has been compiled.
	ln -sf ${WSAEXE} .

WSAEXE = ${BINDIR}/WSA.exe

OBJECTS_WSA= \
	stand_alone_wsa.o

stand_alone_wsa.o: ${OBJECTS}

${WSAEXE}: ${OBJECTS_WSA}
	${LINK.f90} -o ${WSAEXE} ${OBJECTS_WSA} ${OBJECTS} \
		-L${LIBDIR} -lMAGNETOGRAM -lSHARE ${Lflag}

#######################################################################
# Setting BLESS=YES or Y will copy the solution into the reference solution
BLESS=NO

DIFFNUM = ${SCRIPTDIR}/DiffNum.pl -BLESS=${BLESS}



clean: cleanfiles

distclean: clean

