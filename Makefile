###############################################################################
#		Makefile for NRG and associated package
###############################################################################

include Make.sys
 
SRCF	= tki_ga.f90 fourierstm.f90  xystm.f90 smearing_lin.f90 smearing_lin2.f90 smearing_lin3.f90  smearing_lin4.f90  arpes.f90 ef.sh atan.f90 ti1d.f90
SRCM	= tki_ga.f90

SRCF77	= ${SECOND_FILE}

OBJ     = ${SRCF:.f90=.o} ${SRCC:.c=.o} ${SRCF77:.f=.o} 
MODULES = ${SRCM:.f90=.mod}
#EXTRAFILES= Makefile  NRG.in README smearing.f90 green.f90 par.in Models.txt smearing_lin.f90 conv.f90 broad.sh G55.f90 G33.f90 G66.f90 G22.f90 G33_p32.f90 conv_fermi.f90 conv_analytic.f90 green_old.f90 conv_appelbaum.f90 par_appel.in
ALLFILES= ${SRCF} $(SRCF77) ${SRCC} ${EXTRAFILES} ${OTHERFILES} Makefile #Make.sys
 
${TARGET}: ${OBJ}              	    	# First rule is the default
	${FC} ${DATAFL} ${FFLAGS} ${OBJ} ${LIBS} -o $@
tki_ga: tki_ga.f90
	${FC} ${FFLAGS} ${LIBS} tki_ga.f90 -o tki_ga
smearing_lin: smearing_lin.f90
	${FC} ${FFLAGS} ${LIBS} smearing_lin.f90 -o smearing_lin
smearing_lin2: smearing_lin2.f90
	${FC} ${FFLAGS} ${LIBS} smearing_lin2.f90 -o smearing_lin2
smearing_lin3: smearing_lin3.f90
	${FC} ${FFLAGS} ${LIBS} smearing_lin3.f90 -o smearing_lin3
smearing_lin4: smearing_lin4.f90
	${FC} ${FFLAGS} ${LIBS} smearing_lin4.f90 -o smearing_lin4
arpes: arpes.f90
	${FC} ${FFLAGS} ${LIBS} arpes.f90 -o arpes
xystm: xystm.f90
	${FC} ${FFLAGS} ${LIBS} xystm.f90 -o xystm
atan: atan.f90
	${FC} ${FFLAGS} ${LIBS} atan.f90 -o atan
fourierstm: fourierstm.f90
	${FC} ${FFLAGS} ${LIBS} fourierstm.f90 -o fourierstm
tki_ga223: tki_ga223.f90
	${FC} ${FFLAGS} ${LIBS} tki_ga223.f90 -o tki_ga223
read_wannier: read_wannier.f90
	${FC} ${FFLAGS} ${LIBS} read_wannier.f90 -o read_wannier
read_wannier_ryan: read_wannier_ryan.f90
	${FC} ${FFLAGS} ${LIBS} read_wannier_ryan.f90 -o read_wannier_ryan
ti1d: ti1d.f90
	${FC} ${FFLAGS} ${LIBS} ti1d.f90 -o ti1d
g0: g0.f90
	${FC} ${FFLAGS} ${LIBS} g0.f90 -o g0
clean	:
	rm -f ${TARGET} ${OBJ} ${MODULES} *.mod *.d work.*
new	: clean
	make ${TARGET}
package	:	${ALLFILES}
	tar cvf - ${ALLFILES} | gzip -c > tki_ga.tar



