CC=pgf90
CFLAGS=-g -save-temps=obj -Mcuda
OUTPUT=fdtd
SOURCES=fdtd.f90

all:
	${CC} -o ${OUTPUT} ${CFLAGS} ${SOURCES}
    
clean:
	rm ${OUTPUT} *.mod *.o *.s

run:
	pjsub run.sh