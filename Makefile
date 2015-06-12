CC=gfortran-mp-4.9
CFLAGS=
OUTPUT=fdtd_serial
SOURCES=fdtd_serial.f90

all:
	${CC} ${CFLAGS} -o ${OUTPUT} ${SOURCES}
    
clean:
	rm ${OUTPUT} *.mod *.o