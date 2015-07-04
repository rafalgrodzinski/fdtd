CC=gfortran-mp-4.9
CFLAGS=-g -save-temps=obj
OUTPUT=fdtd_serial
SOURCES=fdtd_serial.f90

all:
	${CC} -o ${OUTPUT} ${CFLAGS} ${SOURCES}
    
clean:
	rm ${OUTPUT} *.mod *.o *.s
