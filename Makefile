CC=pgfortran
CFLAGS=-Mcuda -Minfo
OUTPUT=fdtd
SOURCES=fdtd.f90

all:
	${CC} -o ${OUTPUT} ${CFLAGS} ${SOURCES}
    
clean:
	rm ${OUTPUT} *.mod *.o *.s run*.sh.o* run*.sh.e* 2> /dev/null

run:
	pjsub run.sh

runc:
	pjsub run_cuda.sh

interact:
	pjsub --interact -L rscunit=gwacsg
