CC=nvcc
CFLAGS=
OUTPUT=fdtd
SOURCES=fdtd.cu utils.cu fdtd_calculations.cu

fdtd: ${SOURCES}
	${CC} -o ${OUTPUT} ${CFLAGS} ${SOURCES}
    
clean:
	rm -f ${OUTPUT} *.mod *.o *.s run*.sh.o* run*.sh.e* output/* fdtd_prof

run: fdtd
	pjsub run.sh

prof: fdtd
	pjsub run_prof.sh
