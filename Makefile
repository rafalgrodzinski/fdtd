CC=nvcc
CFLAGS=
OUTPUT=fdtd
SOURCES=fdtd.cu utils.cu fdtd_calculations.cu

fdtd: ${SOURCES}
	${CC} -o ${OUTPUT} ${CFLAGS} ${SOURCES}
    
clean:
	rm -rf ${OUTPUT} *.mod *.o *.s run*.sh.o* run*.sh.e* output fdtd_prof

run: fdtd
	if [ ! -d output ]; then mkdir output; fi; pjsub run.sh

prof: fdtd
	if [ ! -d output ]; then mkdir output; fi; pjsub run_prof.sh
