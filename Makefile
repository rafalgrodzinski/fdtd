CC=nvcc
CFLAGS=
OUTPUT=fdtd
SOURCES=fdtd.cu utils.cu fdtd_calculations.cu

fdtd: ${SOURCES}
	${CC} -o ${OUTPUT} ${CFLAGS} ${SOURCES}
    
clean:
	rm ${OUTPUT} *.mod *.o *.s run*.sh.o* run*.sh.e* output/* 2> /dev/null

run: fdtd
	pjsub run.sh

prof: fdtd
	pjsub run_nvprof.sh
