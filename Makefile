CC=nvcc
CFLAGS=
OUTPUT=fdtd
SOURCES=fdtd.cu utils.o

fdtd: utils.o
	${CC} -o ${OUTPUT} ${CFLAGS} ${SOURCES}

utils.o:
	${CC} ${CFLAGS} -c utils.c
    
clean:
	rm ${OUTPUT} *.mod *.o *.s run*.sh.o* run*.sh.e* output/* 2> /dev/null

run: fdtd
	pjsub run.sh

prof: fdtd
	pjsub run_nvprof.sh
