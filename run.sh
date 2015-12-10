#!/bin/sh

#PJM -L rscunit=gwacsg
#PJM -L rscgrp=gpu
#PJM -L vnode-core=6
#PJM -x gpu_per_vnode=1

./fdtd
