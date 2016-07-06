#!/bin/sh

#PJM -L rscunit=gwacsg
#PJM -L rscgrp=gpu
#PJM -L vnode=1
#PJM -L vnode-core=8
#PJM -x gpu_per_vnode=1
#PJM -j
#PJM -L elapse=2:00:00

time /gwfefs/opt/x86_64/cuda/7.5/bin/nvprof --metrics achieved_occupancy,gld_efficiency,gld_throughput,gst_efficiency,gst_throughput,tex_cache_throughput ./fdtd
