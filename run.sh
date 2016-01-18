#!/bin/sh

#PJM -L rscunit=gwacsg
#PJM -L rscgrp=gpu
#PJM -L vnode-core=6
#PJM -x gpu_per_vnode=1

/gwfefs/opt/x86_64/cuda/7.5/bin/nvprof --metrics achieved_occupancy,gld_efficiency,gld_throughput,gst_efficiency,gst_throughput,tex_cache_throughput ./fdtd
/gwfefs/opt/x86_64/cuda/7.5/bin/nvprof ./fdtd
