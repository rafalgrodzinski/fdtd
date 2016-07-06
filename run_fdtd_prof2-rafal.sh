#!/bin/sh

#PJM -L rscunit=gwacsg
#PJM -L rscgrp=gpu
#PJM -L vnode=1
#PJM -L vnode-core=8
#PJM -x gpu_per_vnode=1
#PJM -j
#PJM -L elapse=2:00:00

time /gwfefs/opt/x86_64/cuda/7.5/bin/nvprof --metrics dram_read_throughput,dram_write_throughput,ldst_fu_utilization,cf_fu_utilization,alu_fu_utilization,tex_fu_utilization,achieved_occupancy ./fdtd
