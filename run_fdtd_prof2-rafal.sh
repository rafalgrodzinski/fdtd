#!/bin/sh

#PJM -L rscunit=gwacsg
#PJM -L rscgrp=gpu
#PJM -L vnode-core=6
#PJM -x gpu_per_vnode=1

/gwfefs/opt/x86_64/cuda/7.5/bin/nvprof --metrics dram_read_throughput,dram_write_throughput,ldst_fu_utilization,cf_fu_utilization,alu_fu_utilization,tex_fu_utilization,achieved_occupancy ./fdtd
