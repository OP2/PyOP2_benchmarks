#!/bin/bash
flml=ufl_advection_diffusion
for backend in sequential cuda opencl_amd opencl_intel opencl_nvidia; do
  case $backend in
    opencl_amd)
      ctx=0
      b=opencl
      ;;
    opencl_intel)
      ctx=1
      b=opencl
      ;;
    opencl_nvidia)
      ctx=2:0
      b=opencl
      ;;
    *)
      b=$backend
      ;;
  esac
  for i in 10 20 40 80 160 320; do
    echo Running ${flml} with ${backend} backend and mesh size ${i}
    python -i <<EOF
import os
from parameters import *
with open('${flml}.profile.flml.template') as f1, \
        open(os.path.join('flmls','${flml}.${i}.profile.flml'), 'w') as f2:
    f2.write(f1.read() % {
        'mesh': 'meshes/mesh_${i}',
        'diffusivity': diffusivity,
        'current_time': current_time,
        'dt': dt,
        'endtime': (endtime - dt),
        'backend': '${b}',
        'profile': '${i}.${backend}'
        })
EOF
    PYOPENCL_CTX=${ctx} ${FLUIDITY_DIR}/bin/fluidity flmls/${flml}.${i}.profile.flml
    gprof2dot -f pstats -n 1 ${flml}.${i}.${backend}.cprofile.dat \
      | dot -Tpdf -o ${flml}.${i}.${backend}.cprofile.pdf
    gprof ${FLUIDITY_DIR}/bin/fluidity gmon.out > ${flml}.${i}.${backend}.gprof.dat
    gprof2dot -n 1 ${flml}.${i}.${backend}.gprof.dat \
      | dot -Tpdf -o ${flml}.${i}.${backend}.gprof.pdf
    rm gmon.out
  done
done
