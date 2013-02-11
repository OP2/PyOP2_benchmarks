#!/bin/bash
flml=ufl_advection_diffusion
mesh=meshes/square.0.

mkdir -p profiling

for backend in sequential; do
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
  for meshsize in 000002; do
    echo Running ${flml} with ${backend} backend and mesh size ${meshsize}
    python -i <<EOF
import os
from parameters import *
with open('${flml}.profile.flml.template') as f1, \
        open(os.path.join('flmls','${flml}.${meshsize}.profile.flml'), 'w') as f2:
    f2.write(f1.read() % {
        'mesh': '${mesh}${meshsize}',
        'diffusivity': diffusivity,
        'current_time': current_time,
        'dt': dt,
        'endtime': (endtime - dt),
        'backend': '${b}'
        })
EOF
    PYOPENCL_CTX=${ctx} time mpiexec ${FLUIDITY_DIR}/bin/fluidity flmls/${flml}.${meshsize}.profile.flml
    for j in `seq 0 11`; do
      datfile=profiling/${flml}.${meshsize}.${backend}.${j}.cprofile.dat
      python concat.py "*.$j.cprofile.part" ${datfile}
      #gprof2dot -f pstats -n 1 ${datfile} | dot -Tpdf -o ${datfile}.pdf
      rm *.$j.cprofile.part
      #gprof ${FLUIDITY_DIR}/bin/fluidity gmon.out > ${flml}.${meshsize}.${backend}.gprof.dat
      #gprof2dot -n 1 ${flml}.${meshsize}.${backend}.gprof.dat \
        #| dot -Tpdf -o ${flml}.${meshsize}.${backend}.gprof.pdf
      #rm *.part gmon.out
    done
  done
done
