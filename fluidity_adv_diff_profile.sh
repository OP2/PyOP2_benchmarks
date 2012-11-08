#!/bin/bash
flml=advection_diffusion
for i in 10 20 40 80 160 320; do
  ${FLUIDITY_DIR}/bin/fluidity flmls/${flml}.${i}.flml
  gprof ${FLUIDITY_DIR}/bin/fluidity gmon.out > ${flml}.${i}.gprof.dat
  gprof2dot -n 1 ${flml}.${i}.gprof.dat \
    | dot -Tpdf -o ${flml}.${i}.gprof.pdf
  rm gmon.out
done
