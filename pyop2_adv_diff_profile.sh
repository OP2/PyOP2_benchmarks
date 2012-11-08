#!/bin/bash
for i in 10 20 40 80 160 320; do
  python -m cProfile -o pyop2_adv_diff.${i}.cprofile pyop2_adv_diff.py -m meshes/mesh_${i}
  gprof2dot -f pstats -n 1 pyop2_adv_diff.${i}.cprofile | dot -Tpdf -o pyop2_adv_diff.${i}.cprofile.pdf
done
