all:
	python pyop2_bench.py
clean:
	rm -rf *.flml *.pyc *.edge *.ele *.poly *.geo *.msh *.node *.stat *.vtu *.pvtu input/decomp* input/*checkpoint* input/*.stat
