smallpt_transient: smallpt_transient.cpp
	g++ -O3 -fopenmp smallpt_transient.cpp -o smallpt_transient

.PHONY: render
samples=100
render: smallpt_transient
	./smallpt_transient $(samples)
	convert -delay 8 -loop 0 [0-9]*.ppm image.gif
	rm [0-9]*.ppm
