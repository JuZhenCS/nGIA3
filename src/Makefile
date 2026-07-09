all:
	nvcc -O3 -Xcompiler -fopenmp -arch=sm_120 main.cu -o ngia2
run:
	./ngia2 clustering -p ../protein.packed -r result -i 67
test:
	g++ test_makedb.hpp -lgtest -lgtest_main
	./a.out

