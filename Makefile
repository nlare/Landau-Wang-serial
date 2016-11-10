all:
	g++ -fopenmp -g landau-wang.cpp -o landau-wang -lboost_system -lboost_filesystem -xhost -O3
clean:
	rm landau-wang
