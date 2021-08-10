os :
ubuntu 18.04.5

Compile instruction :
Source.cpp:3個matrix
g++ -O2 Source.cpp -o Source -fopenmp `pkg-config --cflags --libs opencv` -lopenblas
Source.cpp:2個matrix
g++ -O2 Source_v2.cpp -o Source_v2 -fopenmp `pkg-config --cflags --libs opencv` -lopenblas

run :
./Source
./Source_v2


