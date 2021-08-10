ubuntu 18.04.5

Compile instruction :
g++ -O2 gpu.cpp -o gpu -ffast-math -std=c++11 -lOpenCL -mavx -mavx2 

run :
./gpu
