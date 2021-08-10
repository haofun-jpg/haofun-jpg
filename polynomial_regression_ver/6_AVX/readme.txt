os : 
ubuntu 18.04.5

Compile instruction :
g++ -O2 AVX.cpp -o AVX -ffast-math -std=c++11 -pthread -mavx -mavx2

g++ -O2 AVX.cpp -o AVX -ffast-math -std=c++11 -pthread -msse -msse2

run :
./AVX


