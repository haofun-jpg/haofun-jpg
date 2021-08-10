ubuntu 18.04.5

Compile instruction :
最快：但加最多不必要東西
g++ -O2 malloc.cpp -o malloc -ffast-math -std=c++11 -pthread -mavx -mavx2 -fopenmp
次快：一樣有malloc
g++ -O2 malloc_not_use.cpp -o malloc_not_use -ffast-math -std=c++11 -pthread -mavx -mavx2 -fopenmp
最慢：但最精簡
g++ -O2 malloc_thread.cpp -o malloc_thread -ffast-math -std=c++11 -pthread -mavx -mavx2 -fopenmp
run :
./malloc
./malloc_not_use
./malloc_thread


