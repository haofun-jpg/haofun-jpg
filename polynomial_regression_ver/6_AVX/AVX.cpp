#include <iostream>
#include <thread>
#include<fstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <iterator>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <sstream>
#include <pthread.h>
#include <unistd.h>
#include <immintrin.h>
#include <omp.h> 
#include <ctime> // clock 函數所需之標頭檔
using namespace std;

string int2str(int& i) {
    string s;
    stringstream ss(s);
    ss << i;
    return ss.str();
}

void load_image( unsigned short* R_list, unsigned short* G_list, unsigned short* B_list, int size ){
	fstream file;
	file.open("..//1920x1080.txt", ios::in);
	for( int i = 0; i < size ; i++ ){
		file >> R_list[i];
		file >> G_list[i];
		file >> B_list[i];
	}
	file.close();
}


void load_wave(float**& w, int parameter, int nm) {
    int file_nm = nm + 380;
    string directory = "./C_degree4/";
    string filename = directory + int2str(file_nm) + "nm_degree4.txt";
    ifstream file(filename.c_str());
    if (!file) {
        cout << "無法讀入檔案\n";
    }

    for (int j = 0; j < parameter; j++) {
        file >> w[nm][j];
    }

    file.close();
}

void StoreMAE_style(float** image_pixel, int row, int col, int ch ) {


	FILE* fp;
	fp = fopen("regression_output.txt", "w");

	for (int k = 0; k < ch; k++ ) {
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++ ) {
				fprintf( fp, "%10.4f", image_pixel[j + i * col][k] );
				fprintf(fp, "\n");
			}
		}
	}
	fclose(fp);
}

void spectrum( int start, int end, unsigned short wave_size, unsigned short* &R_list,  unsigned short* &G_list, unsigned short* &B_list,
 float** &RGB_var, float** &w, float** &image_pixel ){
   float CSE0;
   float CSE1;
   float CSE2;
   float CSE3;
   float CSE4;
   float CSE5;
   float CSE6;
   float CSE7;
   float CSE8;
   float CSE9;
   float CSE10;
   float CSE11;
   float CSE12;
   float CSE13;
   float CSE14;
   float CSE15;
   float temp ;
   float sum[8] __attribute__((aligned(sizeof(__m256))));
   //float* vRGB[36] __attribute__((aligned(sizeof(__m256))));
   float* vRGB = (float*)aligned_alloc(sizeof(__m256), 36 * sizeof(float));
   register __m256 vRGB_0, vRGB_1, vRGB_2, vRGB_3, vRGB_4;
   register __m256 w0, w1, w2, w3, w4;
   float* pRGB;
   float* pw;
   //#pragma omp parallel for
   for( int i = start; i < end; i++ ){
     vRGB[0] = 1;
     vRGB[1] = 1;
     vRGB[2] = R_list[i];
     vRGB[3] = G_list[i];
     vRGB[4] = B_list[i];
     vRGB[5] = R_list[i]*R_list[i];  // CSE12
     vRGB[6] = R_list[i]*G_list[i]; // CSE13
     vRGB[7] = R_list[i]*B_list[i]; // CSE14
     vRGB[8] = G_list[i]*G_list[i]; // CSE10
     vRGB[9] = G_list[i]*B_list[i]; // CSE11
     vRGB[10] = B_list[i]*B_list[i]; // CSE15
     vRGB[11] = vRGB[5]*R_list[i];  // CSE5
     vRGB[12] = vRGB[5]*G_list[i];  // CSE9
     vRGB[13] = vRGB[5]*B_list[i];  // CSE4
     vRGB[14] = vRGB[6]*G_list[i];  // CSE7
     vRGB[15] = vRGB[6]*B_list[i];  // CSE0
     vRGB[16] = vRGB[7]*B_list[i];  // CSE1
     vRGB[17] = vRGB[8]*G_list[i];  // CSE2
     vRGB[18] = vRGB[8]*B_list[i];  // CSE8
     vRGB[19] = vRGB[9]*B_list[i];  // CSE6
     vRGB[20] = vRGB[10]*B_list[i]; // CSE3
     vRGB[21] = vRGB[11]*R_list[i];
     vRGB[22] = vRGB[11]*G_list[i];
     vRGB[23] = vRGB[11]*B_list[i];
     vRGB[24] = vRGB[12]*G_list[i];
     vRGB[25] = vRGB[12]*B_list[i];
     vRGB[26] = vRGB[13]*B_list[i];
     vRGB[27] = vRGB[14]*G_list[i];
     vRGB[28] = vRGB[14]*B_list[i];
     vRGB[29] = vRGB[15]*B_list[i];
     vRGB[30] = vRGB[16]*B_list[i];
     vRGB[31] = vRGB[17]*G_list[i];
     vRGB[32] = vRGB[17]*B_list[i];
     vRGB[33] = vRGB[18]*B_list[i];
     vRGB[34] = vRGB[19]*B_list[i];
     vRGB[35] = vRGB[20]*B_list[i];
     pRGB = vRGB;
     vRGB_0 = _mm256_load_ps(pRGB);
     pRGB+=8;
     vRGB_1 = _mm256_load_ps(pRGB);
     pRGB+=8;
     vRGB_2 = _mm256_load_ps(pRGB);
     pRGB+=8;
     vRGB_3 = _mm256_load_ps(pRGB);
     pRGB+=8;
     vRGB_4 = _mm256_load_ps(pRGB);
     pRGB+=4;
     for( int j = 0; j < wave_size; j++ ){
         pw = w[j];
         w0 = _mm256_load_ps(pw);
	 pw+=8;
	 w1 = _mm256_load_ps(pw);
	 pw+=8;
	 w2 = _mm256_load_ps(pw);
	 pw+=8;
	 w3 = _mm256_load_ps(pw);
	 pw+=8;
	 w4 = _mm256_load_ps(pw);
	 pw+=4;
         w0 = w0 * vRGB_0;
         w1 = w1 * vRGB_1;
	 w2 = w2 * vRGB_2;
	 w3 = w3 * vRGB_3;
	 w4 = w4 * vRGB_4;
	 w0 = w0 + w1 + w2 + w3 + w4;
	 _mm256_store_ps(sum, w0);
	 image_pixel[i][j] = sum[0]+sum[1]+sum[2]+sum[3]+sum[4]+sum[5]+sum[6]+sum[7];
     }

   }
  
}

int main(){
   int row = 1920;
   int col = 1080;
   int numthread = 18;
   unsigned short parameter = 36;
   unsigned short wave_size = 401;
   int image_size = row * col;

  // 儲存時間用的變數
  clock_t start, end;
  double cpu_time_used;
  // 計算開始時間
  start = clock();

   unsigned short* R_list = (unsigned short*)aligned_alloc(sizeof(unsigned short),image_size * sizeof(unsigned short));
   unsigned short* G_list = (unsigned short*)aligned_alloc(sizeof(unsigned short),image_size * sizeof(unsigned short));
   unsigned short* B_list = (unsigned short*)aligned_alloc(sizeof(unsigned short),image_size * sizeof(unsigned short));


   float** image_pixel = (float**)aligned_alloc( sizeof(float), image_size * sizeof(float*));
   for( int i = 0; i < image_size; i++ ){
	image_pixel[i] = (float*)aligned_alloc( sizeof(float) ,wave_size * sizeof(float));
   }

   float** RGB_var = (float**)aligned_alloc(sizeof(float), image_size * sizeof(float*));
   for( int i = 0; i < image_size; i++ ) 
	RGB_var[i] = (float*)aligned_alloc(sizeof(float), parameter * sizeof(float));
  // 計算結束時間
  end = clock();
  // 計算實際花費時間
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("Time = %f\n", cpu_time_used);

   load_image( R_list, G_list, B_list, image_size );
   cout << "end load_image" << endl;

   float** w = (float**)aligned_alloc(sizeof(float), wave_size * sizeof(float*));
   for (int i = 0; i < wave_size; i++)
	w[i] = (float*)aligned_alloc(sizeof(float), parameter * sizeof(float));

   for (int i = 0; i < wave_size; i++) {
       load_wave(w, parameter, i);
   }
   cout << "end load_wave" << endl;


   std::thread myThreads[numthread];

   unsigned int part = image_size / numthread;

   for( int i = 0; i < numthread; i++ ){
	   myThreads[i] = std::thread(spectrum, part*i,part*(i+1), wave_size, ref(R_list), ref(G_list), ref(B_list), 
	   ref(RGB_var), ref(w), ref(image_pixel));
   }


   for( unsigned int i = 0; i < numthread; i++ ){
	myThreads[i].join();
   }
   //Store_data
   //StoreMAE_style( image_pixel, row, col , wave_size);

}
