#include <iostream>
#include <thread>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <iterator>
#include <sstream>
#include <pthread.h>
#include <unistd.h>
#include <omp.h> 
#include <ctime> // clock 函數所需之標頭檔
#include <arm_neon.h> //neon
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


void load_wave(float*& w, int wave_size, int parameter) {

    for( int i = 0; i < wave_size; i++ ){

      int file_nm = i + 380;
      string directory = "./C_degree4/";
      string filename = directory + int2str(file_nm) + "nm_degree4.txt";
      ifstream file(filename.c_str());
      if (!file) {
          cout << "無法讀入檔案\n";
      }

      for( int j = 0; j < parameter; j++ ){
	  file >> w[ i*parameter + j ];
      }

      file.close();

    }
}


void StoreData(float* image_pixel, int row, int col, int ch ) {


	FILE* fp;
	fp = fopen("regression_ch_out.txt", "w");
	int image_size = row*col;
	for( int k = 0; k < ch; k++ ){
		fprintf( fp, "val(:,:,%d) =\n\n", k + 1 );
		for( int i = 0; i < image_size; i++ ){
			fprintf( fp, "%10.4f", image_pixel[i*ch+k] );
			if( (i+1)%col == 0 ){
				fprintf(fp, "\n");
			}
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}


void StoreMAE_style(float* image_pixel, int row, int col, int ch ) {


	FILE* fp;
	fp = fopen("regression_output.txt", "w");
	int image_size = row*col;
	for (int k = 0; k < ch; k++ ) {
		for( int i = 0; i < image_size; i++ ){
			fprintf( fp, "%10.4f", image_pixel[i*ch+k] );
			fprintf(fp, "\n");
		}
	}
	fclose(fp);
}

void spectrum( int start, int end, unsigned short wave_size, unsigned short parameter,unsigned short* &R_list,  unsigned short* &G_list, unsigned short* &B_list, float* &w, float* &image_pixel, int index ){
   
   float* vRGB = (float*)malloc(sizeof(float)*36);
   float* pRGB;
   float* pw;
   float32x4_t vRGB_0, vRGB_1, vRGB_2, vRGB_3, vRGB_4, vRGB_5, vRGB_6, vRGB_7, vRGB_8;
   float32x4_t w0, w1, w2, w3, w4, w5, w6, w7, w8;
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
     vRGB_0 = vld1q_f32(pRGB);
     pRGB+=4;
     vRGB_1 = vld1q_f32(pRGB);
     pRGB+=4;
     vRGB_2 = vld1q_f32(pRGB);
     pRGB+=4;
     vRGB_3 = vld1q_f32(pRGB);
     pRGB+=4;
     vRGB_4 = vld1q_f32(pRGB);
     pRGB+=4;
     vRGB_5 = vld1q_f32(pRGB);
     pRGB+=4;
     vRGB_6 = vld1q_f32(pRGB);
     pRGB+=4;
     vRGB_7 = vld1q_f32(pRGB);
     pRGB+=4;
     vRGB_8 = vld1q_f32(pRGB);
     pRGB+=4;
     for( int j = 0; j < wave_size; j++ ){
         //pw = &w[j*parameter];
	 pw = (w+j*parameter);
         w0 = vld1q_f32(pw);
	 pw+=4;
	 w1 = vld1q_f32(pw);
	 pw+=4;
	 w2 = vld1q_f32(pw);
	 pw+=4;
	 w3 = vld1q_f32(pw);
	 pw+=4;
	 w4 = vld1q_f32(pw);
	 pw+=4;
	 w5 = vld1q_f32(pw);
	 pw+=4;
	 w6 = vld1q_f32(pw);
	 pw+=4;
	 w7 = vld1q_f32(pw);
	 pw+=4;
	 w8 = vld1q_f32(pw);
	 pw+=4;
         w0 = w0 * vRGB_0;
         w1 = w1 * vRGB_1;
	 w2 = w2 * vRGB_2;
	 w3 = w3 * vRGB_3;
	 w4 = w4 * vRGB_4;
         w5 = w5 * vRGB_5;
         w6 = w6 * vRGB_6;
	 w7 = w7 * vRGB_7;
	 w8 = w8 * vRGB_8;
	 w0 = w0 + w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8;
	 
	 image_pixel[ index*wave_size + j ] = w0[0]+w0[1]+w0[2]+w0[3];
     }

   }
  
}

int main(){
   int row = 1920;
   int col = 1080;
   int numthread = 8;
   unsigned short parameter = 36;
   unsigned short wave_size = 401;
   int image_size = row * col;

   unsigned short* R_list = (unsigned short*)malloc(image_size * sizeof(unsigned short));
   unsigned short* G_list = (unsigned short*)malloc(image_size * sizeof(unsigned short));
   unsigned short* B_list = (unsigned short*)malloc(image_size * sizeof(unsigned short));

   long int image_pixel_size = image_size * wave_size;
   //1920x1080 too big
   float* image_pixel = (float*)malloc(image_pixel_size/2 * sizeof(float));
   float* image_pixel_other = (float*)malloc(image_pixel_size/2 + image_pixel_size %  numthread * sizeof(float));

   load_image( R_list, G_list, B_list, image_size );
   cout << "end load_image" << endl;
   
   float* w = (float*)malloc(wave_size * parameter * sizeof(float));


   load_wave(w, wave_size, parameter);
   cout << "end load_wave" << endl;


   if( image_size < numthread ){
      //image_size < numthread
      //thread num = 1 , calc 0-image_size
	 std::thread myThread;
	 std::thread myThread_other;
	 myThread = std::thread(spectrum, 0, image_size/2, wave_size, parameter, ref(R_list), ref(G_list), ref(B_list), ref(w), ref(image_pixel), 0 );
	 myThread_other = std::thread(spectrum, image_size/2, image_size, wave_size, parameter, ref(R_list), ref(G_list), ref(B_list), ref(w), ref(image_pixel_other), image_size/2);
	 myThread.join();
	 myThread_other.join();
      
   }
   else{
      //image_size >= numthread
      //thread num = numthread , part(i) calc part*i-part*(i+1)
      if( image_size % numthread == 0 ){
	 
	 std::thread myThreads[numthread];
	 unsigned int part = image_size / numthread;
	 
	 int i = 0;
	 for( ; i < numthread/2; i++ ){
	    myThreads[i] = std::thread(spectrum, part*i,part*(i+1), wave_size, parameter, ref(R_list), ref(G_list), ref(B_list), ref(w), ref(image_pixel), part*i );
	 }
	 for( ; i < numthread; i++  ){
	    myThreads[i] = std::thread(spectrum, part*i,part*(i+1), wave_size, parameter, ref(R_list), ref(G_list), ref(B_list), ref(w), ref(image_pixel_other), part*i - ( numthread/2 * part ) );
	 }
	    
	 i = 0;
	 for( ; i < numthread; i++  ){
	    myThreads[i].join();
	 }
      }p
      //thread num = numthread , part(i) calc part*i-part*(i+1)
      //Ex_Threads calc part*(numthread)-image_size
      else{
	 std::thread myThreads[numthread];
	 std::thread Ex_Threads;
	 unsigned int part = image_size / numthread;
	 
	 int i = 0;
	    
	 for( ; i < numthread/2; i++ ){
	    myThreads[i] = std::thread(spectrum, part*i,part*(i+1), wave_size, parameter, ref(R_list), ref(G_list), ref(B_list), ref(w), ref(image_pixel), part*i );
	 }
	 for( ; i < numthread; i++  ){
	    myThreads[i] = std::thread(spectrum, part*i,part*(i+1), wave_size, parameter, ref(R_list), ref(G_list), ref(B_list), ref(w), ref(image_pixel_other), part*i - ( numthread/2 * part ) );
	 }
	 Ex_Threads = std::thread(spectrum, part*numthread, image_size, wave_size, parameter, ref(R_list), ref(G_list), ref(B_list), ref(w), ref(image_pixel_other), part*numthread - ( numthread/2 * part ) );
	    
	 for( int i = 0; i < numthread; i++ ){
	    myThreads[i].join();
	 }
	 Ex_Threads.join();


      }
   
   }
   

   //StoreData(image_pixel, row, col, wave_size );
   //StoreMAE_style(image_pixel, row, col, wave_size );

}
