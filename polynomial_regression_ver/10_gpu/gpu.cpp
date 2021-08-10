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
#include <vector>
#include <CL/cl.hpp>
#include <ctime> // clock 函數所需之標頭檔
using namespace std;


string int2str(int& i) {
    string s;
    stringstream ss(s);
    ss << i;
    return ss.str();
}

void load_image( int* R_list, int* G_list, int* B_list, int size, string filename ){
	fstream file;
	
	file.open(filename.c_str(), ios::in);
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


void StoreCSV(float* image_pixel, int row, int col, int ch, string filename ) {

	filename = filename + ".csv";

	FILE* fp;
	fp = fopen(filename.c_str(), "w");
	int image_size = row*col;
	for (int k = 0; k < ch; k++ ) {
		for( int i = 0; i < image_size; i++ ){
			fprintf( fp, "%10.4f", image_pixel[i*ch+k] );
			fprintf(fp, "\n");
		}
	}

	fclose(fp);
}




string opencl_kernel = 
R"(
__kernel void vecspectrum( 
__global int* parameter_list, 
__global int* R_list,  
__global int* G_list, 
__global int* B_list, 
__global float* w, 
__global float* image_pixel ){



        int i = get_global_id(0);
	int vRGB[36];
	int parameter = parameter_list[0];
	int wave_size = parameter_list[1];


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

	for( int j = 0; j < wave_size; j++ ){
		image_pixel[ i*wave_size + j ] = vRGB[0]*w[0+j*parameter] + vRGB[1]*w[1+j*parameter] + vRGB[2]*w[2+j*parameter] + vRGB[3]*w[3+j*parameter] + 
						 vRGB[4]*w[4+j*parameter] + vRGB[5]*w[5+j*parameter] + vRGB[6]*w[6+j*parameter] + vRGB[7]*w[7+j*parameter] + 
						 vRGB[8]*w[8+j*parameter] + vRGB[9]*w[9+j*parameter] + vRGB[10]*w[10+j*parameter] + vRGB[11]*w[11+j*parameter] + 
						 vRGB[12]*w[12+j*parameter] + vRGB[13]*w[13+j*parameter] + vRGB[14]*w[14+j*parameter] + vRGB[15]*w[15+j*parameter] + 
						 vRGB[16]*w[16+j*parameter] + vRGB[17]*w[17+j*parameter] + vRGB[18]*w[18+j*parameter] + vRGB[19]*w[19+j*parameter] + 
						 vRGB[20]*w[20+j*parameter] + vRGB[21]*w[21+j*parameter] + vRGB[22]*w[22+j*parameter] + vRGB[23]*w[23+j*parameter] +
						 vRGB[24]*w[24+j*parameter] + vRGB[25]*w[25+j*parameter] + vRGB[26]*w[26+j*parameter] + vRGB[27]*w[27+j*parameter] + 
						 vRGB[28]*w[28+j*parameter] + vRGB[29]*w[29+j*parameter] + vRGB[30]*w[30+j*parameter] + vRGB[31]*w[31+j*parameter] + 
						 vRGB[32]*w[32+j*parameter] + vRGB[33]*w[33+j*parameter] + vRGB[34]*w[34+j*parameter] + vRGB[35]*w[35+j*parameter];
	}



   }

)";

int main( int argc , char **argv ) {
   string filename = "../1920x1080.txt";
   int row = 1920;
   int col = 1080;
   unsigned short parameter = 36;
   unsigned short wave_size = 401;
   int parameter_list[2];
   parameter_list[0] = parameter; // parameter
   parameter_list[1] = wave_size;// wave_size
   int image_size = row * col;
   int* R_list = (int*)aligned_alloc(sizeof(int),image_size * sizeof(int));
   int* G_list = (int*)aligned_alloc(sizeof(int),image_size * sizeof(int));

   float* rrr_list = (float*)aligned_alloc(sizeof(float),image_size * sizeof(float));

   int* B_list = (int*)aligned_alloc(sizeof(int),image_size * sizeof(int));

   float* image_pixel = (float*)aligned_alloc( sizeof(float), image_size * wave_size * sizeof(float));

   load_image( R_list, G_list, B_list, image_size, filename );
   cout << "end load_image" << endl;

   float* w = (float*)aligned_alloc(sizeof(float), wave_size * parameter * sizeof(float));
   load_wave(w, wave_size, parameter);
   cout << "end load_wave" << endl;

   //set Platform(cpu)
   std::vector<cl::Platform> platforms;
   cl::Platform::get(&platforms);

   if (platforms.empty()) {
   	printf("No platforms!\n");		
		return 1;
	}
   cl::Platform platform = platforms[0];

   //set Device(gpu)
   std::vector<cl::Device> Devices;
   platform.getDevices(CL_DEVICE_TYPE_GPU, &Devices);
   if (Devices.empty()) {
      printf("No Devices!\n");
      return 1;
   }
   cl::Device device = Devices[0];
   std::cout << "Device : " << device.getInfo<CL_DEVICE_NAME>() <<"\n";
   
   //connection cpu&gpu
   cl::Context context({device});

   //compile opencl kernel
   cl::Program program(context, opencl_kernel);
   if (program.build({device}) != CL_SUCCESS) {
      printf("Fail to build\n");
      return 1;
   }

   //set buffer
   cl::Buffer buffer_parameter(context, CL_MEM_READ_WRITE, 2 * sizeof(int)); 
   cl::Buffer buffer_R(context, CL_MEM_READ_WRITE, image_size * sizeof(int)); 
   cl::Buffer buffer_G(context, CL_MEM_READ_WRITE, image_size * sizeof(int));
   cl::Buffer buffer_B(context, CL_MEM_READ_WRITE, image_size * sizeof(int));
   cl::Buffer buffer_w(context, CL_MEM_READ_WRITE, wave_size * parameter * sizeof(float) );
   cl::Buffer buffer_spectral(context, CL_MEM_READ_WRITE, image_size * wave_size * sizeof(float));

   //sent data to device
   cl::CommandQueue queue(context, device);

   //set kernal arg

   cl::Kernel vec_spectrum(program, "vecspectrum");
   vec_spectrum.setArg(0, buffer_R);
   vec_spectrum.setArg(1, buffer_spectral);

   vec_spectrum.setArg(0, buffer_parameter);
   vec_spectrum.setArg(1, buffer_R);
   vec_spectrum.setArg(2, buffer_G);
   vec_spectrum.setArg(3, buffer_B);
   vec_spectrum.setArg(4, buffer_w);
   vec_spectrum.setArg(5, buffer_spectral);

   queue.enqueueWriteBuffer(buffer_parameter, CL_FALSE, 0, 2 * sizeof(int), parameter_list);
   queue.enqueueWriteBuffer(buffer_R, CL_FALSE, 0, image_size * sizeof(int), R_list);
   queue.enqueueWriteBuffer(buffer_G, CL_FALSE, 0, image_size * sizeof(int), G_list);
   queue.enqueueWriteBuffer(buffer_B, CL_FALSE, 0, image_size * sizeof(int), B_list);
   queue.enqueueWriteBuffer(buffer_w, CL_FALSE, 0, wave_size * parameter * sizeof(float), w);

   queue.enqueueNDRangeKernel(vec_spectrum, cl::NullRange, cl::NDRange(image_size), cl::NullRange);
   queue.enqueueReadBuffer(buffer_spectral, CL_FALSE, 0, image_size * wave_size * sizeof(float), image_pixel);

   queue.finish();

}
