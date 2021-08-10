#include<iostream>
#include<fstream>
#include <stdio.h>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <math.h>
#include <iterator>
#include <omp.h> 
#include <sstream>


using namespace std;



string int2str(int& i) {
    string s;
    stringstream ss(s);
    ss << i;
    return ss.str();
}

void load_image( int* &R_list, int* &G_list, int* &B_list, int size ){
	fstream file;
	file.open("..//1920x1080.txt", ios::in);
	for( int i = 0; i < size ; i++ ){
		file >> R_list[i];
                //cout << R_list[i] << endl;
		file >> G_list[i];
                //cout << G_list[i] << endl;
		file >> B_list[i];
                //cout << B_list[i] << endl;
	}
	file.close();
}

void load_wave(float**& w, int parameter, int nm) {
    float temp;
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

void store_data( float** image_pixel, int image_size, int wave_size ){

	fstream file;
	file.open("regression_output.txt", ios::out);
	
	for( int i = 0; i < wave_size ; i++ ){
		for( int j = 0; j < image_size; j++ ){
			file << image_pixel[j][i] << '\t';
		}
		file << endl ;
	}

	file.close();

/*
	FILE* fp;
	fp = fopen("regression_output.txt", "w");
	for( int i = 0; i < wave_size ; i++ ){
		for( int j = 0; j < image_size; j++ ){
			fprintf( fp, "%10.4f", image_pixel[j][i] );
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
*/
}


void image_2_spectrum( int* R_list, int* G_list, int* B_list, float** w, float** &image_pixel, int image_size, int wave_size ){
	int R;
	int G;
	int B;
	for( int i = 0; i < image_size; i++ ){
		R = R_list[i];
		G = G_list[i];
		B = B_list[i];
		for( int j = 0; j < wave_size; j++ ){
        		image_pixel[i][j] = w[j][0] + w[j][1] +
            		//^1
            		w[j][2] * R + w[j][3] * G + w[j][4] * B +
            		//^2
            		w[j][5] * (R*R) + w[j][6] * (R*G) + w[j][7] * (R*B) + w[j][8] * (G*G) + w[j][9] * (G*B) + w[j][10] * (B*B) +
            		//^3
            		w[j][11] * (R*R*R) + w[j][12] * (R*R*G) + w[j][13] * (R*R*B) + w[j][14] * (R*G*G) + w[j][15] * (R*G*B) +
            		w[j][16] * (R*B*B) + w[j][17] * (G*G*G) + w[j][18] * (G*G*B) + w[j][19] * (G*B*B) + w[j][20] * (B*B*B) +
            		//^4
            		w[j][21] * (R*R*R*R) + w[j][22] * (R*R*R*G) + w[j][23] * (R*R*R*B) + w[j][24] * (R*R*G*G) + w[j][25] * (R*R*G*B) +
            		w[j][26] * (R*R*B*B) + w[j][27] * (R*G*G*G) + w[j][28] * (R*G*G*B) + w[j][29] * (R*G*B*B) + w[j][30] * (R*B*B*B) +
            		w[j][31] * (G*G*G*G) + w[j][32] * (G*G*G*B) + w[j][33] * (G*G*B*B) + w[j][34] * (G*B*B*B) + w[j][35] * (B*B*B*B) +
			//^5
            		w[j][21] * (R*R*R*R*R) + w[j][22] * (R*R*R*R*G) + w[j][23] * (R*R*R*R*B) + w[j][24] * (R*R*R*G*G) + w[j][25] * (R*R*R*G*B) +
            		w[j][26] * (R*R*R*B*B) + w[j][27] * (R*R*G*G*G) + w[j][28] * (R*R*G*G*B) + w[j][29] * (R*R*G*B*B) + w[j][30] * (R*R*B*B*B) +
            		w[j][31] * (R*G*G*G*G) + w[j][32] * (R*G*G*G*B) + w[j][33] * (R*G*G*B*B) + w[j][34] * (R*G*B*B*B) + w[j][35] * (R*B*B*B*B) +
      			w[j][31] * (G*G*G*G*G) + w[j][32] * (G*G*G*G*B) + w[j][33] * (G*G*G*B*B) + w[j][34] * (G*G*B*B*B) + w[j][35] * (G*B*B*B*B) +
      			w[j][31] * (B*B*B*B*B) ;


		}
	}
	
}


int main() {
    int row = 1920;
    int col = 1080;
/*
    cout << "-please input row & col-" << endl;
    cout << "input row :";
    cin >> row;
    cout << "input col :";
    cin >> col;
*/

    int wave_size = 401;
    int parameter = 36; //4 degree have 36 parameter
    float** w = (float**)malloc(wave_size * sizeof(float*));

    for (int i = 0; i < wave_size; i++)
        w[i] = (float*)malloc(parameter * sizeof(float));

    for (int i = 0; i < 401; i++) {
        load_wave(w, parameter, i);
    }
    cout << "end load_wave" << endl;

    int image_size = row * col;

    float** image_pixel = (float**)malloc( image_size * sizeof(float*));
    for( int i = 0; i < image_size; i++ ){
	image_pixel[i] = (float*)malloc(wave_size * sizeof(float));
    }

    int* R_list = (int*)malloc(image_size * sizeof(int));
    int* G_list = (int*)malloc(image_size * sizeof(int));
    int* B_list = (int*)malloc(image_size * sizeof(int));
    load_image( R_list, G_list, B_list, image_size );
    double START, END;
    int index = 0;
    START = clock();
    image_2_spectrum( R_list, G_list, B_list, w, image_pixel, image_size, wave_size );

    END = clock();
    cout << "pixel_2_spectrum所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;
    cout << "end pixel_2_spectrum" << endl;
    //store_data( image_pixel, image_size, wave_size );
    cout << "end store_data" << endl;
   
}
