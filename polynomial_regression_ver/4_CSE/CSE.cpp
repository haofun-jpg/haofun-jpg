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
		file >> G_list[i];
		file >> B_list[i];
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
	int R;
	int G;
	int B;
	for( int i = 0; i < image_size; i++ ){
		for( int j = 0; j < wave_size; j++ ){
		R = R_list[i];
		G = G_list[i];
		B = B_list[i];
		CSE13=(R*G);
		CSE0=(CSE13*B);
		CSE14=(R*B);
		CSE1=(CSE14*B);
		CSE10=(G*G);
		CSE2=(CSE10*G);
		CSE15=(B*B);
		CSE3=(CSE15*B);
		CSE12=(R*R);
		CSE4=(CSE12*B);
		CSE5=(CSE12*R);
		CSE11=(G*B);
		CSE6=(CSE11*B);
		CSE7=(CSE13*G);
		CSE8=(CSE10*B);
		CSE9=(CSE12*G);
		image_pixel[i][j] = (((((((((((((((((((((((((((((((((((w[j][0]+w[j][1])+
		//^1
		(w[j][2]*R))+(w[j][3]*G))+(w[j][4]*B))+
		//^2
		(w[j][5]*CSE12))+(w[j][6]*CSE13))+(w[j][7]*CSE14))+(w[j][8]*CSE10))+(w[j][9]*CSE11))+(w[j][10]*CSE15))+
		//^3
		(w[j][11]*CSE5))+(w[j][12]*CSE9))+(w[j][13]*CSE4))+(w[j][14]*CSE7))+(w[j][15]*CSE0))+
		(w[j][16]*CSE1))+(w[j][17]*CSE2))+(w[j][18]*CSE8))+(w[j][19]*CSE6))+(w[j][20]*CSE3))+
		//^4
		(w[j][21]*(CSE5*R)))+(w[j][22]*(CSE5*G)))+(w[j][23]*(CSE5*B)))+(w[j][24]*(CSE9*G)))+(w[j][25]*(CSE9*B)))+
		(w[j][26]*(CSE4*B)))+(w[j][27]*(CSE7*G)))+(w[j][28]*(CSE7*B)))+(w[j][29]*(CSE0*B)))+(w[j][30]*(CSE1*B)))+
		(w[j][31]*(CSE2*G)))+(w[j][32]*(CSE2*B)))+(w[j][33]*(CSE8*B)))+(w[j][34]*(CSE6*B)))+(w[j][35]*(CSE3*B)));
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
