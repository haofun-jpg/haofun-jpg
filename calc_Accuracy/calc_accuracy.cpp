#include<iostream>
#include<fstream>
#include <stdio.h>
#include <vector>
#include <string.h>
#include <sstream>
#include <string>
#include <time.h>
#include <cmath>
#include <iterator>
using namespace std;


string int2str(int& i) {
    string s;
    stringstream ss(s);
    ss << i;
    return ss.str();
}


void load_data( float* &data, int imagesize ,string filename ){
    float str;	
    fstream file;
    file.open(filename, ios::in);

    for( int i = 0; i < imagesize; i++ ){
	file >> data[i] ;
    }

    file.close();
}

void calcR( float* data1, float* data2, float data1_avg,int data_size, float &R_value ){
	//R^2 = 1 - u/v
	float u = 0;
	float v = 0;
	float temp_u;
	float temp_v;
	for( int i = 0; i < data_size; i++ ){
		temp_u = data1[i] - data2[i];
		u += temp_u*temp_u;
		temp_v = data1[i] - data1_avg;
		v += temp_v*temp_v;
	}
	R_value = 1 - u/v;

}

float calc_true_mean( float* data1, int imagesize ){
	float sum = 0;
	float avg = 0 ;
	for( int i = 0; i < imagesize; i++){
		sum += data1[i];
	}
	avg = sum/imagesize;
	return avg;
}

int main(){

	int imagesize = 300*300*401;
	float* data_CPP = (float*)malloc(imagesize * sizeof(float));
	float* data_regression = (float*)malloc(imagesize * sizeof(float));
	string file_cpp;
	string file_regression;

	float cpp_avg = calc_true_mean( data_CPP, imagesize );
	float R_value = 0;
	float total_R_value = 0;
	int image_num = 20;
	string image_type ="Wetland";
	
	for( int i = 10, index = 11; i < image_num; i++, index++ ){
		//load cpp
		file_cpp = "./" + image_type + "/Algorithm_data/" + int2str(index) + ".txt";
		load_data( data_CPP, imagesize, file_cpp );

		//load regression
		file_regression = "./" + image_type + "/Regression_data/" + int2str(index) + ".txt";
		load_data( data_regression, imagesize, file_regression );

		calcR(data_CPP, data_regression, cpp_avg, imagesize, R_value);
		total_R_value += R_value;
		cout << "image" << index << " accuracy :" << R_value << endl;
	}
	cout << image_type << " accuracy :" << total_R_value/image_num << endl;

}



