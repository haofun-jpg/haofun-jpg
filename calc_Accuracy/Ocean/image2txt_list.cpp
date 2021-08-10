#include<iostream>
#include<fstream>
#include <stdio.h>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <string>
#include <opencv2/opencv.hpp>
#include <time.h>
#include <math.h>
#include <iterator>
#include <omp.h> 

using namespace std;
using namespace cv;
string int2str(int& i) {
    string s;
    stringstream ss(s);
    ss << i;
    return ss.str();
}

void Matshape(cv::Mat  img, string text) {
	//img.cols 宽
	//img.rows 高
	cout << text << " (高, 寬, 通道) = ( " << img.rows << ", " << img.cols << ", " << img.channels() << " )" << endl;
}


void Store_data( Mat im, string output_path ){
    	fstream file;
    	file.open(output_path, ios::out);
        for( int i = 0; i < im.rows; i++ ){
		for( int j = 0; j < im.cols; j++ ){
			for( int k = 0; k < im.channels(); k++ ){
                          file << im.at<cv::Vec3i>(i, j)[k] << '\t';
			}
			file << endl;
		}
        }
	file.close();
}

int main() {

	string index;
	string input_path;
	string output_path;
	for( int i = 0, index = 1; i < 10; i++,index++ ){
		input_path = int2str(index) + ".png";
		Mat im = cv::imread( input_path, 3 );

		Matshape(im, input_path);
		cvtColor(im, im, cv::COLOR_BGR2RGB); // opencv原本是bgr表示轉成rgb
		im.convertTo(im, CV_32S);
		output_path = "./imtxt/" + int2str(index) + ".txt";
		Store_data( im, output_path );
	}
	cout << "Generate image success !!" << endl;
}
