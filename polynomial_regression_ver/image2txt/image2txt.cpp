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


void Matshape(cv::Mat  img, string text) {
	//img.cols 宽
	//img.rows 高
	cout << text << " (高, 寬, 通道) = ( " << img.rows << ", " << img.cols << ", " << img.channels() << " )" << endl;
}


void Store_data( Mat im ){
    	fstream file;
    	file.open("image.txt", ios::out);
        for( int i = 0; i < im.rows; i++ ){
		for( int j = 0; j < im.cols; j++ ){
			for( int k = 0; k < im.channels(); k++ ){
                          file << im.at<cv::Vec3i>(i, j)[k] << '\t';
			}
		}
		file << endl;
        }
	file.close();
}

int main() {
	string data_path = "0.jpg";
	Mat im = cv::imread( data_path, 3 );
	cout << dec;
	cout << "圖片讀取成功!!"<< " size = "  << im.size << endl;
	Matshape(im, data_path);

	cvtColor(im, im, cv::COLOR_BGR2RGB); // opencv原本是bgr表示轉成rgb
        //cout << im << endl;
	im.convertTo(im, CV_32S);

	Store_data( im );
	cout << "Generate image success !!" << endl;
}
