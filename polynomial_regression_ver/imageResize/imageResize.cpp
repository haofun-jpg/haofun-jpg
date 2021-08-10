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


int main() {
	string data_path = "100x100.jpg";
	Mat im = cv::imread( data_path, 3 );
	cout << "圖片讀取成功!!"<< " size = "  << im.size << endl;
	Matshape(im, data_path);

	cvtColor(im, im, cv::COLOR_BGR2RGB); // opencv原本是bgr表示轉成rgb

	resize(im,im,Size(im.cols/10,im.rows/10),0,0,INTER_LINEAR); //縮放
	cv::imwrite("10x10.jpg", im);//store


}
