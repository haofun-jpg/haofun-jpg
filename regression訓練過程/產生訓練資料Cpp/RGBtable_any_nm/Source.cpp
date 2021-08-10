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
extern "C"
{
#include <cblas.h> 
}



using namespace std;
using namespace cv;
int pic_rows = 0;
int pic_cols = 0;
int pic_channel = 0;
int pic_element = 0;
int data_id = 0;


float* list_pic_2d;
float* list_mul_result;


void print_vector_2D(vector<vector<float> > );
void print_vector_3D(vector<vector< vector<float> > > );
void reshape_3d_to_2dtranspose(float*, int, int, int );
void reshape_2d_to_3dtranspose(int, int, int, float*&);
void matrix_mul_transpose(vector<vector<float> > , vector<vector<float> > , vector<vector<float> >&);
void print_list_2D(float* , int , int );
void print_list_3D(float*, int, int, int);

inline int LoadData(string fileName, float* & vecData, int Rows = 0, int Cols = 0, int Chns = 0)
{
	ifstream fin( fileName.c_str() );
	if (!fin) {
		cout << "無法讀入檔案\n";
		return 0;
	}
	for (int i = 0; i < Rows * Cols; i++ ) {
		fin >> vecData[i];
	}
	fin.close();
	return 0;
}

void Matshape(cv::Mat  img, string text) {
	//img.cols 宽
	//img.rows 高
	cout << text << "(高, 寬, 通道) = ( " << img.rows << ", " << img.cols << ", " << img.channels() << " )" << endl;
}

string int2str(int& i) {
	string s;
	stringstream ss(s);
	ss << i;
	return ss.str();
}

void StoreData(float* pic_spectrum, int row, int col, int ch, Mat table, int filenumber ) {

	string filename = int2str(filenumber) + "nm.csv";

	FILE* fp;
	fp = fopen(filename.c_str(), "a");

	int k = filenumber-380 ; //波長此為380nm

	for (int i = 0; i < row; i++ ) {
		for (int j = 0; j < col; j++ ) {
			fprintf(fp, "%d ", data_id);
			data_id++;
			fprintf(fp, ",");
			fprintf(fp, "%d ", table.at<cv::Vec3i>(i, j)[0]);
			fprintf(fp, ",");
			fprintf(fp, "%d ", table.at<cv::Vec3i>(i, j)[1]);
			fprintf(fp, ",");
			fprintf(fp, "%d", table.at<cv::Vec3i>(i, j)[2]);
			fprintf(fp, ",");
			fprintf(fp, "%10.4f", pic_spectrum[j + i * col + k * row * col]);
			fprintf(fp, "\n");
		}
	}
	
	
	fclose(fp);
}

void Matlab_RGB2XYZ(cv::Mat im, float* & list_pic) {
	/*----------------------------
	* 功能 : 做出Matlab_RGB2XYZ如下
	* Matlab code :
	* XYZ_D65 = rgb2xyz(pic);
	*----------------------------
	*/
	int row_XYZ = 3;
	int col_XYZ = 3;
	float* list_Ma = (float*)malloc(row_XYZ * col_XYZ * sizeof(float) );
	list_Ma[0] = 0.4124564, list_Ma[1] = 0.3575761, list_Ma[2] = 0.1804375;
	list_Ma[3] = 0.2126729, list_Ma[4] = 0.7151522, list_Ma[5] = 0.0721750;
	list_Ma[6] = 0.0193339, list_Ma[7] = 0.1191920, list_Ma[8] = 0.9503041;
	im.convertTo(im, CV_32F);
	im = im / 255;
	#pragma omp parallel for
	for (int k = 0; k < im.channels(); k++) {
		#pragma omp parallel for
		for (int i = 0; i < im.rows; i++) {
			for (int j = 0; j < im.cols; j++) {
				if (im.at<cv::Vec3f>(i, j)[k] < 0.04045)
					list_pic[j + i * im.cols + k * pic_rows * pic_cols] = ( im.at<cv::Vec3f>(i, j)[k] / 12.92 ) *100;
				else
					list_pic[j + i * im.cols + k * pic_rows * pic_cols] = ( pow(((im.at<cv::Vec3f>(i, j)[k] + 0.055) / 1.055), 2.4) )*100;
			}
		}
	}
	reshape_3d_to_2dtranspose(list_pic, pic_rows, pic_cols, pic_channel );
	int m = row_XYZ;
	int n = pic_rows * pic_cols;
	int k = col_XYZ;
	//A = m x k, B = k x n, C = m x n
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,m, n, k, 1.0, list_Ma, k, list_pic_2d, n, 0.0, list_mul_result, n);
	//mkl_free(list_Ma);
	reshape_2d_to_3dtranspose(pic_rows, pic_cols, pic_channel, list_pic);

}

 void Create_XYZ(float* list_XYZ, int row_XYZ, int col_XYZ, float* & XYZ_D65 ) {
	/*----------------------------
	 * 功能 : 做出XYZ Matlab code如下
	 * Matlab code :
	 * XYZ = reshape(  (inv(Ma)*diag( (Ma*White_light)./(Ma*White_D65) )*Ma*reshape(XYZ_D65,[],3)')' ,size(pic,1), size(pic,2), 3);
	 *----------------------------
	 */

	reshape_3d_to_2dtranspose(XYZ_D65, pic_rows, pic_cols, pic_channel );
	int m = row_XYZ;
	int n = pic_rows * pic_cols;
	int k = col_XYZ;
	//A = m x k, B = k x n, C = m x n
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, list_XYZ, k, list_pic_2d, n, 0.0, list_mul_result, n);
	reshape_2d_to_3dtranspose( pic_rows, pic_cols, pic_channel, XYZ_D65);

}

 void Extend_One(float* list_pic, float* & Extend1) {
	int block = pic_rows * pic_cols;

	#pragma omp parallel for
	for (int i = 0; i < block; i++ ) {
		Extend1[i] = 1.0;  // one
		Extend1[i + block] = list_pic[i];// X
		Extend1[i + block * 2] = list_pic[i + block];//Y
		Extend1[i + block * 3] = list_pic[i + block + block];//Z
		Extend1[i + block * 4] = list_pic[i] * list_pic[i + block];//X_Y
		Extend1[i + block * 5] = list_pic[i + block] * list_pic[i + block + block];//Y_Z
		Extend1[i + block * 6] = list_pic[i]* list_pic[i + block + block];//X_Z
		Extend1[i + block * 7] = pow(list_pic[i], 2); //XYZ_XYZ_MX
		Extend1[i + block * 8] = pow(list_pic[i + block], 2); //XYZ_XYZ_MY
		Extend1[i + block * 9] = pow(list_pic[i + block + block], 2);// XYZ_XYZ_MZ
		Extend1[i + block * 10] = Extend1[i + block * 4] * list_pic[i + block + block]; //X_Y_Z[i];
		Extend1[i + block * 11] = Extend1[i + block * 7] * list_pic[i]; //XYZ_XYZ_XYZ_MX[i];
		Extend1[i + block * 12] = Extend1[i + block * 8] * list_pic[i + block];//XYZ_XYZ_XYZ_MY[i];
		Extend1[i + block * 13] = Extend1[i + block * 9] * list_pic[i + block + block];//XYZ_XYZ_XYZ_MZ
		Extend1[i + block * 14] = Extend1[i + block * 4] * list_pic[i + block];//X_Y_Y[i];
		Extend1[i + block * 15] = Extend1[i + block * 6] * list_pic[i + block + block];//X_Z_Z[i];
		Extend1[i + block * 16] = list_pic[i] * Extend1[i + block * 4];// X_X_Y[i];
		Extend1[i + block * 17] = Extend1[i + block * 5] * list_pic[i + block + block];//Y_Z_Z[i];
		Extend1[i + block * 18] = list_pic[i] * Extend1[i + block * 6]; //X_X_Z[i];
		Extend1[i + block * 19] = list_pic[i + block] * Extend1[i + block * 5]; //Y_Y_Z[i];
	}
	pic_channel = 20;
}

void Create_CorrectXYZ( float* list_C, int row_C, int col_C, float* & extend )  {
	/*----------------------------
	 * 功能 : 做出Create_CorrectXYZ Matlab code如下
	 * Matlab code :
	 * CorrectXYZ = reshape( (C*reshape(extend,[],size(extend,3))')' ,size(pic,1), size(pic,2), 3);
	 *----------------------------
	 */
	

	reshape_3d_to_2dtranspose(extend, pic_rows, pic_cols, pic_channel);
	int m = row_C;
	int n = pic_rows * pic_cols;
	int k = col_C;
	// C 3 x 20
	//extend 20 x 6
	//A = m x k, B = k x n, C = m x n
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, list_C, k, list_pic_2d, n, 0.0, list_mul_result, n);
	pic_channel = 3;
	reshape_2d_to_3dtranspose(pic_rows, pic_cols, pic_channel, extend);

}

 void Extend_Two( float* CorrectXYZ, float* & Extend2) {
	int block = pic_rows * pic_cols;

	#pragma omp parallel for
	for (int i = 0; i < block; i++) {
		Extend2[i] = CorrectXYZ[i];//X[i];
		Extend2[i + block] = CorrectXYZ[i + block];//Y[i];
		Extend2[i + block * 2] = CorrectXYZ[i + block + block];//Z[i];
		Extend2[i + block * 3] = CorrectXYZ[i] * CorrectXYZ[i + block];//X_Y[i];
		Extend2[i + block * 4] = CorrectXYZ[i + block] * CorrectXYZ[i + block + block];//Y_Z[i];
		Extend2[i + block * 5] = CorrectXYZ[i] * CorrectXYZ[i + block + block];//X_Z[i];
		Extend2[i + block * 6] = Extend2[i + block * 3] * CorrectXYZ[i + block + block]; //X_Y_Z[i];
	}
	pic_channel = 7;

}

 void Create_pic_spectrum(float* list_Spectrum, int row_Spectrum, int col_Spectrum, float* Extend2, float* & pic_spectrum) {
	 /*----------------------------
	  * 功能 : 做出pic_spectrum Matlab code如下
	  * Matlab code :
	  * pic_spectrum = reshape(  (EV(spectral_range-379,:)*M*reshape(extend,[],size(extend,3))')' ,size(pic,1), size(pic,2), size(spectral_range,2));
	  *----------------------------
	  */

	 reshape_3d_to_2dtranspose(Extend2, pic_rows, pic_cols, pic_channel);
	 int m = row_Spectrum;
	 int n = pic_rows * pic_cols;
	 int k = col_Spectrum;
	 // list_Spectrum 401 x 7
	 //extend 7 x (row*col)
	 //A = m x k, B = k x n, C = m x n
	 cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, list_Spectrum, k, list_pic_2d, n, 0.0, list_mul_result, n);
	 pic_channel = 401;
	 reshape_2d_to_3dtranspose(pic_rows, pic_cols, pic_channel, pic_spectrum);

 }

 void light_pic_spectrum(float* list_light, int row_Spectrum, int col_Spectrum, float* &pic_spectrum) {
	 //光譜依照光源做調整
 }


void reshape_3d_to_2dtranspose(float* list_3d, int row, int col, int ch ) {
	#pragma omp parallel for
	for (int k = 0; k < ch; k++) {
		for (int j = 0; j < col; j++) {
			for (int i = 0; i < row; i++) {
				list_pic_2d[i + (j * row) + (k * row * col) ] = list_3d[ j + i * col + k * row * col ];
			}
		}
	}

}

void reshape_2d_to_3dtranspose( int row, int col, int ch, float* & list_3d ) {
	#pragma omp parallel for
	for (int k = 0; k < ch; k++) {
		for (int j = 0; j < col; j++) {
			for (int i = 0; i < row; i++) {
				list_3d[j + i * col + k * row * col] = list_mul_result[i + j * row + k * row * col];
			}
		}
	}

}


void print_list_2D( float* list_2d, int row, int col ) {
	//row major
	cout << "list 2d row major" << endl;
	for (int i = 0; i < row * col; i++ ) {
		cout << list_2d[i]<< " ";
		if ((i + 1) % col == 0)
			cout << endl;
	}
}

void print_list_3D(float* list_3d, int row, int col, int ch ) {
	
	cout << "list 3d row major" << endl;
	for ( int i = 0; i < row * col * ch; i++ ) {
		cout << i << ": " << list_3d[i] << " ";
		if ((i + 1) % col == 0)
			cout << endl;
		if ((i + 1) % (col * row) == 0)
			cout << endl;
	}
	
}



 void Set_img_spectrum( Mat table, int filenumber, int RGB_R) {

	//---------------讀取參數---------------
	//load C.txt    3*20 row*col
	int row_C = 3;
	int col_C = 20;
	float* list_C = (float*)malloc(row_C * col_C * sizeof(float));
	LoadData("C.txt", list_C, row_C, col_C, 0);

	
	//load Spectrum.txt 401*7
	int row_Spectrum = 401;
	int col_Spectrum = 7;
	float* list_Spectrum = (float*)malloc(row_Spectrum * col_Spectrum * sizeof(float));
	LoadData("Spectrum.txt", list_Spectrum, row_Spectrum, col_Spectrum, 0);
	
	//load light.txt    401*1
	int row_light = 401;
	int col_light = 1;
	float* list_light = (float*)malloc(row_light * col_light * sizeof(float));
	LoadData("light.txt", list_light, row_light, col_light, 0);


	
	//load CMF.txt  401*3
	int row_CMF = 401;
	int col_CMF = 3;
	float* list_CMF = (float*)malloc(row_CMF * col_CMF * sizeof(float));
	LoadData("CMF.txt", list_CMF, row_CMF, col_CMF, 0);

	
	//load XYZ.txt  3*3
	int row_XYZ = 3;
	int col_XYZ = 3;
	float* list_XYZ = (float*)malloc(row_XYZ * col_XYZ * sizeof(float));
	LoadData("XYZ.txt", list_XYZ, row_XYZ, col_XYZ, 0);


	

	//---------------讀取參數---------------


	


	//計算部分開始"開始計時"
	
	double START, END;
	START = clock();

	//opencv2XYZ2list
	float* list_pic = (float*)malloc(pic_rows * pic_cols * 3 * sizeof(float));
	Matlab_RGB2XYZ(table, list_pic);
	//cout << "end Matlab_RGB2XYZ" << endl;


	Create_XYZ(list_XYZ, row_XYZ, col_XYZ, list_pic);
	//cout << "end Create_XYZ" << endl;

	float* Extend1 = (float*)malloc(pic_rows * pic_cols * 20 * sizeof(float));
	Extend_One(list_pic, Extend1);
	//cout << "end Extend1" << endl;


	Create_CorrectXYZ( list_C, row_C, col_C, Extend1 );
	//cout << "end Create_CorrectXYZ" << endl;
	
	float* Extend2 = (float*)malloc(pic_rows * pic_cols * 7 * sizeof(float));
	Extend_Two(Extend1, Extend2);
	//cout << "end Extend_Two" << endl;


	float* pic_spectrum = (float*)malloc(pic_rows * pic_cols * 401 * sizeof(float));
	Create_pic_spectrum(list_Spectrum, row_Spectrum, col_Spectrum, Extend2, pic_spectrum);
	//cout << "end Create_pic_spectrum" << endl;
	END = clock();
	//cout << endl << "進行運算所花費的時間：" << (END - START) / CLOCKS_PER_SEC << " S" << endl;
	
	StoreData(pic_spectrum, pic_rows, pic_cols, 401, table, filenumber);
	cout << filenumber << "nm: " << "R = " << RGB_R << " ";
	cout << "StoreData所花費的時間：" << (clock() - START) / CLOCKS_PER_SEC << " S" << endl;

}


 void printf_table( Mat table ) {
	 for (int i = 0; i < pic_rows; i++) {
		 for (int j = 0; j < pic_cols; j++) {
			 for (int k = 0; k < 3; k++) {
				 cout << table.at<cv::Vec3i>(i, j)[k] << " ";
			 }
			 cout << endl;
		 }
	 }
 }


char* TestFun( string str ){
    const char* result = str.c_str();
	int tail = 0;
	while (result[tail] != '\0') {
		tail++;
	}
	char* real_result = (char*)malloc(tail+1);
	for (int i = 0; i <= tail; i++) {
		real_result[i] = result[i];
	}
		return real_result;
}



void firstopenfile( int filenumber ) {
	 FILE* fp;
	 string list_msg = int2str(filenumber) + "nm";
	 string filename = int2str(filenumber) + "nm.csv";
	 fp = fopen(filename.c_str(), "w");

	 //第一列資訊
	 fprintf(fp, "id");
	 fprintf(fp, ",");
	 fprintf(fp, "R");
	 fprintf(fp, ",");
	 fprintf(fp, "G");
	 fprintf(fp, ",");
	 fprintf(fp, "B");
	 fprintf(fp, ",");
	 fprintf(fp, "\n");
	 fclose(fp);
 }

int main() {
	omp_set_num_threads(8);


	#define BIG (256LL * 256 * 401 )
	list_pic_2d = (float*)malloc(BIG * sizeof(float));
	list_mul_result = (float*)malloc(BIG * sizeof(float));


	int OneFile = 256 * 256;
	Mat table = cv::Mat(1, OneFile, CV_32SC3, cv::Scalar(0, 0, 0));
	int R = 0;
	int G = 0;
	int B = 0;

	int times = 256;

	for (int spectrum = 380; spectrum <= 780; spectrum++) {
		//先設定第一次file的資訊
		firstopenfile( spectrum );

		//初始化data_id, R,G,B
		data_id = 0;
		R = 0;
		G = 0;
		B = 0;
		for (int k = 0; k < times; k++) {
			pic_rows = 1;
			pic_cols = OneFile;
			pic_channel = 3;

			for (int i = 0; i < pic_rows; i++) {
				for (int j = 0; j < pic_cols; j++) {
					table.at<cv::Vec3i>(i, j)[0] = R;
					table.at<cv::Vec3i>(i, j)[1] = G;
					table.at<cv::Vec3i>(i, j)[2] = B;
					B++;
					if (B > 255) {
						B = 0;
						G++;
						if (G > 255) {
							G = 0;
							R++;
						}
					}
				}
			}

			Set_img_spectrum(table, spectrum, k);
		}
	}

}
