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

int pic_rows = 0;
int pic_cols = 0;
int pic_channel = 0;
int pic_element = 0;

void print_vector_2D(vector<vector<float> > );
void print_vector_3D(vector<vector< vector<float> > > );
void reshape_3d_to_2d(vector<vector< vector<float> > >, vector<vector<float> >&);
void reshape_2d_to_3d(vector<vector<float> >, vector<vector< vector<float> > > & );
void matrix_mul_transpose(vector<vector<float> > , vector<vector<float> > , vector<vector<float> >&);
void mul_count_3d_matrix(vector<vector< vector<float> > >& , float );


/*----------------------------
 * 功能 : 從 .txt 檔案中讀入資料，儲存到 cv::Mat 矩陣
 *      - 預設按 float 格式讀入資料，
 *      - 如果沒有指定矩陣的行、列和通道數，則輸出的矩陣是單通道、N 行 1 列的
 *----------------------------
 * 函式 : LoadData
 * 訪問 : public
 * 返回 : -1：開啟檔案失敗；0：按設定的矩陣引數讀取資料成功；1：按預設的矩陣引數讀取資料
 *
 * 引數 : fileName    [in]    檔名
 * 引數 : matData [out]   矩陣資料
 * 引數 : matRows [in]    矩陣行數，預設為 0
 * 引數 : matCols [in]    矩陣列數，預設為 0
 * 引數 : matChns [in]    矩陣通道數，預設為 0
 */
inline int LoadData(string fileName, vector<vector<float> > & vecData, int Rows = 0, int Cols = 0, int Chns = 0)
{
	ifstream fin( fileName.c_str() );
	if (!fin) {
		cout << "無法讀入檔案\n";
		return 0;
	}
	for (int i = 0; i < Rows; i++)
		for (int j = 0; j < Cols; j++) fin >> vecData[i][j];
	fin.close();
	return 0;
}

void Matshape(cv::Mat  img, string text) {
	//img.cols 宽
	//img.rows 高
	cout << text << "(高, 寬, 通道) = ( " << img.rows << ", " << img.cols << ", " << img.channels() << " )" << endl;
}

inline cv::Mat Matc3(cv::Mat mat)
{
	cv::Mat M0(mat.rows, mat.cols, CV_32F);
	cv::Mat M1(mat.rows, mat.cols, CV_32F);
	cv::Mat M2(mat.rows, mat.cols, CV_32F);
	#pragma omp parallel for
	for (int i = 0; i < mat.rows; i++)
	{
		#pragma omp parallel for
		for (int j = 0; j < mat.cols; j++)
		{
			M0.at<float>(i, j) = mat.at<cv::Vec3f>(i, j)[0];
			M1.at<float>(i, j) = mat.at<cv::Vec3f>(i, j)[1];
			M2.at<float>(i, j) = mat.at<cv::Vec3f>(i, j)[2];
		}
	}
	cv::vconcat(M0, M1, M1);
	cv::vconcat(M1, M2, M2);
	return M2;
}

void StoreData(cv::Mat pic_spectrum ) {
	typedef cv::Vec<float, 401> Vec401f;

	FILE* fp;
	fp = fopen("Output_data.txt", "w");

	for (int k = 0; k < pic_spectrum.channels(); k++) {
		fprintf(fp, "val(:,:,%d) =\n\n", k + 1);
		for (int i = 0; i < pic_spectrum.rows; i++) {
			for (int j = 0; j < pic_spectrum.cols; j++) {

				fprintf(fp, "%10.4f", pic_spectrum.at<Vec401f>(i, j)[k]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void Matlab_RGB2XYZ(cv::Mat im, vector<vector< vector<float> > > & vec_pic ) {
	/*----------------------------
	* 功能 : 做出Matlab_RGB2XYZ如下
	* Matlab code :
	* XYZ_D65 = rgb2xyz(pic);
	*----------------------------
	*/
	float XYZ_M_float[3][3] = { 0.4124564, 0.3575761, 0.1804375, 0.2126729, 0.7151522, 0.0721750, 0.0193339,  0.1191920,  0.9503041 };
	vector<vector<float>> XYZ_M = { { 0.4124564, 0.3575761, 0.1804375 }, { 0.2126729, 0.7151522, 0.0721750 }, { 0.0193339,  0.1191920,  0.9503041 } };
	im.convertTo(im, CV_32F);
	im = im / 255;
	for (int i = 0; i < im.rows; i++)
	{	
		for (int j = 0; j < im.cols; j++)
		{
			//RGB2XYZ(vector)
			for (int k = 0; k < im.channels(); k++) {
				if( im.at<cv::Vec3f>(i, j)[k] < 0.04045 )
					vec_pic[i][j][k] = im.at<cv::Vec3f>(i, j)[k] / 12.92;
				else
					vec_pic[i][j][k] = pow(((im.at<cv::Vec3f>(i, j)[k] + 0.055) / 1.055), 2.4);
			}
		}
	}
	print_vector_3D(vec_pic);
	system("pause");

	
	//XYZ_D65是3維 Mat_XYZ_M是2維(3x3)
	//XYZ_D65必須與Mat_XYZ_M矩陣相乘
	//XYZ_D65先.reshape(1, XYZ_D65.cols * XYZ_D65.rows); 再將reshape後的XYZ_D65.t();
	//做完Mat_XYZ_M * XYZ_D65 再 (Mat_XYZ_M * XYZ_D65).t(); 
	vector<vector<float> > vec_pic_2d(vec_pic.size()* vec_pic[0].size(), vector<float>(vec_pic[0][0].size()));
	vector<vector<float> > vec_2dmul(vec_pic_2d.size(), vector<float>(XYZ_M.size(), 0));
	reshape_3d_to_2d(vec_pic, vec_pic_2d);
	matrix_mul_transpose(XYZ_M, vec_pic_2d, vec_2dmul);
	print_vector_2D(vec_2dmul);
	reshape_2d_to_3d(vec_2dmul, vec_pic);
	cout << endl;
	print_vector_3D(vec_pic);
	system("pause");
	mul_count_3d_matrix(vec_pic, 100.0);
	print_vector_3D(vec_pic);
	system("pause");

	
}

 inline void Create_XYZ(cv::Mat Mat_Ma, cv::Mat Mat_White_light, cv::Mat Mat_White_D65, cv::Mat & XYZ_D65 ) {
	/*----------------------------
	 * 功能 : 做出XYZ Matlab code如下
	 * Matlab code :
	 * XYZ = reshape(  (inv(Ma)*diag( (Ma*White_light)./(Ma*White_D65) )*Ma*reshape(XYZ_D65,[],3)')' ,size(pic,1), size(pic,2), 3);
	 *----------------------------
	 */
	cv::Mat diag_M;

	diag_M = ((Mat_Ma * Mat_White_light) / (Mat_Ma * Mat_White_D65));
	diag_M = diag_M.diag(diag_M); // diag_M 為(Mat_Ma * Mat_White_light) / (Mat_Ma * Mat_White_D65)"對角方陣"

	XYZ_D65 = XYZ_D65.reshape(1, XYZ_D65.cols * XYZ_D65.rows);
	cv::transpose(XYZ_D65, XYZ_D65);

	// inv(Ma)*diag( (Ma*White_light)./(Ma*White_D65) )*Ma*reshape(XYZ_D65,[],3)'
	XYZ_D65 = Mat_Ma.inv() * diag_M * Mat_Ma * XYZ_D65;
	cv::transpose(XYZ_D65, XYZ_D65);
	XYZ_D65 = XYZ_D65.reshape(pic_channel, pic_rows);

}

 inline void Extend_One(cv::Mat XYZ, cv::Mat & Extend1) {
	cv::Mat MX(XYZ.rows, XYZ.cols, CV_32F);//X 2維
	cv::Mat MY(XYZ.rows, XYZ.cols, CV_32F);//Y 2維
	cv::Mat MZ(XYZ.rows, XYZ.cols, CV_32F);//Z 2維
	#pragma omp parallel for
	for (int i = 0; i < XYZ.rows; i++)
	{
		#pragma omp parallel for
		for (int j = 0; j < XYZ.cols; j++)
		{
			MX.at<float>(i, j) = XYZ.at<cv::Vec3f>(i, j)[0];
			MY.at<float>(i, j) = XYZ.at<cv::Vec3f>(i, j)[1];
			MZ.at<float>(i, j) = XYZ.at<cv::Vec3f>(i, j)[2];
		}
	}

	cv::Mat ones = cv::Mat::ones(XYZ.rows, XYZ.cols, CV_32FC1); //ones
	cv::Mat X_Y, Y_Z, X_Z, XYZ_XYZ, X_Y_Z; //一階&&二階
	cv::Mat XYZ_XYZ_XYZ, X_Y_Y, X_Z_Z, X_X_Y, Y_Z_Z, X_X_Z, Y_Y_Z; //三階
	cv::Mat M_merge(XYZ.rows, XYZ.cols, CV_32F);

	X_Y = MX.mul(MY);
	Y_Z = MY.mul(MZ);
	X_Z = MX.mul(MZ);
	XYZ_XYZ = XYZ.mul(XYZ);
	X_Y_Z = X_Y.mul(MZ);
	XYZ_XYZ_XYZ = XYZ_XYZ.mul(XYZ);
	X_Y_Y = X_Y.mul(MY);
	X_Z_Z = X_Z.mul(MZ);
	X_X_Y = MX.mul(X_Y);
	Y_Z_Z = Y_Z.mul(MZ);
	X_X_Z = MX.mul(X_Z);
	Y_Y_Z = MY.mul(Y_Z);

	cv::vconcat(ones, MX, M_merge);
	cv::vconcat(M_merge, MY, M_merge);
	cv::vconcat(M_merge, MZ, M_merge);
	cv::vconcat(M_merge, X_Y, M_merge);
	cv::vconcat(M_merge, Y_Z, M_merge);
	cv::vconcat(M_merge, X_Z, M_merge);
	cv::vconcat(M_merge, Matc3(XYZ_XYZ), M_merge);
	cv::vconcat(M_merge, X_Y_Z, M_merge);
	cv::vconcat(M_merge, Matc3(XYZ_XYZ_XYZ), M_merge);
	cv::vconcat(M_merge, X_Y_Y, M_merge);
	cv::vconcat(M_merge, X_Z_Z, M_merge);
	cv::vconcat(M_merge, X_X_Y, M_merge);
	cv::vconcat(M_merge, Y_Z_Z, M_merge);
	cv::vconcat(M_merge, X_X_Z, M_merge);
	cv::vconcat(M_merge, Y_Y_Z, M_merge);

	typedef cv::Vec<float, 20> Vec20f;
	int divisor = 0, n = 0, n2 = 0, math = 0, buffer = 0;

	for (int i = 0; i < M_merge.rows; i++)
	{
		if (XYZ.rows % 2 != 0)
		{
			if (i >= XYZ.rows)
			{
				divisor = i - (XYZ.rows) * math;
			}
			else
			{
				divisor = i;
			}
		}
		else
		{
			divisor = i;
		}

		for (int j = 0; j < M_merge.cols; j++)
		{
			if (divisor % 2 == 0)
			{
				Extend1.at<Vec20f>(i % XYZ.rows, j)[math] = M_merge.at<float>(i, j);
				n2++;
			}
			else if (divisor % 2 != 0)
			{
				Extend1.at<Vec20f>(i % XYZ.rows, j)[math] = M_merge.at<float>(i, j);
				n++;
			}
			if (XYZ.rows % 2 != 0)
			{
				if (n2 >= ((XYZ.rows / 2) + 1) * M_merge.cols)
				{
					math++;
					if (math >= 20)
					{
						math--;
					}
					n2 = 0;
				}
			}
			else
			{
				if (n >= (XYZ.rows / 2) * M_merge.cols)
				{
					math++;
					if (math >= 20)
					{
						math--;
					}
					n = 0;
				}
			}
		}
		buffer++;
	}
}

inline void Create_CorrectXYZ(cv::Mat & extend, cv::Mat Mat_C ) {
	/*----------------------------
	 * 功能 : 做出Create_CorrectXYZ Matlab code如下
	 * Matlab code :
	 * CorrectXYZ = reshape( (C*reshape(extend,[],size(extend,3))')' ,size(pic,1), size(pic,2), 3);
	 *----------------------------
	 */

	extend = extend.reshape(1, extend.rows * extend.cols);
	cv::transpose(extend, extend);
	extend = Mat_C * extend;
	cv::transpose(extend, extend);
	extend = extend.reshape(pic_channel, pic_rows);

}

 inline void Extend_Two( cv::Mat CorrectXYZ, cv::Mat & Extend2) {
	cv::Mat MX(CorrectXYZ.rows, CorrectXYZ.cols, CV_32F);//X 2維
	cv::Mat MY(CorrectXYZ.rows, CorrectXYZ.cols, CV_32F);//Y 2維
	cv::Mat MZ(CorrectXYZ.rows, CorrectXYZ.cols, CV_32F);//Z 2維

	cv::Mat X_Y, Y_Z, X_Z, X_Y_Z ;
	cv::Mat M_merge(CorrectXYZ.rows, CorrectXYZ.cols, CV_32F);
	#pragma omp parallel for
	for (int i = 0; i < CorrectXYZ.rows; i++)
	{
		#pragma omp parallel for
		for (int j = 0; j < CorrectXYZ.cols; j++)
		{
			MX.at<float>(i, j) = CorrectXYZ.at<cv::Vec3f>(i, j)[0];
			MY.at<float>(i, j) = CorrectXYZ.at<cv::Vec3f>(i, j)[1];
			MZ.at<float>(i, j) = CorrectXYZ.at<cv::Vec3f>(i, j)[2];
		}
	}
	X_Y = MX.mul(MY); //XY
	Y_Z = MY.mul(MZ);//YZ
	X_Z = MX.mul(MZ);//XZ
	X_Y_Z = X_Y.mul(MZ);//XYZ

	cv::vconcat(MX, MY, M_merge);
	cv::vconcat(M_merge, MZ, M_merge);
	cv::vconcat(M_merge, X_Y, M_merge);
	cv::vconcat(M_merge, Y_Z, M_merge);
	cv::vconcat(M_merge, X_Z, M_merge);
	cv::vconcat(M_merge, X_Y_Z, M_merge);

	int n = 0, math = 0;
	typedef cv::Vec<float, 7> Vec7f;

	for (int i = 0; i < M_merge.rows; i++)
	{

		for (int j = 0; j < M_merge.cols; j++)
		{
			if (i % 2 == 0)
			{
				Extend2.at<Vec7f>(i % CorrectXYZ.rows, j)[math] = M_merge.at<float>(i, j);
			}
			else if (i % 2 != 0)
			{
				Extend2.at<Vec7f>(i % CorrectXYZ.rows, j)[math] = M_merge.at<float>(i, j);
				n ++;
			}
			if (n >= (CorrectXYZ.rows / 2) * M_merge.cols)
			{
				math ++;
				if (math >= 7)
				{
					math--;
				}
				n = 0;
			}
		}
	}

}

 inline void Create_pic_spectrum(cv::Mat & Extend2, cv::Mat Mat_EV, cv::Mat Mat_M ) {
	/*----------------------------
	 * 功能 : 做出pic_spectrum Matlab code如下
	 * Matlab code :
	 * pic_spectrum = reshape(  (EV(spectral_range-379,:)*M*reshape(extend,[],size(extend,3))')' ,size(pic,1), size(pic,2), size(spectral_range,2));
	 *----------------------------
	 */

	Extend2 = Extend2.reshape(1, Extend2.rows * Extend2.cols);
	cv::transpose(Extend2, Extend2);
	Extend2 = Mat_EV * Mat_M * Extend2;
	cv::transpose(Extend2, Extend2);
	Extend2 = Extend2.reshape(401, pic_rows);

}

 void print_vector_2D(vector<vector<float> > vec_C ) {
	 for (int i = 0; i < vec_C.size(); i++) {
		 for (int j = 0; j < vec_C[i].size(); j++)
			 cout << vec_C[i][j] << " ";
		 cout << endl;
	 }
 }

void print_vector_3D(vector<vector< vector<float> > > vec_C) {

	for (int k = 0; k < vec_C[0][0].size(); k++) {
		for (int i = 0; i < vec_C.size(); i++) {
			for (int j = 0; j < vec_C[i].size(); j++) {
				cout << vec_C[i][j][k] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}

 }

void reshape_3d_to_2d(vector<vector< vector<float> > > vec_3d, vector<vector<float> >& vec_2d) {

	if (vec_3d.empty()) {
		cout << "reshape_3d_to_2d : vec is empty " << endl;
	}
	else {
		for (int k = 0; k < vec_3d[0][0].size(); k++) {
			for (int j = 0; j < vec_3d[0].size(); j++) {
				for (int i = 0; i < vec_3d.size(); i++) {
					vec_2d[i + j * vec_3d.size()][k] = vec_3d[i][j][k];
				}
			}
		}
	}

}

void reshape_2d_to_3d(vector<vector<float> > vec_2d, vector<vector< vector<float> > > & vec_3d) {
	if (vec_3d.empty()) {
		cout << "reshape_2d_to_3d : vec is empty " << endl;
	}
	else {
		for (int k = 0; k < vec_3d[0][0].size(); k++) {
			for (int j = 0; j < vec_3d[0].size(); j++) {
				for (int i = 0; i < vec_3d.size(); i++) {
					vec_3d[i][j][k] = vec_2d[i + j * vec_3d.size()][k];
				}
			}
		}
	}
}

void matrix_mul_transpose(vector<vector<float> > vec_arrA, vector<vector<float> > vec_arrB, vector<vector<float> >& vec_arrC) {
	//必須乘法的前後包含transpose才可以使用此function
	//非正常matrix_mul
	int rowA = vec_arrA.size();
	int colA = vec_arrA[0].size();
	int rowB = vec_arrB.size();
	int colB = vec_arrB[0].size();

	if (colA != colB) cout << "can't use matrix_mul_transpose" << endl;
	//colA == colB
	else {
		for (int i = 0; i < rowA; i++) {
			for (int j = 0; j < rowB; j++) {
				for (int k = 0; k < colA; k++) {
					vec_arrC[j][i] += vec_arrA[i][k] * vec_arrB[j][k];

				}
			}
		}
	}
}
void mul_count_3d_matrix(vector<vector< vector<float> > >& vec_3d, float num ) {
	for (int k = 0; k < vec_3d[0][0].size(); k++) {
		for (int i = 0; i < vec_3d.size(); i++) {
			for (int j = 0; j < vec_3d[i].size(); j++) {
				vec_3d[i][j][k] = vec_3d[i][j][k] * num;
			}
		}
	}
}

 inline void Set_img_spectrum( string data_path ) {
	clock_t start, end;

	//---------------讀取參數---------------
	//load C.txt    3*20 row*col
	vector<vector<float> > vec_C (3, vector<float>(20, 0));
	LoadData("C.txt", vec_C, 3, 20, 0);


	//load M.txt    12*7
	vector<vector<float> > vec_M(12, vector<float>(7, 0));
	LoadData("M.txt", vec_M, 12, 7, 0);


	//load EV.txt   401*12
	vector<vector<float> > vec_EV(401, vector<float>(12, 0));
	LoadData("EV.txt", vec_EV, 401, 12, 0);


	//load light.txt    401*1
	vector<vector<float> > vec_light(401, vector<float>(1, 0));
	LoadData("light.txt", vec_light, 401, 1, 0);



	//load Camera white point.txt   3*1
	vector<vector<float> > vec_White_D65(3, vector<float>(1, 0));
	LoadData("Camera_white_point.txt", vec_White_D65, 3, 1, 0);


	//load light white point.txt    3*1
	vector<vector<float> > vec_White_light(3, vector<float>(1, 0));
	LoadData("light_white_point.txt", vec_White_light, 3, 1, 0);


	//load CMF.txt  401*3
	vector<vector<float> > vec_CMF(401, vector<float>(3, 0));
	LoadData("CMF.txt", vec_CMF, 401, 3, 0);


	//set Ma 3*3
	vector<vector<float> > vec_Ma = { { 0.40024, 0.70760, -0.08081 }, { -0.22603, 1.16532, 0.04570 }, { 0,  0,  0.91822 } };


	//---------------讀取參數---------------



	//---------------讀取影像---------------
	
	//讀入影象，並將之轉為3通道影象
	cv::Mat im = cv::imread( data_path, 3 );
	pic_rows = im.rows;
	pic_cols = im.cols;
	pic_channel = im.channels();

	cout << "圖片讀取成功!!" << endl;
	Matshape(im, "image");
	cvtColor(im, im, cv::COLOR_BGR2RGB); // opencv原本是bgr表示轉成rgb
	StoreData(im);

	vector<vector< vector<float> > > vec_pic ;
	//初始化3D矩陣
	vec_pic.resize(pic_rows);
	for (int i = 0; i < pic_rows; i++) {
		vec_pic[i].resize(pic_cols);
		for (int j = 0; j < pic_cols; j++) 
			vec_pic[i][j].resize(pic_channel);
	}
	cout <<  "pic_rows: " << vec_pic.size() << endl;
	cout << "pic_cols: " << vec_pic[0].size() << endl;
	cout << "pic_channel: " << vec_pic[0][0].size() << endl;

	//---------------讀取影像---------------
	
	//計算部分開始"開始計時"
	start = clock();

	//opencv2XYZ2vector
	Matlab_RGB2XYZ(im, vec_pic );
	/*

	
	Create_XYZ(Mat_Ma, Mat_White_light, Mat_White_D65, im);
	Mat_Ma.release();
	Mat_White_light.release();
	Mat_White_D65.release();

	
	cv::Mat Extend1 = cv::Mat::zeros(pic_rows, pic_cols, CV_32FC(20));
	Extend_One(im, Extend1);
	im.release();

	
	Create_CorrectXYZ(Extend1, Mat_C );
	Mat_C.release();


	
	cv::Mat Extend2 = cv::Mat::zeros(pic_rows, pic_cols, CV_32FC(7));
	Extend_Two(Extend1, Extend2);
	Extend1.release();//CorrectXYZ釋放


	
	cv::Mat pic_spectrum = cv::Mat::zeros(pic_rows, pic_cols, CV_32FC(401));
	Create_pic_spectrum( Extend2, Mat_EV, Mat_M );
	Mat_EV.release();
	Mat_M.release();
	
	Matshape(pic_spectrum, "***final output***  pic_spectrum ");
	end = clock();
	cout << "time consume: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;

	return Extend2;
	*/
	
}

int main() {

	omp_set_num_threads(4);
	cv::Mat pic_spectrum = cv::Mat::zeros(pic_rows, pic_cols, CV_32FC(401)); //img_spectrum
	string data_path = "..\\small.jpg";

	Set_img_spectrum(data_path);


	cout << "enter to start store data : pic_spectrum" << endl;
	system("pause");
	StoreData(pic_spectrum); //寫檔file name is "Output_data.txt"
	cout << "Finsh store data !!! " << endl;
}
