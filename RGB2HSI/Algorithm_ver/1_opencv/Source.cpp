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
inline int LoadData(string fileName, cv::Mat& matData, int matRows = 0, int matCols = 0, int matChns = 0)
{
	int retVal = 0;
	// 開啟檔案
	ifstream inFile(fileName.c_str(), ios_base::in);
	if (!inFile.is_open())
	{
		cout << "讀取檔案失敗" << endl;
		retVal = -1;
		return (retVal);
	}
	// 載入資料
	istream_iterator<float> begin(inFile);    //按 float 格式取檔案資料流的起始指標
	istream_iterator<float> end;          //取檔案流的終止位置
	vector<float> inData(begin, end);      //將檔案資料儲存至 std::vector 中
	cv::Mat tmpMat = cv::Mat(inData);       //將資料由 std::vector 轉換為 cv::Mat
	// 輸出到命令列視窗
	//copy(vec.begin(),vec.end(),ostream_iterator<double>(cout,"\t")); 
	// 檢查設定的矩陣尺寸和通道數
	size_t dataLength = inData.size();
	//1.通道數
	if (matChns == 0)
	{
		matChns = 1;
	}
	//2.行列數
	if (matRows != 0 && matCols == 0)
	{
		matCols = dataLength / matChns / matRows;
	}
	else if (matCols != 0 && matRows == 0)
	{
		matRows = dataLength / matChns / matCols;
	}
	else if (matCols == 0 && matRows == 0)
	{
		matRows = dataLength / matChns;
		matCols = 1;
	}
	//3.資料總長度
	if (dataLength != (matRows * matCols * matChns))
	{
		cout << "讀入的資料長度 不滿足 設定的矩陣尺寸與通道數要求，將按預設方式輸出矩陣！" << endl;
		retVal = 1;
		matChns = 1;
		matRows = dataLength;
	}
	// 將檔案資料儲存至輸出矩陣
	matData = tmpMat.reshape(matChns, matRows).clone();
	return (retVal);
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

inline void Matlab_RGB2XYZ(cv::Mat & im ) {
	/*----------------------------
	* 功能 : 做出Matlab_RGB2XYZ如下
	* Matlab code :
	* XYZ_D65 = rgb2xyz(pic);
	*----------------------------
	*/
	float XYZ_M[3][3] = { 0.4124564, 0.3575761, 0.1804375, 0.2126729, 0.7151522, 0.0721750, 0.0193339,  0.1191920,  0.9503041 };
	cv::Mat Mat_XYZ_M = cv::Mat(3, 3, CV_32F, XYZ_M);
	float r = 0;
	float g = 0;
	float b = 0;
	im.convertTo(im, CV_32F);
	im = im / 255;
	#pragma omp parallel for
	for (int i = 0; i < im.rows; i++)
	{	
		#pragma omp parallel for
		for (int j = 0; j < im.cols; j++)
		{
			// r = im.at<cv::Vec3f>(i, j)[0]
			// g = im.at<cv::Vec3f>(i, j)[1]
			// b = im.at<cv::Vec3f>(i, j)[2]

			//r
			if (im.at<cv::Vec3f>(i, j)[0] < 0.04045)
				im.at<cv::Vec3f>(i, j)[0] = im.at<cv::Vec3f>(i, j)[0] / 12.92;
			else
				im.at<cv::Vec3f>(i, j)[0] = pow(((im.at<cv::Vec3f>(i, j)[0] + 0.055) / 1.055), 2.4);
			//g
			if (im.at<cv::Vec3f>(i, j)[1] < 0.04045)
				im.at<cv::Vec3f>(i, j)[1] = im.at<cv::Vec3f>(i, j)[1] / 12.92;
			else
				im.at<cv::Vec3f>(i, j)[1] = pow(((im.at<cv::Vec3f>(i, j)[1] + 0.055) / 1.055), 2.4);
			//b
			if (im.at<cv::Vec3f>(i, j)[2] < 0.04045)
				im.at<cv::Vec3f>(i, j)[2] = im.at<cv::Vec3f>(i, j)[2] / 12.92;
			else
				im.at<cv::Vec3f>(i, j)[2] = pow(((im.at<cv::Vec3f>(i, j)[2] + 0.055) / 1.055), 2.4);

		}
	}

	//XYZ_D65是3維 Mat_XYZ_M是2維(3x3)
	//XYZ_D65必須與Mat_XYZ_M矩陣相乘
	//XYZ_D65先.reshape(1, XYZ_D65.cols * XYZ_D65.rows); 再將reshape後的XYZ_D65.t();
	//做完Mat_XYZ_M * XYZ_D65 再 (Mat_XYZ_M * XYZ_D65).t(); 
	cv::Mat reshape_M;
	reshape_M = im.reshape(1, im.cols * im.rows);
	cv::transpose(reshape_M, reshape_M);
	reshape_M = Mat_XYZ_M * reshape_M;
	cv::transpose(reshape_M, reshape_M);
	im = reshape_M.reshape(3, im.rows);
	im = im * 100;
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

 inline cv::Mat Set_img_spectrum( string data_path ) {
	clock_t start, end;

	//---------------讀取參數---------------
	//load C.txt    3*20 row*col
	cv::Mat Mat_C(3, 20, CV_32F);
	LoadData("C.txt", Mat_C, 3, 20, 0);
	//cout << Mat_C << endl;


	//load M.txt    12*7
	cv::Mat Mat_M(12, 7, CV_32F);
	LoadData("M.txt", Mat_M, 12, 7, 0);
	//cout << Mat_M << endl;


	//load EV.txt   401*12
	cv::Mat Mat_EV(401, 12, CV_32F);
	LoadData("EV.txt", Mat_EV, 401, 12, 0);
	//cout << Mat_EV << endl;


	//load light.txt    401*1
	cv::Mat Mat_light(401, 1, CV_32F);
	LoadData("light.txt", Mat_light, 401, 1, 0);
	//cout << Mat_light << endl;


	//load Camera white point.txt   3*1
	cv::Mat Mat_White_D65(3, 1, CV_32F);
	LoadData("Camera_white_point.txt", Mat_White_D65, 3, 1, 0);
	//cout << Mat_White_D65 << endl;


	//load light white point.txt    3*1
	cv::Mat Mat_White_light(3, 1, CV_32F);
	LoadData("light_white_point.txt", Mat_White_light, 3, 1, 0);
	//cout << Mat_White_light << endl;


	//load CMF.txt  401*3
	cv::Mat Mat_CMF(401, 3, CV_32F);
	LoadData("CMF.txt", Mat_CMF, 401, 3, 0);
	//cout << Mat_CMF << endl;


	//set Ma
	float Ma[3][3] = { 0.40024, 0.70760, -0.08081, -0.22603, 1.16532, 0.04570, 0,  0,  0.91822 };
	cv::Mat Mat_Ma = cv::Mat(3, 3, CV_32F, Ma);
	//cout << Mat_Ma << endl;
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

	//---------------讀取影像---------------

	//計算部分開始"開始計時"
	start = clock();


	//---------------影像轉換COLOR_RGB2XYZ---------------
	Matlab_RGB2XYZ(im);
	//---------------影像轉換COLOR_RGB2XYZ---------------

	
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
	
}

int main() {
	omp_set_num_threads(4);
	cv::Mat pic_spectrum = cv::Mat::zeros(pic_rows, pic_cols, CV_32FC(401)); //img_spectrum
	string data_path = "..\\1920x1080.jpg";

	pic_spectrum = Set_img_spectrum(data_path);

	cout << "enter to start store data : pic_spectrum" << endl;
	system("pause");
	StoreData(pic_spectrum); //寫檔file name is "Output_data.txt"
	cout << "Finsh store data !!! " << endl;
}
