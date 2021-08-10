#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <string>
#include <opencv2/opencv.hpp>
#include <time.h>
#include <math.h>
#include <iterator>
using namespace std;
using namespace cv;

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
int LoadData(string fileName, cv::Mat& matData, int matRows = 0, int matCols = 0, int matChns = 0)
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

int StoreData( string fileName, cv::Mat& matData, int matRows = 0, int matCols = 0 ){
	fstream file;
	file.open(fileName.c_str(), ios::out);

	for( int i = 0; i < matRows; i++ ){
		for( int j = 0; j < matCols; j++ ){
			file << matData.at<float>(i,j) << '\t';
		}
		file << endl;
	}
	file.close();
}
  
int main(){
	string fileName = "out.txt";
	cv::Mat A(2, 3, CV_32F);
	LoadData("a.txt", A, 2, 3, 0);
	cout << A << endl;

	cv::Mat B(3, 2, CV_32F);
	LoadData("b.txt", B, 3, 2, 0);
	cout << B << endl;
	
	cv::Mat C(2, 2, CV_32F);
	C = A*B;
	cout << C << endl;
	StoreData( fileName, C, 2, 2 );

}
