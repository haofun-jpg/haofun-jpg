#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include<vector>
#include "math.h"
#include <string> 
#include <fstream>
#include <sstream>
using namespace std;



double calcRMSE(vector<double> Data,vector<double> Data2, int Num)
{

	double fSum = 0;
	for (int i = 0; i < Num; ++i)
	{
		fSum += (Data[i] - Data2[i]) *(Data[i] - Data2[i]);
	}
	return sqrt(fSum / Num);

}


void load_csv( vector<int> &matrix_id, vector<double> &matrix_color, vector<double> &matrix_algorithm, vector<double> &matrix_regression, string filename ){
    ifstream inFile( filename.c_str(), ios::in);
    if (!inFile)
    {
        cout << "開啟檔案失敗！" << endl;
        exit(1);
    }
    int i = 0;
    string line;
    string field;
    getline(inFile, line);
    getline(inFile, line);
    getline(inFile, line);
    getline(inFile, line);
    while (getline(inFile, line))//getline(inFile, line)表示按行讀取CSV檔案中的資料
    {
        istringstream sin(line); //將整行字串line讀入到字串流sin中

        getline(sin, field, ','); //將字串流sin中的字元讀入到field字串中，以逗號為分隔符 
	matrix_id.push_back(atoi(field.c_str()));


        getline(sin, field, ','); //將字串流sin中的字元讀入到field字串中，以逗號為分隔符 
	matrix_color.push_back(atof(field.c_str()));


        getline(sin, field, ','); //將字串流sin中的字元讀入到field字串中，以逗號為分隔符 
	matrix_algorithm.push_back(atof(field.c_str()));

        getline(sin, field); //將字串流sin中的字元讀入到field字串中，以逗號為分隔符 
	matrix_regression.push_back(atof(field.c_str()));
        i++;
    }
    inFile.close();
    cout << "共讀取了：" << i << "行" << endl;
    cout << "讀取資料完成" << endl;
}

int main()
{
  vector<int> matrix_id;
  vector<double> matrix_color;
  vector<double> matrix_algorithm;
  vector<double> matrix_regression;
  string filename = "../24color_data/5.csv";
  //readfile
  load_csv( matrix_id, matrix_color, matrix_algorithm, matrix_regression,  filename );
  double RMSE_vaule_t1 = 0;
  double RMSE_vaule_t2 = 0;
  int num = 0;
  if( matrix_color.size() == matrix_algorithm.size() ){ 
      num = matrix_color.size();
      RMSE_vaule_t1 = calcRMSE(matrix_color, matrix_algorithm, num);
  }
  else cout << " matrix_color size != matrix_algorithm size" << endl;

  if( matrix_color.size() == matrix_regression.size() ){ 
      num = matrix_color.size();
      RMSE_vaule_t2 = calcRMSE(matrix_color, matrix_regression, num);
  }
  else cout << " matrix_color size != matrix_regression size" << endl;
  cout << "color & algorithm RMSE : " << RMSE_vaule_t1 << endl;
  cout << "color & regression RMSE : " << RMSE_vaule_t2 << endl;

}



