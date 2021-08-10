#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include<vector>
#include "math.h"
using namespace std;



void load_data_CPP( vector<float> &data ){
    float str;	
    fstream file;
    file.open("Output_data.txt", ios::in);

    while( !file.eof( )){
      file >> str ;
      data.push_back( str );
    }
    file.close();
    data.pop_back() ; //最後一個
    cout << "vector size : " << data.size()  << endl;
}

void load_data_regression( vector<float> &data ){
    float str;	
    fstream file;
    file.open("regression_output.txt", ios::in);

    while( !file.eof( )){
      file >> str ;
      data.push_back( str );
    }
    file.close();
    data.pop_back() ; //最後一個
    cout << "vector size : " << data.size()  << endl;
}

double calcRMSE(vector<float> Data,vector<float> Data2, int Num)
{

	double fSum = 0;
	//#pragma omp parallel for
	for (int i = 0; i < Num; ++i)
	{
		fSum += (Data[i] - Data2[i]) *(Data[i] - Data2[i]);
	}
	return sqrt(fSum / Num);

}

int main(){
    vector<float> data_CPP;
    vector<float> data_regression;
    load_data_CPP( data_CPP );
    load_data_regression( data_regression );
    int num = 0 ;
    double RMSE_vaule = 0;
    if( data_CPP.size() == data_regression.size() ){ 
      num = data_CPP.size();
      RMSE_vaule = calcRMSE(data_CPP, data_regression, num);
    }
    else cout << "not same size can't calcRMSE !" << endl;

    cout << "the rmse of dataREAL and check is:" << RMSE_vaule << endl;
}



