#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include<vector>
#include<cmath>
#include <string>
#include "math.h"
using namespace std;

void load_data_CPP( vector<float> &data, string filename ){
    float str;	
    fstream file;
    file.open(filename, ios::in);

    while( !file.eof( )){
      file >> str ;
      data.push_back( str );
    }
    file.close();
    data.pop_back() ; //最後一個
    cout << "vector size : " << data.size()  << endl;
}

void load_data_regression( vector<float> &data, string filename ){
    float str;	
    fstream file;
    file.open(filename, ios::in);

    while( !file.eof( )){
      file >> str ;
      data.push_back( str );
    }
    file.close();
    data.pop_back() ; //最後一個
    cout << "vector size : " << data.size()  << endl;
}

void calcMAE_div_y( vector<float> Data1, vector<float> Data2, int data_size, double &Max_MAE_not_zero_pct, double &MAE_vaule_not_zero_pct, 
			double &Max_MAE_zero, double &MAE_vaule_zero ){
	//初始化
	//data != 0
	Max_MAE_not_zero_pct = 0;
	MAE_vaule_not_zero_pct = 0;
	//data == 0
	Max_MAE_zero = 0;
	MAE_vaule_zero = 0;
	//data1[i] 1 != 0 數量
	double error_not_zero = 0;
	int num_not_zero = 0; 	
	//data1[i] 1 == 0 數量
	double error_not = 0;
	int num_zero = 0;
	

	int num_0_5 = 0;
	int num_6_10 = 0;
	int num_11_15 = 0;
	int num_16_20 = 0;
	int num_21_50 = 0;
	int num_50_N = 0;

	for( int i = 0; i < data_size; i++ ){
		if ( Data1[i] != 0 ){
			error_not_zero = fabs(Data1[i] - Data2[i]) /  Data1[i];  // 不能除以0會INF(無窮大)
			if( error_not_zero > Max_MAE_not_zero_pct ){

				cout << "------------Data1[i] != 0------------" << endl;
				cout << "Data1[i] : " << Data1[i] << endl;
				cout << "Data2[i] : " << Data2[i] << endl;
				cout << "-------------------------------------" << endl;

				Max_MAE_not_zero_pct = error_not_zero;

			}

			if( error_not_zero < 0.06 ) num_0_5++; // 0~5%
			else if ( error_not_zero >= 0.06 && error_not_zero < 0.11 ) num_6_10++; //6~10%
			else if ( error_not_zero >= 0.11 && error_not_zero < 0.16 ) num_11_15++; //11~15% 
			else if ( error_not_zero >= 0.16 && error_not_zero < 0.21 ) num_16_20++; //16~20%
			else if ( error_not_zero >= 0.21 && error_not_zero < 0.51 ) num_21_50++; //21~50%
			else if ( error_not_zero >= 0.51 ) num_50_N++; //51%~N%

			MAE_vaule_not_zero_pct += error_not_zero;
			num_not_zero++;
		}
		else{
			error_not = fabs(Data1[i] - Data2[i]);
			if( error_not > Max_MAE_zero ){

				cout << "------------Data1[i] == 0------------" << endl;
				cout << "Data1[i] : " << Data1[i] << endl;
				cout << "Data2[i] : " << Data2[i] << endl;
				cout << "-------------------------------------" << endl;

				Max_MAE_zero = error_not;
			} 
			MAE_vaule_zero += error_not;
			num_zero++;
		}
	}
	
	MAE_vaule_not_zero_pct /= num_not_zero;
	if( num_zero == 0 ){
		Max_MAE_zero = 0;
		MAE_vaule_zero = 0;
	}
	else {
		MAE_vaule_zero /= num_zero;
	}
	cout << "-----not error MAE PCT-----" << endl;
	cout << "0_5 : " << (float)num_0_5/num_not_zero << "%" << endl;
	cout << "6_10 : " << (float)num_6_10/num_not_zero << "%"  << endl;
	cout << "11_15 : " << (float)num_11_15/num_not_zero << "%"  << endl;
	cout << "16_20 : " << (float)num_16_20/num_not_zero << "%"  << endl;
	cout << "21_50 : " << (float)num_21_50/num_not_zero << "%"  << endl;
	cout << "51_N : " << (float)num_50_N/num_not_zero << "%"  << endl;
	cout << "not_zero: " << num_not_zero << endl;
	cout << "zero: " << num_zero << endl;

}

int main(){
	string file_cpp = "Output_data.txt";
	string file_regression = "regression_output.txt";
    	vector<float> data_CPP;
    	vector<float> data_regression;
    	load_data_CPP( data_CPP, file_cpp );
	load_data_regression( data_regression, file_regression );

	int data_size = 0;
	//data != 0
	double Max_MAE_not_zero_pct = 0;
	double MAE_vaule_not_zero_pct = 0;
	//data == 0
	double Max_MAE_zero = 0;
	double MAE_vaule_zero = 0;

	if( data_CPP.size() == data_regression.size() ){
		data_size = data_CPP.size();
		calcMAE_div_y( data_CPP, data_regression, data_size, Max_MAE_not_zero_pct, MAE_vaule_not_zero_pct, Max_MAE_zero, MAE_vaule_zero );
	}
	else cout << "not same size can't calcMAE !" << endl;

	data_CPP.clear();
	data_regression.clear();

	cout << "************************************" << endl;
	cout << "( Not zero ) percent of MAE is:" << MAE_vaule_not_zero_pct << endl;
	cout << "( Not zero ) percent of <MAX> is:" << Max_MAE_not_zero_pct << endl;

	cout << "( Is zero ) of MAE is:" << MAE_vaule_zero << endl;
	cout << "( Is zero ) of <MAX> is:" << Max_MAE_zero << endl;
	cout << "************************************" << endl;


}



