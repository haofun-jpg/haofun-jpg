#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
using namespace std;

int main()
{
    ofstream outFile("Data.csv", ios::out);
    //ios::out：如果沒有檔案，那麼生成空檔案；如果有檔案，清空檔案
    if (!outFile)
    {
        cout << "開啟檔案失敗！" << endl;
        exit(1);
    }
    //寫入3行資料
    for (int i = 0; i < 3; i++)
    {
        outFile << 12 << ",";
        outFile << 13 << ",";
        outFile << 14 << ",";
        outFile << "NaN" << endl;
    }
    outFile.close();
    cout << "寫入資料完成" << endl;
}
