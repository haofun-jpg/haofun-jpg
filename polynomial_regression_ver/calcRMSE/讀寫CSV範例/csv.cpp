#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
using namespace std;

int main()
{
    //ofstream outFile("Data.csv", ios::out);
    ////ios::out：如果沒有檔案，那麼生成空檔案；如果有檔案，清空檔案
    //if (!outFile)
    //{
    //  cout << "開啟檔案失敗！" << endl;
    //  exit(1);
    //}
    ////寫入3行資料
    //for (int i = 0; i < 3; i++)
    //{
    //  outFile << 12 << ",";
    //  outFile << 13 << ",";
    //  outFile << 14 << ",";
    //  outFile << "NaN" << endl;
    //}
    //outFile.close();
    //cout << "寫入資料完成" << endl;

    ifstream inFile("Data.csv", ios::in);
    if (!inFile)
    {
        cout << "開啟檔案失敗！" << endl;
        exit(1);
    }
    int i = 0;
    string line;
    string field;

    while (getline(inFile, line))//getline(inFile, line)表示按行讀取CSV檔案中的資料
    {
        string field;
        istringstream sin(line); //將整行字串line讀入到字串流sin中

        getline(sin, field, ','); //將字串流sin中的字元讀入到field字串中，以逗號為分隔符 
        cout<<atoi(field.c_str())<<" ";//將剛剛讀取的字串轉換成int

        getline(sin, field, ','); //將字串流sin中的字元讀入到field字串中，以逗號為分隔符 
        cout << atoi(field.c_str()) << " ";//將剛剛讀取的字串轉換成int

        getline(sin, field, ','); //將字串流sin中的字元讀入到field字串中，以逗號為分隔符 
        cout << atoi(field.c_str()) << " ";//將剛剛讀取的字串轉換成int

        getline(sin, field); //將字串流sin中的字元讀入到field字串中，以逗號為分隔符 
        cout << field.c_str() << endl;//將剛剛讀取的字串轉換成int
        i++;
    }
    inFile.close();
    cout << "共讀取了：" << i << "行" << endl;
    cout << "讀取資料完成" << endl;
}
