#coding=utf-8
import multiprocessing as mp
import os, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import linear_model
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.datasets import make_regression
from sklearn.model_selection import train_test_split
import threading as td 
import os, time
from queue import Queue
import joblib



model_380 = joblib.load('380nm_degree4.pkl') #load模型
model_381 = joblib.load('381nm_degree4.pkl') #load模型
model_382 = joblib.load('382nm_degree4.pkl') #load模型

test_data = pd.read_csv("1x7.csv") #load圖片
r = test_data['R'].values
g = test_data['G'].values
b = test_data['B'].values
id = test_data['id'].values
series_dict = {'R':r,'G':g,'B':b}
df=pd.DataFrame(series_dict)
X=df[['R','G','B']]
X_test = PolynomialFeatures(4).fit_transform(X) #測試資料正規化





def job_380(x):
  print("pool", x)
  return model_380.predict(X_test)

def job_381(x):
  print("pool", x)
  return model_381.predict(X_test)

def job_382(x):
  print("pool", x)
  return model_382.predict(X_test)

def multicore():
    result=[[0]*7 for i in range(6)]
    cpus = os.cpu_count() # cpu核心數量
    pool = mp.Pool(processes=cpus) 
    
    result[0] = pool.map(job_380, (0,)) # 使用 pool 還可以接到 function 的回傳值 
    #print(result[0]) 

    result[1] = pool.map(job_381, (1,)) 
    #print(result[1]) 

    result[2] = pool.map(job_382, (2,)) 
    #print(result[2]) 

    result[3] = pool.map(job_380, (3,)) 
    #print(result[3]) 

    result[4] = pool.map(job_381, (4,)) 
    #print(result[4]) 

    result[5] = pool.map(job_382, (5,)) 
    #print(result[5])

    return  result




def main():
    #print(model_380.predict(X_test))
    print(multicore())

if __name__ == '__main__':
    main()

