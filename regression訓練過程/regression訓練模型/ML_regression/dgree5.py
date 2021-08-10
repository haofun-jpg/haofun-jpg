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
import multiprocessing as mp
from queue import Queue
import joblib





# 子Process的工作函數

def job_process(num):
  print("process", num)
  model.predict(X_test)

# 子Thread的工作函數

def job_Thread(q):
  print("Thread", q)
  model.predict(X_test)

def pool_job(q):
  print("pool", q)
  return model.predict(X_test)

def ttt():

    cpus = os.cpu_count() # cpu核心數量
    pool = mp.Pool(processes=cpus) 
    
    result = pool.map(pool_job, range(10)) # 使用 pool 還可以接到 function 的回傳值 
    #print(result) # [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]

def multicore():
    p1 = mp.Process(target=job_process, args=(1,)) # 特別注意 這邊的傳入參數只有一個的話，後面要有逗號
    p2 = mp.Process(target=job_process, args=(2,))
    p3 = mp.Process(target=job_process, args=(3,))
    p4 = mp.Process(target=job_process, args=(4,))
    p5 = mp.Process(target=job_process, args=(5,))
    p6 = mp.Process(target=job_process, args=(6,))
    p7 = mp.Process(target=job_process, args=(7,))
    p8 = mp.Process(target=job_process, args=(8,))
    p9 = mp.Process(target=job_process, args=(9,))
    p10 = mp.Process(target=job_process, args=(10,))
    


    p1.start()
    p2.start()
    p3.start()
    p4.start()
    p5.start()
    p6.start()
    p7.start()
    p8.start()
    p9.start()
    p10.start()
    

    p1.join()
    p2.join()
    p3.join()
    p4.join()
    p5.join()
    p6.join()
    p7.join()
    p8.join()
    p9.join()
    p10.join()

def multithread():
    t1 = td.Thread(target=job_Thread, args=(1,)) # 特別注意 這邊的傳入參數只有一個的話，後面要有逗號
    t2 = td.Thread(target=job_Thread, args=(2,))
    t3 = td.Thread(target=job_Thread, args=(3,))
    t4 = td.Thread(target=job_Thread, args=(4,))
    t5 = td.Thread(target=job_Thread, args=(5,))
    t6 = td.Thread(target=job_Thread, args=(6,)) 
    t7 = td.Thread(target=job_Thread, args=(7,))
    t8 = td.Thread(target=job_Thread, args=(8,))
    t9 = td.Thread(target=job_Thread, args=(9,))
    t10 = td.Thread(target=job_Thread, args=(10,))

    t1.start()
    t2.start()
    t3.start()
    t4.start()
    t5.start()
    t6.start()
    t7.start()
    t8.start()
    t9.start()
    t10.start()

    t1.join()
    t2.join()
    t3.join()
    t4.join()
    t5.join()
    t6.join()
    t7.join()
    t8.join()
    t9.join()
    t10.join()


if __name__ == "__main__":

  """# 讀取模型 - 380nm_degree4.pkl"""

  model = joblib.load('380nm_degree4.pkl') #load模型

  test_data = pd.read_csv("1x7.csv") #load圖片

  r = test_data['R'].values
  g = test_data['G'].values
  b = test_data['B'].values
  id = test_data['id'].values
  series_dict = {'R':r,'G':g,'B':b}
  df=pd.DataFrame(series_dict)
  X=df[['R','G','B']]
  X_test = PolynomialFeatures(4).fit_transform(X) #測試資料正規化


  print(model.predict(X_test))

  print(os.cpu_count()) 


  """# 測試時間 """
'''
  start = time.time()
  for i in range(10):
    model.predict(X_test)
  end = time.time()
  print (end-start)
  '''


  # 開始測量

  ttt()
  # 結束測量
 

'''
  start = time.time()
  multicore()
  end = time.time()
  print (end-start)

  start = time.time()
  multithread()
  end = time.clock()
  print (end-start)
  '''