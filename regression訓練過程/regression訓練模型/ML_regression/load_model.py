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

def main( index ):
  file_name = str(index) + 'nm_degree4.pkl'
  model = joblib.load(file_name) #load模型

  path = "C_degree4" + '\\'
  store_name = path + str(index) + 'nm_degree4.txt'
  f = open(store_name, "w")
  print(model.intercept_[0],file = f)
  print(file = f) #換行
  for i in range(len(model.coef_[0])) :
    print(model.coef_[0][i], file = f)
  f.close()
  
if __name__ == "__main__":
	for i in range(380,781):
		main( i )
