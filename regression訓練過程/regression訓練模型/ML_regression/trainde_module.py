import pandas as pd
import joblib



def load():
	
	model_380 = joblib.load('380nm_degree4.pkl') #load模型
	model_381 = joblib.load('381nm_degree4.pkl') #load模型
	model_382 = joblib.load('382nm_degree4.pkl') #load模型

