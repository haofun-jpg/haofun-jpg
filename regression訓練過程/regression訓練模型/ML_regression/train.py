import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import linear_model
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.datasets import make_regression
from sklearn.model_selection import train_test_split
import joblib



def main( index ):

	wavelength = str(index) + 'nm'
	print( "input : " + wavelength )
	path = "D:\RGBtable_any_nm\Project2" + '\\' + wavelength + ".csv"
	# '\\' = \
	data = pd.read_csv( path )
	
	r = data['R'].values
	g = data['G'].values
	b = data['B'].values
	id = data['id'].values
	spectrum = data[wavelength].values
	series_dict = {'R':r,'G':g,'B':b, wavelength:spectrum}
	df=pd.DataFrame(series_dict)
	X=df[['R','G','B']]
	y=df[[wavelength]]
	X_tt = PolynomialFeatures(4).fit_transform(X)
	model = linear_model.LinearRegression()
	model.fit(X_tt, y)
	joblib.dump(model, wavelength + '_degree4.pkl') #模型儲存
	print( "finish stored " + wavelength )
	

if __name__ == '__main__':
	for i in range(380,781):
		main( i )