from sklearn.preprocessing import PolynomialFeatures
from sklearn import linear_model
import numpy as np
import pandas as pd
data = pd.read_csv("smallnm.csv")
print(data)

r = data['R'].values
g = data['G'].values
b = data['B'].values
id = data['id'].values
spectrum = data['380nm'].values
series_dict = {'R':r,'G':g,'B':b,'380nm':spectrum}
df=pd.DataFrame(series_dict)
X=df[['R','G','B']]
y=df[['380nm']]

#train正規化
X_tt = PolynomialFeatures(degree = 2).fit_transform(X)
#模型訓練
model = linear_model.LinearRegression()
model.fit(X_tt, y)

print("model.intercept_: ")
print(model.intercept_)
print("model.coef_: ")
print(model.coef_)

test_data = pd.read_csv("1x7.csv") #load_test_圖片
r = test_data['R'].values
g = test_data['G'].values
b = test_data['B'].values
id = test_data['id'].values
series_dict = {'R':r,'G':g,'B':b}
df=pd.DataFrame(series_dict)
X=df[['R','G','B']]

X_test = PolynomialFeatures(2).fit_transform(X)
print("X_test正規化: ")
print(X_test)
print(model.predict(X_test))

'''
a=[[2,3,4]]
pf=PolynomialFeatures(degree=2)
print(pf.fit_transform(a))
print(pf.get_feature_names('R''G''B'))
pf=PolynomialFeatures(degree=2,include_bias=False)
print(pf.fit_transform(a))
pf=PolynomialFeatures(degree=2,interaction_only=True)
print(pf.fit_transform(a))
'''