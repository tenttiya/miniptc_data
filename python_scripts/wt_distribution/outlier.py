import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from fitter import Fitter, get_common_distributions, get_distributions

with open("diff_wt.txt","r") as file:
    data=file.read()
#    print("The file contents are:")
#    print(data)
get_total = data.split("\n")
#print(get_total)
res = [eval(i) for i in get_total]
#print(res)
arr = np.array(res)
#print(type(arr))
#f = Fitter(arr, distributions= ['cauchy'])
#f.fit()
#f.summary()
#plt.show()
#print(f.fitted_param['cauchy'])

def get_outlier_bounds(l, tol=1.5):
    q1 = np.quantile(l, 0.25)
    q3 = np.quantile(l, 0.75)
   # print(q1, q3)
    iqr = q3 - q1
    high = q3 + tol*iqr
    low = q1 - tol*iqr
    return low, high

def get_outlier_list(l, tol=1.5):
    q1 = np.quantile(l, 0.25)
    q3 = np.quantile(l, 0.75)
   # print(q1, q3)
    iqr = q3 - q1
    high = q3 + tol*iqr
    low = q1 - tol*iqr
    upper = []
    lower = []
    for i in l:
        if i > high: 
            upper.append(i)
        if i < low:
            lower.append(i)
    return upper, lower


print(get_outlier_bounds(arr, tol=3))

print(get_outlier_list(arr, tol=3))