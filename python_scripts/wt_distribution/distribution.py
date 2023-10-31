import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from fitter import Fitter, get_common_distributions, get_distributions

plt.rc('font', family='Arial')

with open("diff_wt.txt","r") as file:
    data=file.read()
#    print("The file contents are:")
#    print(data)
get_total = data.split("\n")
#print(get_total)
res = [eval(i) for i in get_total]
#print(res)
arr = np.array(res)
print(type(arr))
f = Fitter(arr, distributions= get_common_distributions())
f.fit()
f.summary()
plt.show()