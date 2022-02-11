import os
os.chdir('D:/E/博士/R_程序/HGSOC/Data')
import numpy as np
from pandas import DataFrame

from sklearn import svm

#from sklearn.neighbors import KNeighborsClassifier as KNN
#from sklearn.neighbors import KNeighborsRegressor as KNN

#from sklearn.linear_model import LogisticRegression as LR

#from sklearn.ensemble import RandomForestClassifier
#from sklearn.ensemble import RandomForestRegressor as RF

#from sklearn.ensemble import AdaBoostRegressor as AdaBoost

from sklearn.feature_selection import RFE

import pandas as pd
# ¶ÁÈë±í´ï¾ØÕó
dataframe = pd.read_table("matrix_DE.txt")
# »ñÈ¡Ñù±¾Ãû³Æ
samplenames = np.array(dataframe.columns)
# »ñÈ¡Ì½ÕëÃû³Æ
featurenames = np.array(dataframe.index)
# »ñÈ¡±í´ïÖµµÄ¾ØÕó
esetarray = np.array(dataframe)
# ±í´ïÖµµÄ¾ØÕó×ªÖÃ
esetarray = esetarray.transpose()

# ½«Ñù±¾·ÖÀà£¬
# ·Ö³ÉÈýÀà:PB BM LN
# sampletype = [1]*26+[2]*19+[3]*17
# ·Ö³ÉÁ½Àà:PB BMorLN
sampletype = [0]*10+[1]*10


# 支持向量机
# 支持向量机(support vector machine)是一种分类算法，但是也可以做回归，根据输入的数据不同可做不同的模型
#（若输入标签为连续值则做回归，若输入标签为分类值则用SVC()做分类）。
clf = svm.SVC(kernel='linear')  # clf.fit(X, y)
rfe = RFE(clf, n_features_to_select=1)



# 随机森林分类
# =============================================================================
# class RandomForestClassifierWithCoef(RandomForestClassifier):    
#     def fit(self, *args, **kwargs):        
#         super(RandomForestClassifierWithCoef, self).fit(*args, **kwargs)        
#         self.coef_= self.feature_importances_ 
# 
# ##clf = RandomForestClassifierWithCoef(n_estimators=500, min_samples_leaf=5, n_jobs=-1)
# ##rfe = RFECV(estimator=clf, step=1, scoring='accuracy',cv=2)
# #clf = RandomForestClassifierWithCoef()
# rfe = RFE(clf, n_features_to_select=1)
# =============================================================================

# 随机森林回归
# =============================================================================
# clf = RF()
# rfe = RFE(clf, n_features_to_select=1)
# =============================================================================


# K近邻
#clf = kNN(n_neighbors=5, algorithm='auto', weights='distance', n_jobs=1)
# =============================================================================
# clf = KNN()
# rfe = RFE(clf, n_features_to_select=1) 
# =============================================================================


#lf = LR()

# Adaboost 方法

# =============================================================================
# clf = AdaBoost()  # AdaBoost(n_estimators=100)
# rfe = RFE(clf, n_features_to_select=1) 
# =============================================================================


rfe.fit(esetarray, sampletype)

print("Features sorted by their rank:")
print(sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames)))
result = sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames))
resultframe = DataFrame(result)
resultframe.to_csv("ranklist_SVMrfe.txt", sep="\t")
