
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV 

# 省略了训练集，测试集的划分 x_train,y_train,x_test,y_test  'RandomForestClassifier 没有属性 coef，将RandomForestClassifier改进以适用于RFECV所做的工作'

class RandomForestClassifierWithCoef(RandomForestClassifier):    
    def fit(self, *args, **kwargs):        
        super(RandomForestClassifierWithCoef, self).fit(*args, **kwargs)        
        self.coef_= self.feature_importances_ 
        
rf = RandomForestClassifierWithCoef(n_estimators=500, min_samples_leaf=5, n_jobs=-1)
rfecv = RFECV(estimator=rf, step=1, scoring='accuracy',cv=2)
selector = rfecv.fit(x_train, y_train)
print('RFECV 选择出的特征个数 ：' , rfecv.n_features_)  # RFECV选择的特征个数
print('特征优先级 ： ', rfecv.ranking_)        # 1代表选择的特征
x_train_rfecv = rfecv.transform(x_train)
x_test_rfecv = rfecv.transform(x_test)

# 随机森林
rf_clf = RandomForestClassifier(n_jobs=-1, max_depth=100, n_estimators=800)
rf_clf.fit(x_train_rfecv, y_train)
y_pred = rf_clf.predict(x_test_rfecv) 

# 模型评价
p = precision_score(y_test, y_pred, average='macro')
r = recall_score(y_test, y_pred, average='macro')
f1 = f1_score(y_test, y_pred, average='macro')
print('precision ： ' + str(p))
print('recall ： ' + str(r))
print('f1 : ' + str(f1)) 

# 可视化
# 不同特征数量下的评分
plt.figure()
plt.xlabel("Number of features selected")
plt.ylabel("Cross validation score ")
plt.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_)
plt.show()
