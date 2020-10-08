import os
os.chdir('D:/E/博士/R_程序/HGSOC/Data')
import numpy as np
from pandas import DataFrame


from sklearn.ensemble import AdaBoostRegressor as AdaBoost
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


# Adaboost 方法

# =============================================================================
 clf = AdaBoost()  # AdaBoost(n_estimators=100)
 rfe = RFE(clf, n_features_to_select=1) 
# =============================================================================


rfe.fit(esetarray, sampletype)

print("Features sorted by their rank:")
print(sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames)))
result = sorted(zip(map(lambda x: round(x, 4), rfe.ranking_), featurenames))
resultframe = DataFrame(result)
resultframe.to_csv("ranklist_ABrfe.txt", sep="\t")
