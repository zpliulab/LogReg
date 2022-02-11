#LogReg

**LogReg**: A method of regularized logistic regression for biomarker discovery from gene expression data. 

In this work, we also compared **LogReg** with the other five recursive feature elimination (**RFE**) feature selection methods, namely, AB-RFE, NN-RFE, RF-RFE, KNN-RFE and SVM-RFE. 

# Citation
Li, Lingyu, and Zhi-Ping Liu. "Biomarker discovery for predicting spontaneous preterm birth from gene expression data by regularized logistic regression." **Computational and Structural Biotechnology Journal** 18 (2020): 3434-3446.

@article{li2020biomarker,
  title={Biomarker discovery for predicting spontaneous preterm birth from gene expression data by regularized logistic regression},
  author={Li, Lingyu and Liu, Zhi-Ping},
  journal={Computational and Structural Biotechnology Journal},
  volume={18},
  pages={3434--3446},
  year={2020},
  publisher={Elsevier}
}


#Data

In the **Data** file, we only give some necessary files by each R program. Some of these input files only give the first few lines, but this does not affect the results of the work (**LogReg**).


#R code for LogReg

The **serial number (1) (2) ... (10)** represents the order in which the program runs in our work.

(1)GSE59491_expr.R ---- Processing original data of GSE59491.

(2)GSE73685_expr.R ---- Processing original data of GSE73685.

(3)Ttest.R ---- It contains some functions used to select differentially expressed genes (DEGs), and is included in “DEgene.R”. 

(4)DEgene.R ---- Identify the differentially expressed genes (DEGs) on a dataset and find candidates with adjusted P-value < 0.05.

(5)GSE59491_cel37_rep.R ---- Solve the regularized logistic regression with seven effective penalties, i.e., ridge, lasso, elastic net, L0, L1/2, SCAD and MCP.

(6)Feature_select.R ---- Extract data from independent dataset and origina discovery dataset based on identified biomarker.

(7)Class_ROC.R" ---- Verify the identified biomarkers on an independent dataset by the AUC value. (**Figure 5** a in our work)

(8)Box_PRS.R ---- Obtain the boxplot of of preterm risk score (PRS) on the independent dataset. (**Figure 5** b in our work)

(9)Hyper_test.R ---- Calculate the number of overlapping genes selected by hypergeometric test.  (**Table 4** in our work)

(10) stability.R ---- Calculate the stability of selected features/genes by seven different methods.



(8) "RFE_ROC_clear.R" ---- SVM classifier.

(9) "feature_overlap_clear.R" ---- SVM: Intersection of feature subsets obtained by 5 RFE methods to identify biomarkers.

(10) "feature_select_clear.R" ---- SVM: Extract features from the test set and verification set separately.

(11) "class_feature.R" ---- SVM: Perform independent data set verification and draw ROC curve.

(12) "cluster_clear.R" ---- SVM: Enrichment analysis of identified biomarkers.

(13) "ref_knn.R" ---- SVM: KNN-RFE.

(14) "ref_nnet.R.R" ---- SVM: NN-RFE.

(15) "rfsvm.py" ---- SVM: RF-RFE.

(16) "svmrfe.py" ---- SVM: SVM-RFE.

(17) "abrfe.py" ---- SVM: AB-RFE.
