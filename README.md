# LogReg
LogReg: A method of regularized logistic regression for biomarker discovery from gene expression data and compare with other 5 RFE methods. In the Data_SPTB folder, we only give examples of the file format input by each R program. The input file only gives the first few lines, but this does not affect the results of the article.

(1) "Box_PRS.R"  ----  Obtain the boxplot of of preterm risk score (PRS) on the independent dataset.

(2) "Class_ROC.R" ----   Verify the identified biomarkers on an independent dataset by the AUC value.
						 					 
(3) "DEgene.R" ----   Identify the differentially expression genes (DEGs) on a dataset 
						and find candidates with adjusted P-value < 0.05.		

(4) "GSE59491_cel.R" ----  Solve the regularized logistic regression with seven effective penalties,
							i.e., ridge, lasso, elastic net, L0, L1/2, SCAD and MCP. 				  

(5) "GSE59491_expr.R" ----  Processing original data of GSE59491. 	

(6) "GSE73685_expr.R" ----  Processing original data of GSE73685. 	

(7) "Hyper_test.R" ----  Calculate the number of overlapping genes selected by hypergeometric test. 

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
