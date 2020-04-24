# LogReg
LogReg proposed a comparative study of the regularized logistic regression with seven effective penalties, i.e., ridge, lasso, elastic net, L0, L1/2, SCAD and MCP, for the selection of strongly indicative genes from gene expression data. 

(1) "Box_PRS.R"  ----  Obtain the boxplot of of preterm risk score (PRS) on the independent dataset.

(2) "Class_ROC.R" ----   Verify the identified biomarkers on an independent dataset by the AUC value.
						 					 
(3) "DEgene.R" ----   Identify the differentially expression genes (DEGs) in a dataset and find candidates with adjusted P-value < 0.05.		
(4) "GSE59491_cel.R" ----  Solve the regularized logistic regression with seven effective penalties, i.e., ridge, lasso, elastic net, L0, L1/2, SCAD and MCP. 				  

(5) "GSE59491_expr.R" ----  Processing original data of GSE59491. 	

(6) "GSE73685_expr.R" ----  Processing original data of GSE73685. 	

(7) "Hyper_test.R" ----  Calculate the number of overlapping genes selected by hypergeometric test. 	
