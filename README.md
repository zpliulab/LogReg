# [LogReg (regularized logistic regression)](https://github.com/zpliulab/LogReg)

[![Screenshot](https://ars.els-cdn.com/content/image/1-s2.0-S2001037020304505-gr1.jpg)](https://doi.org/10.1016/j.csbj.2020.10.028)

In this work, we provide **a computational method of regularized logistic regression** for discovering biomarkers of **spontaneous preterm birth (SPTB)** from gene expression data. The successful identification of SPTB biomarkers will greatly benefit the interference of infant gestational age for reducing the risks of pregnant women and preemies. **Obviously, the proposed method of discovering biomarkers for SPTB can be easily extended for other complex diseases**.


## LogReg
<!--START_SECTION:news-->
* **LogReg**: A method of regularized logistic regression for biomarker discovery from gene expression data. 
* In this work, we also compared **LogReg** with the other five recursive feature elimination (**RFE**) feature selection methods, namely, AB-RFE, NN-RFE, RF-RFE, KNN-RFE and SVM-RFE. 
* If you have any questions about **LogReg**, please directly contact the corresponding author [Prof. Zhi-Ping Liu](https://scholar.google.com/citations?user=zkBXb_kAAAAJ&hl=zh-CN&oi=ao) with the E-mail: zpliu@sdu.edu.cn

<!--END_SECTION:news-->


## Citation
Li, Lingyu, and Zhi-Ping Liu. "**Biomarker discovery for predicting spontaneous preterm birth from gene expression data by regularized logistic regression**." Computational and Structural Biotechnology Journal 18 (2020): 3434-3446. [LogReg paper website](https://doi.org/10.1016/j.csbj.2020.10.028)


## Data
<!--START_SECTION:news-->
* In the **Data** and **Data37** files, we only give some necessary files by each R program. 
* Some of these input files only give the first few lines, but this does not affect the results of the work (**LogReg**).
<!--END_SECTION:news-->


## R code for LogReg
The **serial number (1) (2) ... (10)** represents the order in which the program runs in our work.
<!--START_SECTION:news-->
* (1) GSE59491_expr.R ---- Processing original data of GSE59491.
* (2) GSE73685_expr.R ---- Processing original data of GSE73685.
* (3) Ttest.R ---- It contains some functions used to select differentially expressed genes (DEGs), and is included in “DEgene.R”. 
* (4) DEgene.R ---- Identify the differentially expressed genes (DEGs) on a dataset and find candidates with adjusted P-value < 0.05.
* (5) GSE59491_cel37_rep.R ---- Solve the regularized logistic regression with seven effective penalties, i.e., ridge, lasso, elastic net, L0, L1/2, SCAD and MCP.
* (6) Feature_select.R ---- Extract data from independent dataset and origina discovery dataset based on identified biomarker.
* (7) Class_ROC.R" ---- Verify the identified biomarkers on an independent dataset by the AUC value. (**Figure 5** a in our work)
* (8) Box_PRS.R ---- Obtain the boxplot of of preterm risk score (PRS) on the independent dataset. (**Figure 5** b in our work)
* (9) Hyper_test.R ---- Calculate the number of overlapping genes selected by hypergeometric test.  (**Table 4** in our work)
* (10)  Stability.R ---- Calculate the stability of selected features/genes by seven different methods.
<!--END_SECTION:news-->


## Stability of Feature Selection Techniques for Bioinformatics Dat
<!--START_SECTION:news-->
* [Feature selection is one of the most fundamental problems in data analysis, machine learning, and data mining](https://doi.org/10.1007/978-3-030-64583-0_19). Especially in domains where the chosen features are subject to further experimental research, the stability of the feature selection is very important. Stable feature selection means that the set of selected features is robust with respect to different data sets from the same data generating distribution.
* For data sets with similar features, the evaluation of feature selection stability is more difficult. An example of such data sets is gene expression data sets, where genes of the same biological processes are often highly positively correlated.  Here, stability measures that take into account the similarities between feature subsets are defined as **Hamming stability**.
* We identified **all combinations** containing 4 sets from the 7 candidate sets, and calculated the **stability index (stabilityHamming)** of 35 combinations, the results are shown in the barplot obtained by R code **StabilityLogReg.R**.
* The **red bar** is the stability index of our chosen subset of features (2 3 5 6, i.e. **Lasso, Elastic net, L1/2, SCAD**), which ranks **fifth out of all 35 combinations**.
* This shows that the subset of features we choose not only has **high accuracy/AUC**, but also has **high stability**.
<!--END_SECTION:news-->


## RcodeForRFE
R code for five **RFE** feature selection methods.
<!--START_SECTION:news-->
* (1)"svmrfe.py" ---- SVM: SVM-RFE.
* (2) "abrfe.py" ---- SVM: AB-RFE.
* (3)"rfe_nnet.R" ---- SVM: NN-RFE.
* (4)"rfsvm.py" ---- SVM: RF-RFE.
* (5)"ref_knn.R" ---- SVM: KNN-RFE.
* (6)"feature_overlap_clear.R" ---- SVM: Intersection of feature subsets obtained by 5 RFE methods to identify biomarkers.
* (7)"feature_select_clear.R" ---- SVM: Extract features from the test set and verification set separately.
* (8)"RFE_ROC_clear.R" ---- SVM classifier.
* (9)"class_feature_clear.R" ---- SVM: Perform independent data set verification and draw ROC curve.
* (10) "cluster_clear.R" ---- SVM: Enrichment analysis of identified biomarkers.
<!--END_SECTION:news-->


## LogReg (2020), Zhi-Ping Liu all rights reserved
This program package is supported by the copyright owners and coders "as is" and without warranty of any kind, express or implied, including, but not limited to, the implied warranties of merchantability and fitness for a particular purpose. In no event shall the copyright owner or contributor be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, without limitation, procurement of substitute goods or services; loss of use, data, or profits; or business interruption), regardless of the theory of liability, whether in contract, strict liability or tort (including negligence or otherwise) for any use of the software, even if advised of the possibility of such damages.
