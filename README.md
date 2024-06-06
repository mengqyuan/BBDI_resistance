# Dissection of BET Bromodomain Inhibitor Resistance in Triple-Negative Breast Cancer

BET bromodomain inhibitor (BBDI) holds promise for cancer therapy, yet drug resistance poses a significant hurdle. We conducted an extensive investigation into BBDI resistance in Triple-Negative Breast Cancer (TNBC) by integrating single-cell RNA sequencing (scRNA-seq) and multi-omics analysis. Our results revealed cellular dynamic changes in the evolution of BBDI resistance. scRNA-seq analysis indicated that ferroptosis inhibition contributed to BBDI resistance. In vitro experiments demonstrated that GPX4 overexpression led to decreased inhibition rates and elevated IC50 in TNBC cells. We also proposed a FR20 algorithm to quantify BBDI resistance by integrating 20 ferroptosis regulators’ gene expression. Further validation by nine independent datasets confirmed the accuracy and scalability of FR20. Re-clustering analysis of scRNA-seq data revealed five resistant sub-clusters, illustrating heterogeneity with a temporal evolution and diverse drug responses to BBDI. Utilizing computational screening (D-FR20 method) and experimental validation, the small molecule filgotinib showed the ability to re-sensitize BBDI and effectively eliminate resistant TNBC cells. Collectively, this study provided novel insights into BBDI resistance, offering a promising avenue for TNBC treatment.

### ** Data availability **

The data sets analyzed in this study are available from the Gene Expression Omnibus (GEO) repository, Cancer Therapeutics Response Portal (CTRP), Cancer Cell Line Encyclopedia (CCLE), Genomics of Drug Sensitivity in Cancer (GDSC) and LINCS L1000 . GEO data sets are under the following accession numbers: GSE131135, GSE164813, GSE103082 available from GEO repository (https://www.ncbi.nlm.nih.gov/geo/). CTRP data sets include drug responses of JQ1 for cell lines available from CTRP (https://portals.broadinstitute.org/ctrp/). CCLE and GDSC data sets include gene expression profiles of cell lines available from depmap portal (https://depmap.org/portal/data_page/?tab=allData) and GDSC (https://www.cancerrxgene.org/). LINCS L1000 includes drug-perturbed gene expression profiles of 720216 small molecule instances available from NIH LINCS (https://lincsproject.org/LINCS/tools/workflows/find-the-best-place-to-obtain-the-lincs-l1000-data).

### ** Function of codes **

Fig2.r - Code to replicate the plots in Figure 2a-g.  
Fig3.r - Code to replicate the plots in Figure 3a and 3b.  
Fig4.r - Code to replicate the plots in Figure 4a-i.  
Fig5.r - Code to replicate the plots in Figure 5a, b, d-h.  
CaDRReS.py - Code to replicate the plots in Figure 5c.  
D-FR20.py - Code of D-FR20 method to identify small molecules to reverse BBDI resistance, as shown in Figure 6a.  
Fig7a.py - Code to replicate the plots in Figure 7a.  
Fig7b.r - Code to replicate the plots in Figure 7b.  
FigS1.r - Code to replicate the plots in Figure S1.  
FigS2.r - Code to replicate the plots in Figure S2.  
FigS4.r - Code to replicate the plots in Figure S4c-e.  
FigS5.r - Code to replicate the plots in Figure S5.  
FigS6.r - Code to replicate the plots in Figure S6.  

### Contact	

All comments, questions and suggestions: weijiang@nuaa.edu.cn
