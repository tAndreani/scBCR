# scBCR Benchmarking
Benchmarking computational methods for B cell receptor reconstruction from single cell RNA-seq data. This repository is an extension of [scBCR](https://gitlab.com/tAndreani/scBCR), a platform with all the methods currently available and ready to be installed and used to reconstruct paired haeavy and light chain in individual B cells.  

<p align="center">
<img src="https://user-images.githubusercontent.com/6462162/150326571-4ac5952c-b291-4a5d-9e9a-1a35110b51da.PNG" alt="Image" width="800" height="550" style="display: block; margin: 0 auto" />
  
**Fig 1 Benchmark framework**: **(A)** Available datasets with ground truth (Canzar et al. 2017, 
Upadhyay et al. 2018) consisting of plasmablasts with unknown antigen specificity were obtained from the corresponding publications. In addition, a scRNAseq dataset including tetanus toxoid-specific B-cells and B-cells with unknown antigen specificity was generated in this work (Leiden). All the datasets were downsampled to achieve different read coverage and length. **(B)** An additional dataset was simulated to investigate the effects of different levels of somatic hypermutations in the variable domains 
of BCRs on the performance of each method. Sensitivity and accuracy were used as metrics to evaluate each method depending on the type of used data. Antibodies were produced using a subset of clonotype-forming patient-specific BCRs and their specificity was experimentally validated. Finally, the execution time was investigated, and a final score was calculated to give a final recommendation on the method choice.  


[Fig 2 Performance methods on experimental data: Leiden, Canzar, Updahyay](https://github.com/tAndreani/scBCR/blob/main/Scripts/Plot_Sensitivity.r)  

[Fig 3 Pairwise overlap of V-(D)-J genes](https://github.com/tAndreani/scBCR/blob/main/Scripts/Overlap_Genes/Figure_3_Lower_Heatmap.r)  
  
[Fig 4 Performance methods on simulated data: 4 datasets with different amount of somatic hypermutations (from 15 to 60)](https://github.com/tAndreani/scBCR/blob/main/Scripts/PLOT_accuracy_SHMs.r)

[Fig 6 Time of execution for each method in the Leiden dataset](https://github.com/tAndreani/scBCR/blob/main/Scripts/Plot_Time.r)

<p align="center"><img src="Images/Fig.6-A.PNG" alt="Image" width="600" height="450">
<p align="center"><img src="Images/Fig.6-B.PNG" alt="Image" width="700" height="450">

[Fig 7 Heatmap and decision tree to evaluate the performance of each method and recommendation on different datasets](https://github.com/tAndreani/scBCR/blob/main/Scripts/HeatMap_All_Tools_Evaluation.r)

#### Credits:  
Tommaso Andreani, Linda Slot, Samuel Gabillard, Carsten Struebing, Claus Reimertz, Yaligara Veeranagouda, Aleida Bakker, Reza Olfati-Saber, Rene E. M. Toes, Franck Auge´, Deimante Simaite.  [Benchmarking computational methods for B-cell receptor reconstruction from single-cell RNA-seq data](https://www.biorxiv.org/content/10.1101/2022.03.24.485600v2) in _bioRxiv_.
