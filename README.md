# Dominguez-Cell-2023-DataSources
Data and scripts for analyzing early cardiac development, to be published as "Graded mesoderm assembly governs cell fate and morphogenesis of the early mammalian heart".
<br><br>

## References to outside datasets analyzed herein
de Soysa, T.Y., Ranade, S.S., Okawa, S., Ravichandran, S., Huang, Y., Salunga, H.T., Schricker, A., Del Sol, A., Gifford, C.A., and Srivastava, D. (2019). Single-cell analysis of cardiogenesis reveals basis for organ-level developmental defects. Nature 572, 120–124.
<br><br>
Tyser, R.C.V., Ibarra-Soria, X., McDole, K., Arcot Jayaram, S., Godwin, J., van den Brand, T.A.H., Miranda, A.M.A., Scialdone, A., Keller, P.J., Marioni, J.C., et al. (2021). Characterization of a common progenitor pool of the epicardium and myocardium.
<br><br>
Gao, R., Liang, X., Cheedipudi, S., Cordero, J., Jiang, X., Zhang, Q., Caputo, L., Günther, S., Kuenne, C., Ren, Y., et al. (2019). Pioneering function of Isl1 in the epigenetic control of cardiomyocyte cell fate. Cell Res. 29, 486–501.
<br><br>
Quaranta, R., Fell, J., Rühle, F., Rao, J., Piccini, I., Araúzo-Bravo, M.J., Verkerk, A.O., Stoll, M., and Greber, B. (2018). Revised roles of ISL1 in a hES cell-based model of human heart chamber specification. eLife 7.
<br><br>

## R sessionInfo()
`R version 4.0.2 (2020-06-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scater_1.22.0               scuttle_1.4.0               SingleCellExperiment_1.16.0 SummarizedExperiment_1.24.0 GenomicRanges_1.46.1       
 [6] GenomeInfoDb_1.30.1         MatrixGenerics_1.6.0        matrixStats_0.62.0          biomaRt_2.50.3              cowplot_1.1.1              
[11] igraph_1.3.4                future_1.27.0               readxl_1.3.1                data.table_1.14.2           scales_1.2.1               
[16] GOplot_1.0.2                gridExtra_2.3               ggdendro_0.1.23             org.Mm.eg.db_3.14.0         topGO_2.46.0               
[21] SparseM_1.81                GO.db_3.14.0                AnnotationDbi_1.56.2        IRanges_2.28.0              S4Vectors_0.32.4           
[26] Biobase_2.54.0              graph_1.72.0                BiocGenerics_0.40.0         stringi_1.7.8               sp_1.5-0                   
[31] SeuratObject_4.1.0          Seurat_4.1.1                ggpubr_0.4.0                circular_0.4-95             ggridges_0.5.3             
[36] extrafont_0.17              RColorBrewer_1.1-3          reshape2_1.4.4              patchwork_1.1.2             viridis_0.6.2              
[41] viridisLite_0.4.0           forcats_0.5.1               stringr_1.4.0               dplyr_1.0.9                 purrr_0.3.4                
[46] readr_2.1.2                 tidyr_1.2.0                 tibble_3.1.8                tidyverse_1.3.1             readODS_1.7.0              
[51] ggplot2_3.3.6              

loaded via a namespace (and not attached):
  [1] utf8_1.2.2                reticulate_1.25           tidyselect_1.1.2          RSQLite_2.2.16            htmlwidgets_1.5.4         BiocParallel_1.28.3      
  [7] grid_4.2.2                Rtsne_0.16                ScaledMatrix_1.2.0        munsell_0.5.0             codetools_0.2-18          ica_1.0-3                
 [13] miniUI_0.1.1.1            withr_2.5.0               spatstat.random_2.2-0     colorspace_2.0-3          progressr_0.10.1          filelock_1.0.2           
 [19] rstudioapi_0.13           ROCR_1.0-11               ggsignif_0.6.3            tensor_1.5                Rttf2pt1_1.3.10           listenv_0.8.0            
 [25] GenomeInfoDbData_1.2.7    polyclip_1.10-0           bit64_4.0.5               parallelly_1.32.1         vctrs_0.4.1               generics_0.1.3           
 [31] BiocFileCache_2.2.1       R6_2.5.1                  ggbeeswarm_0.6.0          rsvd_1.0.5                DelayedArray_0.20.0       bitops_1.0-7             
 [37] spatstat.utils_2.3-1      cachem_1.0.6              assertthat_0.2.1          promises_1.2.0.1          beeswarm_0.4.0            rgeos_0.5-9              
 [43] gtable_0.3.0              beachmat_2.10.0           globals_0.16.0            goftest_1.2-3             rlang_1.0.4               splines_4.2.2            
 [49] rstatix_0.7.0             extrafontdb_1.0           lazyeval_0.2.2            spatstat.geom_2.4-0       broom_0.7.12              abind_1.4-5              
 [55] modelr_0.1.8              backports_1.4.1           httpuv_1.6.5              tools_4.2.2               ellipsis_0.3.2            spatstat.core_2.4-4      
 [61] Rcpp_1.0.9                plyr_1.8.7                sparseMatrixStats_1.6.0   progress_1.2.2            zlibbioc_1.40.0           RCurl_1.98-1.8           
 [67] prettyunits_1.1.1         rpart_4.1.19              deldir_1.0-6              pbapply_1.5-0             zoo_1.8-10                haven_2.4.3              
 [73] ggrepel_0.9.1             cluster_2.1.4             fs_1.5.2                  magrittr_2.0.3            scattermore_0.8           lmtest_0.9-40            
 [79] reprex_2.0.1              RANN_2.6.1                mvtnorm_1.1-3             fitdistrplus_1.1-8        hms_1.1.2                 mime_0.12                
 [85] xtable_1.8-4              XML_3.99-0.10             compiler_4.2.2            KernSmooth_2.23-20        crayon_1.5.1              htmltools_0.5.3          
 [91] mgcv_1.8-41               later_1.3.0               tzdb_0.3.0                lubridate_1.8.0           DBI_1.1.3                 dbplyr_2.1.1             
 [97] rappdirs_0.3.3            MASS_7.3-58               boot_1.3-28               Matrix_1.5-1              car_3.1-0                 cli_3.3.0                
[103] parallel_4.2.2            pkgconfig_2.0.3           plotly_4.10.0             spatstat.sparse_2.1-1     xml2_1.3.3                vipor_0.4.5              
[109] XVector_0.34.0            rvest_1.0.2               digest_0.6.29             sctransform_0.3.3         RcppAnnoy_0.0.19          spatstat.data_2.2-0      
[115] Biostrings_2.62.0         cellranger_1.1.0          leiden_0.4.2              uwot_0.1.13               DelayedMatrixStats_1.16.0 curl_4.3.2               
[121] shiny_1.7.2               lifecycle_1.0.1           nlme_3.1-160              jsonlite_1.8.0            BiocNeighbors_1.12.0      carData_3.0-5            
[127] fansi_1.0.3               pillar_1.8.1              lattice_0.20-45           KEGGREST_1.34.0           fastmap_1.1.0             httr_1.4.4               
[133] survival_3.4-0            glue_1.6.2                png_0.1-7                 bit_4.0.4                 blob_1.2.3                BiocSingular_1.10.0      
[139] memoise_2.0.1             irlba_2.3.5               future.apply_1.9.0
`

<sub>R version 4.0.2 (2020-06-22)</sub><br>
<sub>Platform: x86_64-pc-linux-gnu (64-bit)</sub><br>
<sub>Running under: Ubuntu 20.04.5 LTS</sub><br>
<sub></sub><br>
<sub>Matrix products: default</sub><br>
<sub>BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0</sub><br>
<sub>LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0</sub><br>
<sub></sub><br>
<sub>locale:</sub><br>
<sub> [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   </sub><br>
<sub> [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            </sub><br>
<sub>[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       </sub><br>
<sub></sub><br>
<sub>attached base packages:</sub><br>
<sub>[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     </sub><br>
<sub></sub><br>
<sub>other attached packages:</sub><br>
<sub> [1] scater_1.22.0               scuttle_1.4.0               SingleCellExperiment_1.16.0 SummarizedExperiment_1.24.0 GenomicRanges_1.46.1       </sub><br>
<sub> [6] GenomeInfoDb_1.30.1         MatrixGenerics_1.6.0        matrixStats_0.62.0          biomaRt_2.50.3              cowplot_1.1.1              </sub><br>
<sub>[11] igraph_1.3.4                future_1.27.0               readxl_1.3.1                data.table_1.14.2           scales_1.2.1               </sub><br>
<sub>[16] GOplot_1.0.2                gridExtra_2.3               ggdendro_0.1.23             org.Mm.eg.db_3.14.0         topGO_2.46.0               </sub><br>
<sub>[21] SparseM_1.81                GO.db_3.14.0                AnnotationDbi_1.56.2        IRanges_2.28.0              S4Vectors_0.32.4           </sub><br>
<sub>[26] Biobase_2.54.0              graph_1.72.0                BiocGenerics_0.40.0         stringi_1.7.8               sp_1.5-0                   </sub><br>
<sub>[31] SeuratObject_4.1.0          Seurat_4.1.1                ggpubr_0.4.0                circular_0.4-95             ggridges_0.5.3             </sub><br>
<sub>[36] extrafont_0.17              RColorBrewer_1.1-3          reshape2_1.4.4              patchwork_1.1.2             viridis_0.6.2              </sub><br>
<sub>[41] viridisLite_0.4.0           forcats_0.5.1               stringr_1.4.0               dplyr_1.0.9                 purrr_0.3.4                </sub><br>
<sub>[46] readr_2.1.2                 tidyr_1.2.0                 tibble_3.1.8                tidyverse_1.3.1             readODS_1.7.0              </sub><br>
<sub>[51] ggplot2_3.3.6              </sub><br>
<sub></sub><br>
<sub>loaded via a namespace (and not attached):</sub><br>
<sub>  [1] utf8_1.2.2                reticulate_1.25           tidyselect_1.1.2          RSQLite_2.2.16            htmlwidgets_1.5.4         BiocParallel_1.28.3      </sub><br>
<sub>  [7] grid_4.2.2                Rtsne_0.16                ScaledMatrix_1.2.0        munsell_0.5.0             codetools_0.2-18          ica_1.0-3                </sub><br>
<sub> [13] miniUI_0.1.1.1            withr_2.5.0               spatstat.random_2.2-0     colorspace_2.0-3          progressr_0.10.1          filelock_1.0.2           </sub><br>
<sub> [19] rstudioapi_0.13           ROCR_1.0-11               ggsignif_0.6.3            tensor_1.5                Rttf2pt1_1.3.10           listenv_0.8.0            </sub><br>
<sub> [25] GenomeInfoDbData_1.2.7    polyclip_1.10-0           bit64_4.0.5               parallelly_1.32.1         vctrs_0.4.1               generics_0.1.3           </sub><br>
<sub> [31] BiocFileCache_2.2.1       R6_2.5.1                  ggbeeswarm_0.6.0          rsvd_1.0.5                DelayedArray_0.20.0       bitops_1.0-7             </sub><br>
<sub> [37] spatstat.utils_2.3-1      cachem_1.0.6              assertthat_0.2.1          promises_1.2.0.1          beeswarm_0.4.0            rgeos_0.5-9              </sub><br>
<sub> [43] gtable_0.3.0              beachmat_2.10.0           globals_0.16.0            goftest_1.2-3             rlang_1.0.4               splines_4.2.2            </sub><br>
<sub> [49] rstatix_0.7.0             extrafontdb_1.0           lazyeval_0.2.2            spatstat.geom_2.4-0       broom_0.7.12              abind_1.4-5              </sub><br>
<sub> [55] modelr_0.1.8              backports_1.4.1           httpuv_1.6.5              tools_4.2.2               ellipsis_0.3.2            spatstat.core_2.4-4      </sub><br>
<sub> [61] Rcpp_1.0.9                plyr_1.8.7                sparseMatrixStats_1.6.0   progress_1.2.2            zlibbioc_1.40.0           RCurl_1.98-1.8           </sub><br>
<sub> [67] prettyunits_1.1.1         rpart_4.1.19              deldir_1.0-6              pbapply_1.5-0             zoo_1.8-10                haven_2.4.3              </sub><br>
<sub> [73] ggrepel_0.9.1             cluster_2.1.4             fs_1.5.2                  magrittr_2.0.3            scattermore_0.8           lmtest_0.9-40            </sub><br>
<sub> [79] reprex_2.0.1              RANN_2.6.1                mvtnorm_1.1-3             fitdistrplus_1.1-8        hms_1.1.2                 mime_0.12                </sub><br>
<sub> [85] xtable_1.8-4              XML_3.99-0.10             compiler_4.2.2            KernSmooth_2.23-20        crayon_1.5.1              htmltools_0.5.3          </sub><br>
<sub> [91] mgcv_1.8-41               later_1.3.0               tzdb_0.3.0                lubridate_1.8.0           DBI_1.1.3                 dbplyr_2.1.1             </sub><br>
<sub> [97] rappdirs_0.3.3            MASS_7.3-58               boot_1.3-28               Matrix_1.5-1              car_3.1-0                 cli_3.3.0                </sub><br>
<sub>[103] parallel_4.2.2            pkgconfig_2.0.3           plotly_4.10.0             spatstat.sparse_2.1-1     xml2_1.3.3                vipor_0.4.5              </sub><br>
<sub>[109] XVector_0.34.0            rvest_1.0.2               digest_0.6.29             sctransform_0.3.3         RcppAnnoy_0.0.19          spatstat.data_2.2-0      </sub><br>
<sub>[115] Biostrings_2.62.0         cellranger_1.1.0          leiden_0.4.2              uwot_0.1.13               DelayedMatrixStats_1.16.0 curl_4.3.2               </sub><br>
<sub>[121] shiny_1.7.2               lifecycle_1.0.1           nlme_3.1-160              jsonlite_1.8.0            BiocNeighbors_1.12.0      carData_3.0-5            </sub><br>
<sub>[127] fansi_1.0.3               pillar_1.8.1              lattice_0.20-45           KEGGREST_1.34.0           fastmap_1.1.0             httr_1.4.4               </sub><br>
<sub>[133] survival_3.4-0            glue_1.6.2                png_0.1-7                 bit_4.0.4                 blob_1.2.3                BiocSingular_1.10.0      </sub><br>
<sub>[139] memoise_2.0.1             irlba_2.3.5               future.apply_1.9.0</sub><br>
