---
layout: post
comments: true
tags: dementia R ETL
title: Loading Gene Expression Data
---
I've recently been working with the openly-available [*Aging, Dementia, and TBI Study* data](http://aging.brain-map.org/overview/home) from the Allen Institute for Brain Science as part of a Master's degree thesis project. The goal of the project is to construct models of dementia status using the gene expression, pathological, and medical history data provided. In the spirit of "open science", I'm going to try to do an "open thesis", which means I'll be blogging along as I work on the project in posts with the tag, `dementia`. There's also [a GitHub repository for this project](https://github.com/brilliantFire/Allen-aging-dementia-TBI).

As part of exploring this data, one of the tasks I have set for myself is a [differential gene expression (DGE) analysis](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0190152) to see which genes have different patterns of expression under different conditions (dementia versus no dementia, for example). The scientists at the Allen Institute have already done an analysis like this. (Be sure to [check out their awesome interactive tool](http://aging.brain-map.org/rnaseq/searches?%7B%22search_type%22%3A%22differential%22%2C%22target_features%22%3A%5B%7B%22value%22%3A%2210294%22%7D%5D%2C%22contrast_features%22%3A%5B%7B%22value%22%3A%2210294%22%7D%5D%2C%22target_tumors%22%3A%5B%7B%22value%22%3A%22309335438%22%7D%2C%7B%22value%22%3A%22309335440%22%7D%2C%7B%22value%22%3A%22309335441%22%7D%2C%7B%22value%22%3A%22309335443%22%7D%2C%7B%22value%22%3A%22309335448%22%7D%2C%7B%22value%22%3A%22309335451%22%7D%2C%7B%22value%22%3A%22309335452%22%7D%2C%7B%22value%22%3A%22309335453%22%7D%2C%7B%22value%22%3A%22309335454%22%7D%2C%7B%22value%22%3A%22309335456%22%7D%2C%7B%22value%22%3A%22309335457%22%7D%2C%7B%22value%22%3A%22309335460%22%7D%2C%7B%22value%22%3A%22309335461%22%7D%2C%7B%22value%22%3A%22309335464%22%7D%2C%7B%22value%22%3A%22309335465%22%7D%2C%7B%22value%22%3A%22309335467%22%7D%2C%7B%22value%22%3A%22309335469%22%7D%2C%7B%22value%22%3A%22309335474%22%7D%2C%7B%22value%22%3A%22309335475%22%7D%2C%7B%22value%22%3A%22309335476%22%7D%2C%7B%22value%22%3A%22309335478%22%7D%2C%7B%22value%22%3A%22309335482%22%7D%2C%7B%22value%22%3A%22309335486%22%7D%2C%7B%22value%22%3A%22309335487%22%7D%2C%7B%22value%22%3A%22309335488%22%7D%2C%7B%22value%22%3A%22309335489%22%7D%2C%7B%22value%22%3A%22309335490%22%7D%2C%7B%22value%22%3A%22309335494%22%7D%2C%7B%22value%22%3A%22309335495%22%7D%2C%7B%22value%22%3A%22309335497%22%7D%2C%7B%22value%22%3A%22326765649%22%7D%2C%7B%22value%22%3A%22326765652%22%7D%2C%7B%22value%22%3A%22326765656%22%7D%2C%7B%22value%22%3A%22326765657%22%7D%2C%7B%22value%22%3A%22326765659%22%7D%2C%7B%22value%22%3A%22326765661%22%7D%2C%7B%22value%22%3A%22326765663%22%7D%2C%7B%22value%22%3A%22326765665%22%7D%2C%7B%22value%22%3A%22326765666%22%7D%2C%7B%22value%22%3A%22326765667%22%7D%2C%7B%22value%22%3A%22326765668%22%7D%2C%7B%22value%22%3A%22326765671%22%7D%2C%7B%22value%22%3A%22326765672%22%7D%2C%7B%22value%22%3A%22326765677%22%7D%2C%7B%22value%22%3A%22326765681%22%7D%2C%7B%22value%22%3A%22326765682%22%7D%2C%7B%22value%22%3A%22326765683%22%7D%2C%7B%22value%22%3A%22326765684%22%7D%2C%7B%22value%22%3A%22326765688%22%7D%2C%7B%22value%22%3A%22326765689%22%7D%2C%7B%22value%22%3A%22467056391%22%7D%2C%7B%22value%22%3A%22467056397%22%7D%2C%7B%22value%22%3A%22467056405%22%7D%2C%7B%22value%22%3A%22467056407%22%7D%5D%2C%22contrast_tumors%22%3A%5B%7B%22value%22%3A%22309335439%22%7D%2C%7B%22value%22%3A%22309335442%22%7D%2C%7B%22value%22%3A%22309335444%22%7D%2C%7B%22value%22%3A%22309335445%22%7D%2C%7B%22value%22%3A%22309335446%22%7D%2C%7B%22value%22%3A%22309335447%22%7D%2C%7B%22value%22%3A%22309335449%22%7D%2C%7B%22value%22%3A%22309335450%22%7D%2C%7B%22value%22%3A%22309335455%22%7D%2C%7B%22value%22%3A%22309335458%22%7D%2C%7B%22value%22%3A%22309335459%22%7D%2C%7B%22value%22%3A%22309335462%22%7D%2C%7B%22value%22%3A%22309335463%22%7D%2C%7B%22value%22%3A%22309335466%22%7D%2C%7B%22value%22%3A%22309335468%22%7D%2C%7B%22value%22%3A%22309335470%22%7D%2C%7B%22value%22%3A%22309335471%22%7D%2C%7B%22value%22%3A%22309335477%22%7D%2C%7B%22value%22%3A%22309335479%22%7D%2C%7B%22value%22%3A%22309335480%22%7D%2C%7B%22value%22%3A%22309335481%22%7D%2C%7B%22value%22%3A%22309335484%22%7D%2C%7B%22value%22%3A%22309335485%22%7D%2C%7B%22value%22%3A%22309335491%22%7D%2C%7B%22value%22%3A%22309335492%22%7D%2C%7B%22value%22%3A%22309335493%22%7D%2C%7B%22value%22%3A%22326765648%22%7D%2C%7B%22value%22%3A%22326765650%22%7D%2C%7B%22value%22%3A%22326765651%22%7D%2C%7B%22value%22%3A%22326765653%22%7D%2C%7B%22value%22%3A%22326765654%22%7D%2C%7B%22value%22%3A%22326765655%22%7D%2C%7B%22value%22%3A%22326765658%22%7D%2C%7B%22value%22%3A%22326765660%22%7D%2C%7B%22value%22%3A%22326765662%22%7D%2C%7B%22value%22%3A%22326765669%22%7D%2C%7B%22value%22%3A%22326765673%22%7D%2C%7B%22value%22%3A%22326765674%22%7D%2C%7B%22value%22%3A%22326765675%22%7D%2C%7B%22value%22%3A%22326765676%22%7D%2C%7B%22value%22%3A%22326765678%22%7D%2C%7B%22value%22%3A%22326765680%22%7D%2C%7B%22value%22%3A%22326765685%22%7D%2C%7B%22value%22%3A%22326765686%22%7D%2C%7B%22value%22%3A%22326765687%22%7D%2C%7B%22value%22%3A%22467056406%22%7D%2C%7B%22value%22%3A%22467056408%22%7D%2C%7B%22value%22%3A%22467056409%22%7D%5D%2C%22page_num%22%3A0%2C%22column_grouping%22%3A%22act_demented%2Cnincds_arda_diagnosis%22%7D)).  

There are several R packages available ([edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), others) for performing differential analysis on RNA-seq data, like we have here. [*RNA-seq*](https://en.wikipedia.org/wiki/RNA-Seq) is a means of quantifying the amount of mRNA transcripts produced from the genes in a biological sample.   

In this post, I'll discuss:  
1. A brief data overview; and,    
2. Some R code that loads the gene expression data from URLs.  

There's a [Python version of this load script](https://github.com/brilliantFire/Allen-aging-dementia-TBI/blob/master/ETL/expected_count_TPM_FPKM_data_load.py) in the project repo, too. 

### A Little About the Data
I'll discuss the data and its origins in greater detail in a future post. For now, I'd just like to touch on some major points:  

* The data consist of gene expression levels and neuropathological measurements from 377 brain tissue samples, as well as demographic and medical history information, from 107 participants in the [*Adult Changes in Thought (ACT)*](https://www.kpwashingtonresearch.org/our-research/research-areas/aging-geriatrics/act-study-long-running-study-aging-examines-changes-kaiser-permanente-patients-over-time/) study. This massive longitudinal study aims to track the brain function of a large cohort of individuals as they age.  
* The initial patient sample was n=110 individuals, 55 who had a history of at least one TBI and 55 age-and-sex-matched controls with no history of TBI. Three subjects were excluded because their data failed to meet quality control standards (see [TBI Overview](http://help.brain-map.org/download/attachments/9895983/TBI_Overview.pdf?version=2&modificationDate=1492728684044&api=v2) from the study website for more details).    
* All of the data is [freely available](http://aging.brain-map.org/download/index) from the Allen Institute for Brain Science (big thanks to the Allen Institute for this resource!).  

### Loading Gene Expression Measurements
The first thing I have to do is get the data. What I need is a dataframe where each column contains the raw expression counts from the RNA-seq experiment for one of the 377 brain tissue samples and each row is a gene. I'll be using the `fread()` function from the `data.table` package to access files on [the study data download site](http://aging.brain-map.org/download/index) via their URLs.


```R
library(data.table)
```

The first file I want to put into a dataframe is a .csv that contains URLs for files for each of the 377 samples. 


```R
data_files <- data.frame(fread('http://aging.brain-map.org/data/tbi_data_files.csv'))
head(data_files)
```


<table>
<thead><tr><th scope=col>donor_id</th><th scope=col>donor_name</th><th scope=col>specimen_id</th><th scope=col>specimen_name</th><th scope=col>rna_well</th><th scope=col>rna_integrity_number</th><th scope=col>structure_id</th><th scope=col>structure_acronym</th><th scope=col>structure_name</th><th scope=col>rnaseq_profile_id</th><th scope=col>rnaseq_total_reads</th><th scope=col>rnaseq_percent_reads_aligned_to_mrna</th><th scope=col>rnaseq_percent_reads_aligned_to_ncrna</th><th scope=col>rnaseq_percent_reads_aligned_to_genome_only</th><th scope=col>gene_level_fpkm_file_link</th><th scope=col>anonymized_bam_file_link</th><th scope=col>anonymized_bam_index_file_link</th><th scope=col>bigwig_file_link</th></tr></thead>
<tbody>
	<tr><td>309335438                                 </td><td>H14.09.001                                </td><td>309357595                                 </td><td>H14.09.001.HIP.05                         </td><td>320630832                                 </td><td>7.3                                       </td><td>10294                                     </td><td>HIP                                       </td><td>hippocampus (hippocampal formation)       </td><td>496100314                                 </td><td>32275545                                  </td><td>31.7                                      </td><td>7.11                                      </td><td>48.1                                      </td><td>/api/v2/well_known_file_download/496128434</td><td>/api/v2/well_known_file_download/501023521</td><td>/api/v2/well_known_file_download/501023519</td><td>/api/v2/well_known_file_download/501056759</td></tr>
	<tr><td>309335438                                 </td><td>H14.09.001                                </td><td>309357596                                 </td><td>H14.09.001.PCx.01                         </td><td>320630834                                 </td><td>7.2                                       </td><td>10557                                     </td><td>FWM                                       </td><td>white matter of forebrain                 </td><td>496100278                                 </td><td>32515376                                  </td><td>29.0                                      </td><td>8.21                                      </td><td>49.5                                      </td><td>/api/v2/well_known_file_download/496106975</td><td>/api/v2/well_known_file_download/500938472</td><td>/api/v2/well_known_file_download/500938470</td><td>/api/v2/well_known_file_download/500941271</td></tr>
	<tr><td>309335438                                 </td><td>H14.09.001                                </td><td>309357596                                 </td><td>H14.09.001.PCx.01                         </td><td>320630836                                 </td><td>7.1                                       </td><td>10208                                     </td><td>PCx                                       </td><td>parietal neocortex                        </td><td>496100290                                 </td><td>34426215                                  </td><td>29.1                                      </td><td>6.59                                      </td><td>52.1                                      </td><td>/api/v2/well_known_file_download/496555481</td><td>/api/v2/well_known_file_download/500941225</td><td>/api/v2/well_known_file_download/500941223</td><td>/api/v2/well_known_file_download/500941630</td></tr>
	<tr><td>309335438                                 </td><td>H14.09.001                                </td><td>309357599                                 </td><td>H14.09.001.TCx.01                         </td><td>320630838                                 </td><td>7.3                                       </td><td>10235                                     </td><td>TCx                                       </td><td>temporal neocortex                        </td><td>496100279                                 </td><td>31714711                                  </td><td>31.4                                      </td><td>6.97                                      </td><td>48.7                                      </td><td>/api/v2/well_known_file_download/496106814</td><td>/api/v2/well_known_file_download/500936841</td><td>/api/v2/well_known_file_download/500936839</td><td>/api/v2/well_known_file_download/500941005</td></tr>
	<tr><td>309335439                                 </td><td>H14.09.002                                </td><td>309357603                                 </td><td>H14.09.002.HIP.01                         </td><td>320630842                                 </td><td>6.4                                       </td><td>10294                                     </td><td>HIP                                       </td><td>hippocampus (hippocampal formation)       </td><td>496100281                                 </td><td>33402591                                  </td><td>29.5                                      </td><td>7.21                                      </td><td>50.6                                      </td><td>/api/v2/well_known_file_download/496106950</td><td>/api/v2/well_known_file_download/500938896</td><td>/api/v2/well_known_file_download/500938894</td><td>/api/v2/well_known_file_download/500941300</td></tr>
	<tr><td>309335439                                 </td><td>H14.09.002                                </td><td>309357607                                 </td><td>H14.09.002.PCx.01                         </td><td>320630848                                 </td><td>6.5                                       </td><td>10557                                     </td><td>FWM                                       </td><td>white matter of forebrain                 </td><td>496100284                                 </td><td>31652453                                  </td><td>27.7                                      </td><td>6.66                                      </td><td>53.0                                      </td><td>/api/v2/well_known_file_download/496106772</td><td>/api/v2/well_known_file_download/500940428</td><td>/api/v2/well_known_file_download/500940426</td><td>/api/v2/well_known_file_download/500941557</td></tr>
</tbody>
</table>



The URLs for the gene expression data are in the 'gene_level_fpkm_file_link' column. I'll use the 'rnaseq_profile_id'


```R
# Series of links to gene expression files
data_links <- data_files['gene_level_fpkm_file_link']

# Grab the list of sample IDs
sample_IDs <- data_files['rnaseq_profile_id']
```

Behind each of the links in `data_links` is a file that contains three different measures of gene expression for each of 50,283 genes:  
1. 'expected_count' - this is the raw estimated count of the number of transcripts; 
2. 'TPM' - transcripts per million; this is count normalized by the length of the gene; and,  
3. 'FPKM' - fragments per kilobase million; normalized a different way than TPM. 

Here's a peak at the file for the first sample:


```R
sample000_url <- paste('http://aging.brain-map.org', toString(data_links[1,1]),sep='')
start <- Sys.time()    # timing how long it takes to get the file
sample000 <- data.frame(fread(sample000_url))
stop <- Sys.time()
head(sample000)
```


<table>
<thead><tr><th scope=col>gene_id</th><th scope=col>transcript_id.s.</th><th scope=col>length</th><th scope=col>effective_length</th><th scope=col>expected_count</th><th scope=col>TPM</th><th scope=col>FPKM</th></tr></thead>
<tbody>
	<tr><td>        1                                                                                                                                                     </td><td>NM_130786.3                                                                                                                                                   </td><td>1766.00                                                                                                                                                       </td><td>1726.13                                                                                                                                                       </td><td>  16.00                                                                                                                                                       </td><td>  0.97                                                                                                                                                        </td><td>  0.68                                                                                                                                                        </td></tr>
	<tr><td>       10                                                                                                                                                     </td><td>NM_000015.2,XM_011544358.1                                                                                                                                    </td><td>1317.00                                                                                                                                                       </td><td>1277.13                                                                                                                                                       </td><td>   1.00                                                                                                                                                       </td><td>  0.08                                                                                                                                                        </td><td>  0.06                                                                                                                                                        </td></tr>
	<tr><td>      100                                                                                                                                                     </td><td>NM_000022.2,XM_005260236.2,XM_011528478.1,XM_011528479.1,XR_244129.1                                                                                          </td><td>1493.40                                                                                                                                                       </td><td>1453.54                                                                                                                                                       </td><td>  19.00                                                                                                                                                       </td><td>  1.37                                                                                                                                                        </td><td>  0.97                                                                                                                                                        </td></tr>
	<tr><td>     1000                                                                                                                                                     </td><td>NM_001792.3,XM_005258181.2,XM_005258182.1,XM_011525787.1,XM_011525788.1                                                                                       </td><td>3838.31                                                                                                                                                       </td><td>3798.45                                                                                                                                                       </td><td>1341.00                                                                                                                                                       </td><td> 36.96                                                                                                                                                        </td><td> 25.99                                                                                                                                                        </td></tr>
	<tr><td>    10000                                                                                                                                                     </td><td>NM_001206729.1,NM_005465.4,NM_181690.2,XM_005272994.3,XM_005272995.2,XM_005272997.3,XM_006711726.2,XM_011544011.1,XM_011544012.1,XM_011544013.1,XM_011544014.1</td><td>2773.40                                                                                                                                                       </td><td>2733.53                                                                                                                                                       </td><td>4572.14                                                                                                                                                       </td><td>175.26                                                                                                                                                        </td><td>123.24                                                                                                                                                        </td></tr>
	<tr><td>100009613                                                                                                                                                     </td><td>NR_103835.1                                                                                                                                                   </td><td> 804.00                                                                                                                                                       </td><td> 764.14                                                                                                                                                       </td><td>   0.00                                                                                                                                                       </td><td>  0.00                                                                                                                                                        </td><td>  0.00                                                                                                                                                        </td></tr>
</tbody>
</table>



The 'gene_id' column provides the [Entrez Gene database number](https://www.ncbi.nlm.nih.gov/gene) for the gene. Here are the dimensions of each individual table.


```R
dim(sample000)
```


<ol class=list-inline>
	<li>50283</li>
	<li>7</li>
</ol>



Let's see how long it took to get the file. 


```R
stop-start
```


    Time difference of 8.970883 secs


Not surprisingly, it takes a while to get a file that big over the interwebz (at least with the connection I have here).  Regardless, to get all the data, we'll loop through all 377 URLs in `data_links` and add the sample columns to dataframes. Even though I only need the 'expected_count' column from each file for the differential analysis, I'll make 'TPM' and 'FPKM' dataframes, too.


```R
# Initialize dataframes
TPM <- data.frame(matrix(nrow = 50283))
FPKM <- data.frame(matrix(nrow = 50283))
raw_read_counts <- data.frame(matrix(nrow = 50283))
```

Now, we populate the dataframes with gene expression data and save them to csv.


```R
# start a timer for DF construction
start <- Sys.time()

# loops through sample files and adds data for each sample to a dataframe with the 
# 'rnaseq_profile_id' as header
for (sample in 1:nrow(data_links)){
    # I limited the print updates in this notebook.
    # Take out the 'if' to print an update with each of the 377 samples
    if (sample >= 370){
        print(paste('Sample #', sample, 'being added.'))
        flush.console()     # Forces print in loop
        }
    url <- paste('http://aging.brain-map.org', toString(data_links[sample,1]),sep='')
    sample_data <- data.frame(fread(url))
    # Add 'gene_id' from the first sample to each dataframe (it's the same for each file)
    if (sample == 1){
        TPM['gene_id'] <- sample_data['gene_id']
        FPKM['gene_id'] <- sample_data['gene_id']
        raw_read_counts['gene_id'] <- sample_data['gene_id']
        }
    TPM[toString(sample_IDs[sample,1])] <- sample_data['TPM']
    FPKM[toString(sample_IDs[sample,1])] <- sample_data['FPKM']
    raw_read_counts[toString(sample_IDs[sample,1])] <- sample_data['expected_count']
    }

# stop timer, calculate running time
stop <- Sys.time()
duration <- stop-start

# drop extra column
TPM <- TPM[,-1]
FPKM <- FPKM[,-1]
raw_read_counts <- raw_read_counts[,-1]

# Write to csv for later
write.csv(raw_read_counts,'raw_read_counts.csv')
write.csv(TPM,'TPM.csv')
write.csv(FPKM,'FPKM.csv')

print(paste('All samples loaded & saved! Loop duration:', round(duration, 2), 'minutes'))
```

    [1] "Sample # 370 being added."
    [1] "Sample # 371 being added."
    [1] "Sample # 372 being added."
    [1] "Sample # 373 being added."
    [1] "Sample # 374 being added."
    [1] "Sample # 375 being added."
    [1] "Sample # 376 being added."
    [1] "Sample # 377 being added."
    [1] "All samples loaded & saved! Loop duration: 46.26 minutes"
    

Here's a look at the raw_read_counts dataframe.


```R
head(raw_read_counts)
```


<table>
<thead><tr><th scope=col>gene_id</th><th scope=col>496100314</th><th scope=col>496100278</th><th scope=col>496100290</th><th scope=col>496100279</th><th scope=col>496100281</th><th scope=col>496100284</th><th scope=col>496100283</th><th scope=col>496100285</th><th scope=col>496100288</th><th scope=col>...</th><th scope=col>496100648</th><th scope=col>496100643</th><th scope=col>496100666</th><th scope=col>496100651</th><th scope=col>496100652</th><th scope=col>496100650</th><th scope=col>496100664</th><th scope=col>496100657</th><th scope=col>496100653</th><th scope=col>496100646</th></tr></thead>
<tbody>
	<tr><td>        1</td><td>  16.00  </td><td>  34.17  </td><td>  13.15  </td><td>  13.26  </td><td>  12.63  </td><td>  14.73  </td><td>  16.32  </td><td>   8.00  </td><td>  13.29  </td><td>...      </td><td>  11.17  </td><td>   8.17  </td><td>   9.03  </td><td>  10.26  </td><td>   8.24  </td><td>   8.00  </td><td>  16.37  </td><td>  12.33  </td><td>   4.08  </td><td>  10.81  </td></tr>
	<tr><td>       10</td><td>   1.00  </td><td>   0.00  </td><td>   0.00  </td><td>   1.00  </td><td>   7.00  </td><td>   1.00  </td><td>   3.00  </td><td>   0.00  </td><td>   0.00  </td><td>...      </td><td>   0.00  </td><td>   2.00  </td><td>   2.00  </td><td>   3.00  </td><td>   0.00  </td><td>   1.00  </td><td>   2.00  </td><td>   1.00  </td><td>   1.00  </td><td>   3.00  </td></tr>
	<tr><td>      100</td><td>  19.00  </td><td>  76.00  </td><td>  16.00  </td><td>  22.00  </td><td>  25.00  </td><td>  35.00  </td><td>  23.00  </td><td>  30.00  </td><td>  41.00  </td><td>...      </td><td>  46.00  </td><td>  30.34  </td><td>  32.00  </td><td>  81.00  </td><td>  33.00  </td><td>  45.00  </td><td>  39.00  </td><td>  22.00  </td><td>  20.00  </td><td>  21.00  </td></tr>
	<tr><td>     1000</td><td>1341.00  </td><td> 571.00  </td><td>1232.00  </td><td>1340.06  </td><td>1293.00  </td><td>1167.00  </td><td>1001.00  </td><td>1352.00  </td><td>1157.00  </td><td>...      </td><td>1418.00  </td><td>1160.00  </td><td>1380.00  </td><td>1047.00  </td><td>1302.00  </td><td>1652.00  </td><td>1510.00  </td><td>1285.00  </td><td>1570.00  </td><td>1799.00  </td></tr>
	<tr><td>    10000</td><td>4572.14  </td><td>2965.36  </td><td>5881.98  </td><td>5978.63  </td><td>4225.08  </td><td>4781.67  </td><td>4512.38  </td><td>6087.73  </td><td>5130.29  </td><td>...      </td><td>5618.93  </td><td>4326.01  </td><td>4420.92  </td><td>4664.40  </td><td>3952.09  </td><td>6153.32  </td><td>5903.16  </td><td>5362.97  </td><td>6598.29  </td><td>6838.52  </td></tr>
	<tr><td>100009613</td><td>   0.00  </td><td>   0.00  </td><td>   0.00  </td><td>   1.00  </td><td>   0.00  </td><td>   0.00  </td><td>   0.00  </td><td>   1.00  </td><td>   0.00  </td><td>...      </td><td>   1.00  </td><td>   0.00  </td><td>   0.00  </td><td>   0.00  </td><td>   0.00  </td><td>   1.00  </td><td>   0.00  </td><td>   0.00  </td><td>   0.00  </td><td>   1.00  </td></tr>
</tbody>
</table>



Looks good! 

#### **Next up:** Differential gene expression analysis
