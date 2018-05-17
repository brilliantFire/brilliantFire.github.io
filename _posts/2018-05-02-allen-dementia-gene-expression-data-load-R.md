---
layout: post
comments: true
tags: dementia R-stats ETL
title: Loading Gene Expression Data
---
I've recently been working with the openly-available [*Aging, Dementia, and TBI Study* data](http://aging.brain-map.org/overview/home) from the Allen Institute for Brain Science as part of a Master's degree thesis project. The goal of the project is to construct models of dementia status using the gene expression, pathological, and medical history data provided. In the spirit of "open science", I'm going to try to do an "open thesis", which means I'll be blogging along as I work on the project in posts with the tag, `dementia`. There's also [a GitHub repository for this project](https://github.com/brilliantFire/Allen-aging-dementia-TBI).

As part of exploring this data, one of the tasks I have set for myself is a [differential gene expression (DGE) analysis](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0190152) to see which genes have different patterns of expression under different conditions (dementia versus no dementia, for example). The scientists at the Allen Institute have already done an analysis like this. (Be sure to [check out their awesome interactive tool](http://aging.brain-map.org/rnaseq/searches?%7B%22search_type%22%3A%22differential%22%2C%22target_features%22%3A%5B%7B%22value%22%3A%2210294%22%7D%5D%2C%22contrast_features%22%3A%5B%7B%22value%22%3A%2210294%22%7D%5D%2C%22target_tumors%22%3A%5B%7B%22value%22%3A%22309335438%22%7D%2C%7B%22value%22%3A%22309335440%22%7D%2C%7B%22value%22%3A%22309335441%22%7D%2C%7B%22value%22%3A%22309335443%22%7D%2C%7B%22value%22%3A%22309335448%22%7D%2C%7B%22value%22%3A%22309335451%22%7D%2C%7B%22value%22%3A%22309335452%22%7D%2C%7B%22value%22%3A%22309335453%22%7D%2C%7B%22value%22%3A%22309335454%22%7D%2C%7B%22value%22%3A%22309335456%22%7D%2C%7B%22value%22%3A%22309335457%22%7D%2C%7B%22value%22%3A%22309335460%22%7D%2C%7B%22value%22%3A%22309335461%22%7D%2C%7B%22value%22%3A%22309335464%22%7D%2C%7B%22value%22%3A%22309335465%22%7D%2C%7B%22value%22%3A%22309335467%22%7D%2C%7B%22value%22%3A%22309335469%22%7D%2C%7B%22value%22%3A%22309335474%22%7D%2C%7B%22value%22%3A%22309335475%22%7D%2C%7B%22value%22%3A%22309335476%22%7D%2C%7B%22value%22%3A%22309335478%22%7D%2C%7B%22value%22%3A%22309335482%22%7D%2C%7B%22value%22%3A%22309335486%22%7D%2C%7B%22value%22%3A%22309335487%22%7D%2C%7B%22value%22%3A%22309335488%22%7D%2C%7B%22value%22%3A%22309335489%22%7D%2C%7B%22value%22%3A%22309335490%22%7D%2C%7B%22value%22%3A%22309335494%22%7D%2C%7B%22value%22%3A%22309335495%22%7D%2C%7B%22value%22%3A%22309335497%22%7D%2C%7B%22value%22%3A%22326765649%22%7D%2C%7B%22value%22%3A%22326765652%22%7D%2C%7B%22value%22%3A%22326765656%22%7D%2C%7B%22value%22%3A%22326765657%22%7D%2C%7B%22value%22%3A%22326765659%22%7D%2C%7B%22value%22%3A%22326765661%22%7D%2C%7B%22value%22%3A%22326765663%22%7D%2C%7B%22value%22%3A%22326765665%22%7D%2C%7B%22value%22%3A%22326765666%22%7D%2C%7B%22value%22%3A%22326765667%22%7D%2C%7B%22value%22%3A%22326765668%22%7D%2C%7B%22value%22%3A%22326765671%22%7D%2C%7B%22value%22%3A%22326765672%22%7D%2C%7B%22value%22%3A%22326765677%22%7D%2C%7B%22value%22%3A%22326765681%22%7D%2C%7B%22value%22%3A%22326765682%22%7D%2C%7B%22value%22%3A%22326765683%22%7D%2C%7B%22value%22%3A%22326765684%22%7D%2C%7B%22value%22%3A%22326765688%22%7D%2C%7B%22value%22%3A%22326765689%22%7D%2C%7B%22value%22%3A%22467056391%22%7D%2C%7B%22value%22%3A%22467056397%22%7D%2C%7B%22value%22%3A%22467056405%22%7D%2C%7B%22value%22%3A%22467056407%22%7D%5D%2C%22contrast_tumors%22%3A%5B%7B%22value%22%3A%22309335439%22%7D%2C%7B%22value%22%3A%22309335442%22%7D%2C%7B%22value%22%3A%22309335444%22%7D%2C%7B%22value%22%3A%22309335445%22%7D%2C%7B%22value%22%3A%22309335446%22%7D%2C%7B%22value%22%3A%22309335447%22%7D%2C%7B%22value%22%3A%22309335449%22%7D%2C%7B%22value%22%3A%22309335450%22%7D%2C%7B%22value%22%3A%22309335455%22%7D%2C%7B%22value%22%3A%22309335458%22%7D%2C%7B%22value%22%3A%22309335459%22%7D%2C%7B%22value%22%3A%22309335462%22%7D%2C%7B%22value%22%3A%22309335463%22%7D%2C%7B%22value%22%3A%22309335466%22%7D%2C%7B%22value%22%3A%22309335468%22%7D%2C%7B%22value%22%3A%22309335470%22%7D%2C%7B%22value%22%3A%22309335471%22%7D%2C%7B%22value%22%3A%22309335477%22%7D%2C%7B%22value%22%3A%22309335479%22%7D%2C%7B%22value%22%3A%22309335480%22%7D%2C%7B%22value%22%3A%22309335481%22%7D%2C%7B%22value%22%3A%22309335484%22%7D%2C%7B%22value%22%3A%22309335485%22%7D%2C%7B%22value%22%3A%22309335491%22%7D%2C%7B%22value%22%3A%22309335492%22%7D%2C%7B%22value%22%3A%22309335493%22%7D%2C%7B%22value%22%3A%22326765648%22%7D%2C%7B%22value%22%3A%22326765650%22%7D%2C%7B%22value%22%3A%22326765651%22%7D%2C%7B%22value%22%3A%22326765653%22%7D%2C%7B%22value%22%3A%22326765654%22%7D%2C%7B%22value%22%3A%22326765655%22%7D%2C%7B%22value%22%3A%22326765658%22%7D%2C%7B%22value%22%3A%22326765660%22%7D%2C%7B%22value%22%3A%22326765662%22%7D%2C%7B%22value%22%3A%22326765669%22%7D%2C%7B%22value%22%3A%22326765673%22%7D%2C%7B%22value%22%3A%22326765674%22%7D%2C%7B%22value%22%3A%22326765675%22%7D%2C%7B%22value%22%3A%22326765676%22%7D%2C%7B%22value%22%3A%22326765678%22%7D%2C%7B%22value%22%3A%22326765680%22%7D%2C%7B%22value%22%3A%22326765685%22%7D%2C%7B%22value%22%3A%22326765686%22%7D%2C%7B%22value%22%3A%22326765687%22%7D%2C%7B%22value%22%3A%22467056406%22%7D%2C%7B%22value%22%3A%22467056408%22%7D%2C%7B%22value%22%3A%22467056409%22%7D%5D%2C%22page_num%22%3A0%2C%22column_grouping%22%3A%22act_demented%2Cnincds_arda_diagnosis%22%7D)).  

There are several R packages available ([edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html), [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), others) for performing differential analysis on RNA-seq data, like we have here. [*RNA-seq*](https://en.wikipedia.org/wiki/RNA-Seq) is a means of quantifying the amount of mRNA transcripts produced from the genes in a biological sample.  More on that later.

For now, in this post, I'll discuss:  
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
names(data_files)
```

![`data_files` dataframe column names]({{ site.baseurl }}/images/2018-05-02-image-01.png)


The URLs for the gene expression data are in the 'gene_level_fpkm_file_link' column. I'll use the 'rnaseq_profile_id' to label the columns/samples.


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
names(sample000)
```
![RNA-seq sample dataframe columns]({{ site.baseurl }}/images/2018-05-02-image-02.png)

The 'gene_id' column provides the [Entrez Gene database number](https://www.ncbi.nlm.nih.gov/gene) for the gene. Here are the dimensions of each individual table.

```R
dim(sample000)
```
![size of each sample file]({{ site.baseurl }}/images/2018-05-02-image-03.png)

Let's see how long it took to get the file.

```R
stop-start
```
![time for one sample file to load]({{ site.baseurl }}/images/2018-05-02-image-04.png)


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
![running output]({{ site.baseurl }}/images/2018-05-02-image-05.png)

Here's a look at the raw_read_counts dataframe.

```R
head(raw_read_counts)
```
![count matrix]({{ site.baseurl }}/images/2018-05-02-image-06.png)

Looks good!

#### **Next up:** Differential gene expression analysis
