# TMB
The tumor microbiome analysis pipeline

Taken from the supplementary methods section of the paper (citation needed):

**16S data analysis.**  
Relative abundances were converted to read counts by multiplying by the total number of reads. Sample reads were normalized within each sequencing library, by a factor representing the ratio of the average number of reads/sample in the specific library to the overall average across libraries. Samples with fewer than 1000 normalized reads (including negative controls) and species with relative abundances of less than 10-4 were discounted from further analysis. 
We then apply a series of six filters to detect and remove contamination:
**Filter 1**: We removed species that had a high prevalence in negative controls; 113 species that appeared in over 7.5% (this number represents the inflection point of the bimodal distribution between most species that appear in very few controls, to species that exhibit a higher range of prevalences (Fig. S9)) of control samples (n=636) were assumed to come from lab handling procedures, biological reagents and kits and were removed from further analysis. 139 species that appeared in over 7.5% of empty paraffin control samples that had more than 1,000 normalized reads (n=165) were assumed to come from global paraffin processing or represent global hospital environment contamination (Table S2). 84 of them overlapped with the previous list so only an additional 55 and were removed from further analysis. Note that overall 168 species were removed as general contaminants (113+55) while in figure 3B it looks like only 167 species were removed as general contaminants (9,190 - 9,023). This difference stems from the fact that one general contaminant that was found in paraffin did not appear at all in our tumor/normal samples and is thus not included in the 9,190 species presented on figure 3B. Many of these general contaminant species corroborate with previously reported lists of contaminating species (Table S2, column 2395) (21, 23). Next, the relative abundances and normalized reads values of higher level taxonomies such as genera and families were calculated based on the sums of those measures in their lower taxonomic members.  
**Filter 2-4**: In order to remove the less prevalent contaminations (that could originate from rarer contamination events that occur during processing or cross-contaminations between samples) we compared taxa prevalences in samples to their prevalences in controls. We took into account that different contamination spectrums can vary across processing batches and time and that different tissue types may have different microbial biomass levels, efficiencies of DNA extraction and PCR amplification so we shifted to a 'per-condition' (a combined tissue+tumor/NAT/normal status) contamination analysis. Therefore, our comparisons between taxon prevalence in samples and controls were done on a per-condition, per-batch basis, comparing samples from the same condition to the controls that were processed alongside them. This was done for the three different processing batch types; DNA extraction batch (filter 2), PCR amplification batch (filter 3), and sequencing library batch (filter 4). We employed the non-parametric, exact, one-tailed binomial test on each taxon’s prevalence, x=number of successes in the sample and  n=sample size to the background p= taxa's prevalence in the relevant batch controls. Only taxa that passed a p-value cutoff of 0.05 in all batch comparisons were considered for the center filter (filter 6). Specifically for the colon tumor and NAT conditions, only one center was sampled so an FDR<0.2 was applied to the p-values of the sequencing library run filter (filter 4) and the center filter was not applied.
**Filter 5**:   We compared prevalence of taxa in all samples per-condition to their prevalence in a set of empty paraffin controls originating from the same centers as the samples, requiring a one-tailed binomial p-value<=0.05 to pass.
**Filter 6**: To account for a center-specific batch effect, we compared prevalence of taxa in samples, again per-condition, this time per-medical-center to the same-center specific controls (DNA extraction, PCR, and library controls that ran with these specific center samples).  This was performed via binomial test (as described above) on each center per-condition. Only taxa with p-value<=0.05  in at least two centers, and a Fisher's combined p-value across all centers, FDR corrected < 0.20, passed this filter.


Requirements:
-------------
Python 2.7 interpreter (python 3+ will fail on some new syntax issues).
Python libraries:
-	pandas
-	numpy
-	scipy
-	os
-	itertools
-	matplotlib.pyplot
-	mne.stats
-	pickle

Files description
-----------------
Infrastructure files

- taxonomy.py	Contains some basic taxonomy definitions, used primarily for filtering for species

Core pipeline files

- analysis_pipeline_16S.py	This is the heart of the pipeline, contains 3 classes that take care of different stages of the pipeline from loading the data into data structures (dataLoader) to hit calling (hitCaller) to controlling the flow between these two classes and logging/plotting different information along the way (analysisPipeline)

- metadata_class.py	This is an infrastructural file that handles the loading, curating and mapping of the clinical metadata to the samples.

Example
---------
Example files are also provided:
- metadata.xlsx - tabular row data about each sample (each column is a sample). There are some required rows of information in order for the pipeline run, e.g. extraction batch numbers which are used by the batch filters.
- reconstruction_files - these are the files produced by the 5R analysis tool (gituhub, NoamShental/5R). They represent the relative read abundances for all species detected in the samples by the 5R method. These files will need to be unzipped and untarred. 

The code is ready to run the example if the example files are stored in the same directory as the code is.
To set different directory paths for the example files, or for your own reconstruction files, please change the parameters you sent to the load_data() function of the dataLoader class object (see the main() function in analysisPipeline for an example).

*Note: 
When using reconstruction files produced by the 5R, a mild tweaking of the file is required in order to input into this analysis pipeline as follows:

The reconstruction files need to be:
1)	Unzipped, and untarred
2)	The first line contains the number of reads per sample, this needs to be cut out of the file into a separate file that I usually name as the <reconstruction filename>_first_line. Then you need to delete the ‘Total number of reads’ prefix and all trailing tabs.
a.	This can be done using cat <reconstruction file> | head -1 > <reconstruction file>_first_line
b.	Then delete the ‘Total number of reads’ prefix and all trailing tabs
c.	In the reconstruction file remove the first line.
Now you’re ready to use the infrastructure for reading these files into python datastructures, hit calling (deciding what’s noise and what’s signal) 


Run
-------------
Run the main() function in Analysis_Pipeline_16S_TMB.py

Contact
-------------
For questions please email:
Ilana Livyatan: ilana.livyatan@weizmann.ac.il
Ravid Straussman: ravidst@weizmann.ac.il
