# TMB
The tumor microbiome analysis pipeline

Requirements:
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

Infrastructure files

- taxonomy.py	Contains some basic taxonomy definitions, used primarily for filtering for species

Core pipeline files

- analysis_pipeline_16S.py	This is the heart of the pipeline, contains 3 classes that take care of different stages of the pipeline from loading the data into data structures (dataLoader) to hit calling (hitCaller) to controlling the flow between these two classes and logging/plotting different information along the way (analysisPipeline)

- metadata_class.py	This is an infrastructural file that handles the loading, curating and mapping of the clinical metadata to the samples.


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
