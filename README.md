# FAST - FindAmpliconSTructure

FAST is a tool for analyzing the structure and characteristics of amplicons, based on Copy Number Variations and Structural Variations data.

## Prerequisites
Python 2.7.x
<br>
Packages:
<br>
- numpy >= 1.8.2
- pandas == 0.17.0<br>
- xlsxwriter >= 0.6.6<br>
- xlwt >= 0.7.5<br>
- xlrd<br>

## Usage
`python Fast.py [-h] --config CONFIG [--excel] [--remove_cent]`
<br>
Parameters:
- --help, -h				show this help message and exit
-	--config CONFIG, -c CONFIG	Configuration file
-	--excel, -e           			generated colored excel file
-	--remove_unmappble, -r     	remove unmappable areas, specifically centromeres and telomeres data

## Inputs
### Input files
FAST requires the following input data:
1. Structural Variations (SVs) in the genome of the analyzed, as generated by BreakDancerMax.
2. Copy Number Variations (CNVs) in the genome of the analyzed sample, as generated by Control-FREEC.
The current version of FAST relies on the file formats of these tools. Thus, it is recommended to use these tools for generating inputs to FAST. However, other tools can be used as well as long as the file format is kept.

### Configuration File
The configuration file consists of two groups; the first is [PARAMS], in which certain parameters and constants are set. The user may use the default values, or update them. The second group is [FILES], in which the input files are supplied.

#### Configuration File [PARAMS] parameters:
-	WINDOW = 15000
Window size used for digitization.
-	WINDOW_MULT = 3
The tool calculates threshold as WINDOW_MULT * WINDOW, and uses it for classifying the breakpoints to the corresponding segments.
-	CNV_MINIMAL_CN = 6
Minimal copy number value for considering an area as amplified.
-	SV_MINIMAL_SCORE = 90
Breakpoints with a lower score will be removed from the analysis.
-	SV_MINIMAL_NUM_READS = 2
Breakpoints with a lower number of supporting reads will be removed from the analysis.
-	SV_MINIMAL_DELETIONS_REMOVE = 1000
Deletions which are shorter than this value will be removed from the analysis.

#### Configuration File [FILES] parameters:
-	UNMAPPBLE_INP
A comma delimited file which contains centromeres and telomeres info. Assumed to contain the following fields:
chr
chromStart
chromEnd
type (centromere / telomere)
-	SAMPLE_ID
-	FAST_OUTPUT_DIR
-	CNV_INP
Copy Number Variations (CNVs) file, generated by Control-FREEC
-	SV_INP
Structural Variations (SVs) file, generated by BreakDancerMax.

## Output
FAST output is stored in the configured output directory (FAST_OUTPUT_DIR). 

The tool outputs the following files:
-	<SAMPLE_ID>_cnv_FAST.csv
<br>CNV file after Fast processing
-	<SAMPLE_ID>_sv_FAST.csv
<br>SV file after Fast processing 
-	<SAMPLE_ID>_amplicons_FAST.csv
<br>Amplicon Analysis
-	<SAMPLE_ID>_segments_FAST.csv
<br>Segment Analysis
-	<SAMPLE_ID>_focused_genes_FAST.csv
<br>Analysis of the focused genes loci
-	<SAMPLE_ID>_FAST.xlsx
<br>File summarizing the data of all above files. In this file colors are assigned to cells, such that each amplicon has its own color.  This makes it easier for the human eye to see the relations between the different records.

## Running FAST on test data
Test data of HCC1954 is supplied under `example\data\HCC1954`
<br>
Run by:  `python Fast.py -c example\data\HCC1954\HCC1954_config.txt -e –r`
