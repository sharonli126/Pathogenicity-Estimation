# Description
This software is a pathogenicity socring system for genetic variants in patient, particularly autoimmune disease and kidney disease patients.
It aims to enhance the accuracy of ranking potential pathogenic variants, particularly missense variants, in kidney and autoimmune disease patients.

# How to run the software
⚠️ IMPORTANT: The required data file is not included in the GitHub repository due to its size.

🔹 Step 1: Download the Data File
    
    👉 Download the required data file [here](https://drive.google.com/file/d/1Lo4_9DeWrOemetByDw0wGuHLoQRTEcPM/view?usp=sharing)

🔹 Step 2: Add the Data File to the Project
    
    After downloading, move the file into the cloned repository's directory.

🔹 Step 3: Open and Run the Jupyter Notebook
    
    Open main.ipynb in Jupyter Notebook.
    
    Follow the steps in the notebook to:
    
    Run the program with test data OR Replace the test data with patient data for analysis.

# How to find potential pathogenic variant
Sort the table by total score after you obtain the ouput table, the higher score it is, the more likely it is pathogenic variant. It was validated to have all targeted candidates within top 3%.

# Attributes description 
These are attributes you will obtain in final dataset from the output of the software.

There will be 40 columns, * are annotation added to the patient dataset after running the program

- INDEX: row number

- genomic_nomencalture: the genetic variant id with its location (chr:location original base > mutated base)

- Consequence: type of mutation 

- IMPACT: based on the type of mutation, reflecting the level of impact of particular to gene function

- Amino acid: change of amino acid (orignal/ mutated)

- Codons: the original codon and the mutated codon

- SYMBOL: human gene symbol -> which gene the genetic variant is likely to affect

- AF: minor allele frequency from Ensembl

- gnomAD3_AF: allele frequency from Gnomad version 3

- gnomAD_AF: allele frequency from Gnomad version 1, is used when gnomAD3_AF have invalid values, like "."

- CADD_Phred: A combined measurement that assesses the deleteriousness of genetic variants derived from a machine learning model.

- SpliceAI_pred_DS_AG: SpliceAI predicted effect on splicing. Delta score for acceptor gain

- SpliceAI_pred_DS_AL: 'SpliceAI predicted effect on splicing. Delta score for acceptor loss

- SpliceAI_pred_DS_DG: SpliceAI predicted effect on splicing. Delta score for donor gain

- SpliceAI_pred_DS_DL: SpliceAI predicted effect on splicing. Delta score for donor loss 

- revel*: REVEL score; matched by genetic variant; a measurement built by random forest; reflecting pathogenecity of the missense variant

- term*: gene otonlogy term which identify variants involved in production or regulatory functions of protein such as interferon, cytokines which are potentially deleterious when mutated, it is for later filter, instead of the scoring system.

- mouse_symbol*: mouse gene sysmbol'

- mgi_note*: phenotype of mutation at a particular gene

-  LOEUF_upper*: matched by gene; to assess the intolerance of genes to loss-of-function (LoF) variants, representing the upper confidence bound of the observed/expected (o/e) ratio for LoF mutations.

- lof_z_4_1*: macthed by gene; reflecting the distribution of o/e of loss of function variant but with z score 

- mis_oe_4_1*: matched by gene; observed missense mutation/ expected mutation -> the higher the more likely it is tolerant to missense mutation; 

- mis_z_4_1*: macthed by gene; reflecting the distribution of o/e of missense variant but with z score 

- syn_oe_4_1*: matched by gene; observed synonymous mutation/ expected mutation -> the higher the more likely it is tolerant to synonymous mutation

- syn_z_4_1*:macthed by gene; reflecting the distribution of o/e of synonymous variant but with z score 

- PolyPhe: showing the pathogenecity prediction by PolyPhe; for comparison of genetic score

- ClinVar_CLNSIG: showing the pathogenecity prediction by ClinVar; for comparison of genetic score

- am_pathogenicity*: Alpha missense score from hg 38, it considers protein sequence conservation and structural context of missense variants. Since mutations occurring in highly conserved sequences and affecting protein structure, such as folding and binding sites, are likely to be pathogenic, it is a crucial component for the estimation.

- distinct_patientDatabaseId_count*: The number of patient in the APF database with the particular genetic variants

- variant_frequency*: allele frequency calculated from the APF database; it is served as a attribute to identify false positive in two allele frequency annotation by distinct_patientDatabaseId_count/ total of patient, which is 297

- consequence_score*: calculated based on IMPACT

- splice_score*: calculated by all SpiceAI; only applied for splice variant in the scoring system

- maf_score*: calculated by all the allele frequency annotation and  variant_frequency from APF database

- caddscore*: calculated by CADD PHRED  

- am_score*: calculated by am_pathogenecity 

- mgiscore*: calculted by mgi_notes, when there is keyword -> 1 score 

- revelscore*: calculated by revel  

- genoterm_presence*: to identify whether there is target genoterm from term; presence = 1, absence = 0

- total score*: obtained from sum of score/ maximum score of each variant type (missense, splice, and others)

