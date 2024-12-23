# Retrovirus-Analyzer
Program made for analyzing several sequences of a virus of interest and finding similarity between them on a specific gene. In this case, the example and the analysis report, was made with retrovirus in mind. 
This is the repository created for the second project of the subject Programming fundamentals for biological sciences. National University of Colombia (UNAL).

## The steps for gathering the sequences correctly from genbank are easy and as follows:
 1) Head to NCBI through the link https://www.ncbi.nlm.nih.gov/
 2) Type on the search bar "<ins>_human immunodeficiency virus complete genome_</ins>" or something similar
![image](https://github.com/user-attachments/assets/b40673a8-e757-423a-a66d-8d45f2cb193f)
 3) Click on the one you prefer the most and head over to the "Send to" option, click on coding sequences and download
![image](https://github.com/user-attachments/assets/8d96d763-236b-485d-8e76-f2106630c043)

## IMPORTANT NOTE: It's very important to notice that the programm as it is will only work from the terminal/console because of the sys.argv codes. The parameters must be inserted in the following order:

 1) file_name_1, sequence downloaded in cds format from Genbank
 2) file_name_2, sequence downloaded in cds format from Genbank
 3) file_name_3, sequence downloaded in cds format from Genbank
 4) file_name_4, sequence downloaded in cds format from Genbank
 5) file_name_5, sequence downloaded in cds format from Genbank
 6) file_name_6, sequence downloaded in cds format from Genbank
 7) number_of_sequences, in this case, the previous 6
 8) gene_name, the exact gene name as in the sequence; for example: [gene=gag]

### How it should look:
####  python .\retrovirus_analysis_plotter.py  .\sequence.txt .\sequence (1).txt .\sequence (2).txt .\sequence (3).txt .\sequence (4).txt .\sequence (5).txt 6 [gene=gag]

#### From line 131-140 there is a commented section that can be used for local running on the text editor or IDLE, where the information can be entered manually.

### Also, the graphics must be personalized based on the average length of the sequences entered. Specifically, in the lines: 251, 295, 308

## Results: 
 1) The output from the code in the terminal is:
![image](https://github.com/user-attachments/assets/ae34624d-611c-4ff8-b3ba-b658749b8caa)
 2) The code generates each of the graphics in jpg format
 3) The analysis report is created with pylatex and pdflatex tools and generated in pdf format



