import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import random
import sys
from scipy.stats import spearmanr
from pylatex import Document, Section, Subsection, Figure, Command
from pylatex.utils import NoEscape


"""
Please, read the README.md and LICENSE files beforehand in my GitHub repository:
https://github.com/monitoxx/Retrovirus-Analyzer

The following code was made by me (https://github.com/monitoxx) , with some help from Gemini Assistant and ChatGpt for 
finding errors and helping develop the LaTeX part of my code.
The comments provide a really self-explanatory experience from now on.
Reach me at my email for any doubts: ljimenezbe@unal.edu.co
"""

'''
IMPORTANT NOTE: It´s very important to notice that the programm as it is will only work from the terminal/console because of the sys.argv codes.
The parameters must be inserted in the following order:

        1. file_name_1, sequences downloaded in cds format from Genbank
        2. file_name_2, sequences downloaded in cds format from Genbank
        3. file_name_3, sequences downloaded in cds format from Genbank
        4. file_name_4, sequences downloaded in cds format from Genbank
        5. file_name_5, sequences downloaded in cds format from Genbank
        6. file_name_6, sequences downloaded in cds format from Genbank
        7. number_of_sequences, in this case, the previous 6
        8. gene_name, the exact gene name as in the sequence; [gene=gag]

How it should look:
'''
#  python .\retrovirus_analysis_plotter.py  .\sequence.txt .\sequence (1).txt .\sequence (2).txt .\sequence (3).txt .\sequence (4).txt .\sequence (5).txt 6 [gene=gag]
'''
From line 131-140 there is a commented section that can be used for local running on the text editor or IDLE, where the 
information can be entered manually.
'''

####Input from command line on the terminal####
# -----You can modify the following depending on how many sequences you want to analyze

file_name_1 = sys.argv[1]
file_name_2 = sys.argv[2]
file_name_3 = sys.argv[3]
file_name_4 = sys.argv[4]
file_name_5 = sys.argv[5]
file_name_6 = sys.argv[6]
file_names = [
    file_name_1,
    file_name_2,
    file_name_3,
    file_name_4,
    file_name_5,
    file_name_6
]
number_of_sequences = int(sys.argv[7])
gene_name = sys.argv[8]


def gc_content(sequence):
    sequence = sequence.rstrip("\n")
    sequence = sequence.upper()
    G_count = sequence.count("G")
    C_count = sequence.count("C")
    GC_content = G_count + C_count
    GC_content = GC_content / len(sequence)
    GC_content = round(100 * GC_content, 2)
    return float(GC_content), int(len(sequence))
    pass

#The following function identifies our gene of interest and extracts the correspondent sequence
def sequence_reader(file_name, gene_name):
    sequence = ''
    sequence_name = None

    try:
        with open(file_name, 'r') as file:
            include_sequence = False  # Flag for including the sequence
            for line in file:
                #include_sequence = False # Flag for including the sequence
                if line.startswith('>'):
                    # If header contains our gene of interest, activates the flag
                    include_sequence = gene_name in line
                    if include_sequence:
                        # Busca el código GenBank
                        match = re.search(r'\|(.+?)\_', line)
                        if match:
                            sequence_name = match.group(1)  # Extracts code of Genbank
                elif include_sequence:  #Acumulates lines of sequence only when our flag is activated
                    line_ = line.strip()
                    sequence += line.strip()

    except FileNotFoundError:
        print(f"File {file_name} not found.")


    return sequence_name, sequence

def gc_list_generator(sequence, sequence_name):
    position_list = []
    gc_list = []
    gc_list_random = []
    code_list = []

    # Generates a random sequences once
    seq_random_list = list(sequence)
    random.shuffle(seq_random_list)
    seq_random = ''.join(seq_random_list)

    for x in range(0, len(sequence), 70):
        # Fragments of the original and random sequences
        temp = sequence[x:x + 70]
        temp2 = seq_random[x:x + 70]

        # Calculates GC content
        gc_temp, _ = gc_content(temp)
        gc_temp2, _ = gc_content(temp2)

        # Adds dats to lists
        position_list.append(x)
        gc_list.append(gc_temp)
        gc_list_random.append(gc_temp2)
        code_list.append(sequence_name)

    return position_list, gc_list, gc_list_random, code_list

'''
file_names = [
    'sequence.txt',
    'sequence (1).txt',
    'sequence (2).txt',
    'sequence (3).txt',
    'sequence (4).txt',
    'sequence (5).txt'
]
number_of_sequences = 6
gene_name = '[gene=gag]'
'''

def df_creator(file_names, gene_name):
    # Create an empty dataFrame
    df_sequences = pd.DataFrame(columns=['GenBank Code', 'Sequence', 'Length'])
    df_gc_content = pd.DataFrame(columns=['GenBank Code', '% GC Content', 'Length'])
    df_gc_content_lists = pd.DataFrame(columns=['GenBank Code', '% GC Content', '% GC Content Random', 'Position'])

    # Iterate the files and add the data to the dataFrames
    for file_name in file_names:
        sequence_name, sequence = sequence_reader(file_name, gene_name)
        if sequence_name and sequence:  # Only adds up data if they are valid
            var_gc_content, len_sequence = gc_content(sequence)
            position_list, gc_list, gc_list_random, code_list = gc_list_generator(sequence, sequence_name)

            df_sequences = pd.concat([df_sequences, pd.DataFrame(
                {'GenBank Code': [sequence_name],
                 'Sequence': [sequence],
                 'Length': [len_sequence]
                 })], ignore_index=True)
            df_gc_content = pd.concat([df_gc_content, pd.DataFrame(
                {'GenBank Code': [sequence_name],
                 '% GC Content': [var_gc_content],
                 'Length': [len_sequence]
                 })], ignore_index=True)
            df_gc_content_lists = pd.concat([df_gc_content_lists, pd.DataFrame(
                {'GenBank Code': [code_list],
                 '% GC Content': [gc_list],
                 '% GC Content Random': [gc_list_random],
                 'Position': [position_list]
                 })], ignore_index=True)

    #This will later help me for my spearman analysis, so that indexes match
    length = min(df_sequences['Length'])


    # Unfolds lists in lines
    df_gc_content_lists_exploded = df_gc_content_lists.explode(['% GC Content', '% GC Content Random', 'Position', 'GenBank Code'])
    # Makes sure that numeric columns are in the correct format:
    df_gc_content_lists_exploded['% GC Content'] = pd.to_numeric(df_gc_content_lists_exploded['% GC Content'])
    df_gc_content_lists_exploded['% GC Content Random'] = pd.to_numeric(df_gc_content_lists_exploded['% GC Content Random'])
    df_gc_content_lists_exploded['Position'] = pd.to_numeric(df_gc_content_lists_exploded['Position'])


    # Group data by GenBank Code
    grouped_data = df_gc_content_lists_exploded.groupby('GenBank Code')

    # Calculate quartiles for each group
    df_gc_content_lists_exploded['Quartile'] = grouped_data['% GC Content'].transform(pd.qcut, q=4, labels=["Q1", "Q2", "Q3", "Q4"])
    # Repeat for random GC content
    df_gc_content_lists_exploded['Quartile Random'] = grouped_data['% GC Content Random'].transform(pd.qcut, q=4, labels=["Q1", "Q2", "Q3", "Q4"])

    # Creates a pivot matrix for heatmap:
    heatmap_data = df_gc_content_lists_exploded.pivot_table(
        index='GenBank Code', columns='Position', values='% GC Content', aggfunc='mean')

    # Creates a long dataFrame for plotting the comparisons
    df_long = df_gc_content_lists_exploded.melt(
        id_vars=['Position', 'GenBank Code'],
        value_vars=['% GC Content', '% GC Content Random'],
        var_name='Metric', value_name='Value')

    #print(df_gc_content_lists_exploded)

    return df_sequences, df_gc_content, df_gc_content_lists, df_gc_content_lists_exploded, length, heatmap_data, df_long


df_sequences, df_gc_content, df_gc_content_lists, df_gc_content_lists_exploded, length, heatmap_data, df_long = df_creator(file_names, gene_name)

# Verify data
#print(df_sequences)
#print(df_gc_content)
#print(df_gc_content_lists)
#print(df_gc_content_lists_exploded)
gene_name= gene_name.split('=')[1][:-1]
genbank_code_list = ', '.join(df_sequences['GenBank Code'])



#Next is the Spearman analysis between the gc contents of the sequences
for i in range(0, number_of_sequences-1):
    index = round((length/70)-1) #70 Because this is the length of each line of the normal format of the cds.txt files
    try:
        corr, p_value = spearmanr(df_gc_content_lists.loc[i,'% GC Content'][0:index], df_gc_content_lists.loc[i+1,'% GC Content'][0:index])
        print(f"Spearman's correlation of GC Content for the fragmented sequences of "
              f"{df_gc_content_lists.loc[i,'GenBank Code'][0]} and {df_gc_content_lists.loc[i+1,'GenBank Code'][0]} : "
              f"{corr}, p-value: {p_value}")
    except IndexError:
        pass




######GRAPHICS#############


sns.barplot(data=df_long, x='Position', y='Value', hue='Metric')
plt.title('Comparison of % GC Content and GC Content Random')
plt.xticks(rotation = 20, size=7)
plt.xlabel('Position')
plt.ylabel('GC Content')
plt.legend(loc='best')
plt.savefig('barplot_comparison_gc.jpg', dpi=600)
plt.close()

sns.lineplot(data=df_long, x='Position', y='Value', hue='Metric')
plt.title('Comparison of % GC Content and GC Content Random')
plt.xticks(rotation = 20, size=7)
plt.xlabel('Position')
plt.xlim(0,1400)
plt.ylabel('GC Content')
plt.legend(loc='best')
plt.savefig('lineplot_comparison_gc.jpg', dpi=600)
plt.close()


# Histogram of % GC Content
sns.histplot(data=df_gc_content_lists_exploded, x='% GC Content', hue='GenBank Code', kde=True, bins=30, fill=True)
plt.title('Histogram of % GC Content')
plt.xlabel('% GC Content')
plt.ylabel('Frequency')
plt.savefig('histogram_gc_content.jpg', dpi=600)
plt.close()


# Graphics of density
plt.figure(1, figsize=(11,6))

plt.subplot(1,2,1)
sns.kdeplot(data=df_gc_content_lists_exploded, x='% GC Content', hue='GenBank Code', fill=True)
plt.title('Density of % GC Content by GenBank Code')
plt.xlabel('% GC Content')
plt.ylabel('Density')
plt.grid(True)

plt.subplot(1,2,2)
sns.kdeplot(data=df_gc_content_lists_exploded, x='% GC Content Random', hue='GenBank Code', fill=True)
plt.title('Density of % GC Content Random by GenBank Code')
plt.xlabel('% GC Content Random')
plt.ylabel('Density')
plt.grid(True)
plt.tight_layout()
plt.savefig('density_gc.jpg', dpi=600)
plt.close()


# Lineplot of % GC Content by position in the sequence
plt.figure(2,figsize=(11,6))

plt.subplot(1,2,1)
sns.lineplot(data=df_gc_content_lists_exploded, x='Position', y='% GC Content', hue='GenBank Code', marker='o')
plt.legend(title="GenBank Code", fontsize="7", loc='best')
plt.title("Tendency of % GC Content by position in the sequence")
plt.xlim(0,1400)
plt.xlabel('Position')
plt.ylabel('% GC Content')
plt.grid(True)

plt.subplot(1,2,2)
# Lineplot of % GC Content by position in the sequence
sns.lineplot(data=df_gc_content_lists_exploded, x='Position', y='% GC Content', hue='GenBank Code', marker='h')
plt.legend(title='GenBank Code', fontsize='7', loc='best')
plt.title('Tendency of % GC Content with % GC for each whole sequence')
plt.xlabel('Position')
plt.ylabel('% GC Content')
plt.grid(True)
plt.xlim(0,1400)      #Adjust in accordance to the avg length of your sequences.
for i in range(0, number_of_sequences): #Since there's a clear tendency on the % GC, I did all the labels in the same color.
    plt.axhline(df_gc_content.loc[i,'% GC Content'],linestyle='-', color='darkcyan', label=df_gc_content.loc[i, 'GenBank Code'])
plt.tight_layout()
plt.savefig('lineplot_gc.jpg', dpi=600)
plt.close()


#Violin Plot
plt.figure(3, figsize=(14,6))

plt.subplot(1,2,1)
sns.violinplot(data=df_gc_content_lists_exploded, x='GenBank Code', y='% GC Content', hue='GenBank Code', split=True)
plt.title('Distribution of % GC Content by GenBank Code')
plt.xticks(rotation = 20, size=7)
plt.xlabel('GenBank Code')
plt.ylabel('% GC Content')

plt.subplot(1,2,2)
sns.violinplot(data=df_gc_content_lists_exploded, x='GenBank Code', y='Quartile', hue='GenBank Code', split=True)
plt.title('Distribution of sequences by Quartiles')
plt.xticks(rotation = 20, size=7)
plt.xlabel('GenBank Code')
plt.ylabel('Quartiles')
plt.tight_layout()
plt.savefig('violinplot_gc_quartiles.jpg', dpi=600)
plt.close()


# Boxplot for % GC Content by quartiles
plt.figure(4, facecolor='white', figsize=(11,6))

plt.subplot(1,2,1)
sns.boxplot(x='Quartile', y='% GC Content', data=df_gc_content_lists_exploded, hue='GenBank Code')
plt.title('Distribution of % GC Content by Quartiles')
plt.xlabel('Quartiles')
plt.ylabel('% GC Content')
plt.legend(fontsize='9', loc ='best')
plt.grid()

plt.subplot(1,2,2)

sns.boxplot(x='Quartile Random', y='% GC Content Random', data=df_gc_content_lists_exploded, hue='GenBank Code')
plt.title('Distribution of % GC Content Random by Quartiles')
plt.xlabel('Quartiles')
plt.ylabel('% GC Content Random')
plt.legend(fontsize='9', loc ='best')
plt.grid()
plt.tight_layout()
plt.savefig('boxplot_gc_random.jpg', dpi=600)
plt.close()


# Scatterplot of % GC Content by position and GenBank Code with regression lines
sns.lmplot(x='Position', y='% GC Content', data=df_gc_content_lists_exploded, hue='GenBank Code',
           scatter_kws={'alpha':0.5})
#plt.title('Dispersion and regression lines for % GC Content by position')
plt.xlabel('Position')
plt.ylabel('% GC Content')
plt.grid()
plt.savefig('scatterplot_gc.jpg', dpi=600)
plt.close()


# Heatmap for % GC Content shows correlation
sns.heatmap(heatmap_data, cmap='coolwarm', annot=False)
plt.title('Heatmap of % GC Content between sequences by position')
plt.xlabel('Position')
plt.xticks(rotation = 20, size=7)
plt.ylabel('GenBank Code')
plt.yticks(rotation = 0, size=8)
plt.savefig('heatmap_gc.jpg', dpi=600)
plt.close()


####################################### Creates the LaTeX Document#########################################
doc = Document()

# Document title
doc.preamble.append(NoEscape(r'\usepackage[headheight=20pt, top=2cm]{geometry}'))
doc.preamble.append(Command('title', f'Analysis report for the "{gene_name}" gene in {number_of_sequences} supposedly alike sequences'))
doc.preamble.append(Command('author', 'Luis Jimenez'))
doc.preamble.append(Command('date', NoEscape(r'\today')))
doc.append(NoEscape(r'\maketitle'))

# Section
with doc.create(Section('Introduction')):
    doc.append('The following analysis is made by calculating the guanine-cytosine percentage '
               'of a specific gene in several viral sequences and comparing them to hopefully '
               'observe any similarities. The Spearman test for correlations was also made between '
               'some sequences. It is know that viruses are classified in accordance '
               'to their characteristics , such as capside type, DNA chain type, infective '
               f'capacity, etc. In this case, the {gene_name} gene in retrovirus is key in the '
               "replication cycle, starting at the formation of new virions until its release and maturation. "
               'It codifies for the structural proteins needed for assembly and keeping the integrity '
               'of the viral particle. The analyzed sequences extracted from GeneBank are the following: '
               f'{genbank_code_list}.')


# Graphics
with doc.create(Section('Results')):
    with doc.create(Subsection('Histogram of % GC Content')):
        with doc.create(Figure(position='h!')) as graphic:
            graphic.add_image('histogram_gc_content.jpg', width=NoEscape(r'0.8\textwidth'))
            graphic.add_caption('Distribution of % GC content in the sequences analyzed by fragments at a time.')

    with doc.create(Subsection('Density of % GC Content and % GC Content for the randomized sequences')):
        with doc.create(Figure(position='!htbp')) as graphic1:
            graphic1.add_image('density_gc.jpg', width=NoEscape(r'1\textwidth'))
            graphic1.add_caption('Density of % GC content in the sequences analyzed by fragments at a time.')

    with doc.create(Subsection('Lineplots of % GC Content by position in the sequence')):
        with doc.create(Figure(position='!htbp')) as graphic2:
            graphic2.add_image('lineplot_gc.jpg', width=NoEscape(r'1\textwidth'))
            graphic2.add_caption('Distribution of % GC content in the sequences analyzed by fragments at a time.')

    with doc.create(Subsection('Violinplots for the distribution of % GC Content by sequence and quartiles')):
        with doc.create(Figure(position='!htbp')) as graphic2:
            graphic2.add_image('violinplot_gc_quartiles.jpg', width=NoEscape(r'1\textwidth'))
            graphic2.add_caption('Distribution of % GC content by GenBank Code and quartiles.')

    with doc.create(Subsection('Scatterplot of % GC Content by position and GenBank Code with regression lines')):
        with doc.create(Figure(position='!htbp')) as graphic3:
            graphic3.add_image('scatterplot_gc.jpg', width=NoEscape(r'0.8\textwidth'))
            graphic3.add_caption('Dispersion and regression lines for % GC Content by position')

    with doc.create(Subsection('Boxplots for % GC Content by quartiles')):
        with doc.create(Figure(position='h!')) as graphic4:
            graphic4.add_image('boxplot_gc_random.jpg', width=NoEscape(r'1\textwidth'))
            graphic4.add_caption('Distribution of % GC Content and Random GC Content by Quartiles')

    with doc.create(Subsection('Heatmap for % GC Content shows similarities')):
        with doc.create(Figure(position='h!')) as graphic6:
            graphic6.add_image('heatmap_gc.jpg', width=NoEscape(r'0.8\textwidth'))
            graphic6.add_caption('Comparison of % GC Content and GC Content Random')

    with doc.create(Subsection('Bar and line plots for % GC Content and Random GC Content for whole sequences')):
        with doc.create(Figure(position='h!')) as graphic5:
            graphic5.add_image('barplot_comparison_gc.jpg', width=NoEscape(r'0.8\textwidth'))
            graphic5.add_caption('Comparison of % GC Content and GC Content Random - Barplot')

        with doc.create(Figure(position='h!')) as graphic7:
            graphic7.add_image('lineplot_comparison_gc.jpg', width=NoEscape(r'0.8\textwidth'))
            graphic7.add_caption('Comparison of % GC Content and GC Content Random - Lineplot')




# Compile the document
doc.generate_pdf('analysis_report' ,compiler='pdflatex', clean_tex=False)
print('PDF generated as analysis_report.pdf')

