from Bio import SeqIO
from fuzzywuzzy import fuzz
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
from geneFinder import Orfs
import matplotlib.pyplot as plt
import numpy as np
import argparse

'''
SequenceComparison.py
Author: Roni Altshuler
This class generates various graphs/figures comparing the sequence similarity scores and open reading frame (ORF)
similarity scores using FuzzyWuzzy string comparison and uses the Orfs class from Derfel's geneFinder.py file
'''

class SequenceSimilarity:
    def __init__(self, filename):
        '''
        Constructor for SequenceSimilarity class that takes a file name as input and performs some operations on the
        contents of the fast file
        Reads the contents of the file with the filename input and uses SeqIO.parse from the Bio.SeqIO module to
        parse the file as a FASTA file. It then stores the resulting list of records in the self.records instance
        variable.
        '''
        self.records = list(SeqIO.parse(filename, "fasta"))
        self.headerList = [r.id for r in self.records] # Creates a list of all the sequences in the file
        self.sequenceList = [str(r.seq) for r in self.records] # Creates a list of all the headers in the file
        self.combineList = list(zip(self.headerList, self.sequenceList)) # Combining the header and sequence lists into a single list of tuples
        self.result = []
        for i in range(len(self.sequenceList)): # Iterating over each pair of sequences
            result = []
            for j in range(len(self.sequenceList)):
                result.append(fuzz.ratio(self.sequenceList[i], self.sequenceList[j])) # Adding to result list the similarity scores of the sequences
            self.result.append(result)

    def generate_scatter_plot(self):
        '''
        Generates a scatter plot showing the FuzzyWuzzy sequence similarity scores between sequences in a fasta file.
        Each dot in the plot represents a pair of sequences, with the x-axis representing the similarity score for
        the first sequence and the y-axis representing the similarity score for the second sequence.
        '''
        font = {'family': 'serif', 'color': 'black', 'size': 17}
        plt.style.use('seaborn-v0_8')
        fig, ax = plt.subplots(figsize=(15, 8))  # Increase the figure size

        # Define a dictionary to hold the shape for each sequence
        shape_dict = {}
        shapes = ['o', '^', 's', 'D', '*', 'X', 'P', 'p', 'H', 'h', 'v', '<', '>', 'd']
        for i, header in enumerate(self.headerList):
            shape_dict[header] = shapes[i % len(shapes)]  # Cycle through the shapes list for each new sequence

        # Define a color map based on the number of sequences
        cmap = plt.get_cmap('tab10')
        colors = [cmap(i) for i in np.linspace(0, 1, len(self.headerList))]

        # Create a dictionary to map header names to colors
        header_colors = dict(zip(self.headerList, colors))

        for i in range(len(self.sequenceList)):
            plt.scatter(np.arange(1, len(self.sequenceList) + 1), self.result[i], label=self.headerList[i],
                        linewidths=1.5,
                        c=[header_colors[self.headerList[i]]] * len(self.sequenceList), s=100, edgecolors='black',
                        marker=shape_dict[self.headerList[i]])  # Use the shape_dict to get the shape for each sequence

        plt.title('FuzzyWuzzy Sequence Similarity Scores Scatter Plot', font, fontweight='bold')
        plt.xlabel('Sequence Index', font, fontweight='bold')
        plt.ylabel('FuzzyWuzzy Score', font, fontweight='bold')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=14)  # Increase legend font size
        plt.tight_layout()  # Automatically adjust subplot parameters to give specified padding
        plt.savefig('FuzzyWuzzySequenceSimilarityScoresScatterPlot.png', dpi=600)

    def generate_violin_plots(self):
        '''
        Generates a violin plot to display the ORF similarity scores between sequences in a fasta file.
        The violin plot displays the ORF similarity scores for each sequence.
        '''
        font = {'family': 'serif', 'color': 'black', 'size': 17}
        plt.style.use('seaborn-v0_8')
        fig, ax1 = plt.subplots(nrows=1, ncols=1,
                                figsize=(8, 8))  # Increase the figure size and create a single subplot

        # Define a color map based on the number of sequences
        cmap = plt.get_cmap('tab10')
        colors = [cmap(i) for i in np.linspace(0, 1, len(self.headerList))]

        # Violin plot
        orf_scores = []
        for i in range(len(self.sequenceList)):
            score_list = []
            for j in range(len(self.sequenceList)):
                score_list.append(self.calculate_similarity_score(self.sequenceList[i], self.sequenceList[j]))
            orf_scores.append(score_list)
        ax1.violinplot(orf_scores, showmeans=True, showextrema=True, showmedians=True)
        ax1.set_title('ORF Similarity Scores Violin Plot', font, fontweight='bold')
        ax1.set_xlabel('Sequences', font, fontweight='bold')
        ax1.set_ylabel('ORF Score', font, fontweight='bold')
        ax1.set_xticks(np.arange(1, len(self.headerList) + 1))
        ax1.tick_params(axis='y', which='major', labelsize=14)  # Increase tick font size

        plt.tight_layout()  # Automatically adjust subplot parameters to give specified padding
        plt.savefig('ORFSimilarityScoresViolinPlot.png', dpi=600)

    def generate_histogram_plot(self):
        '''
        Generates a histogram plot of the FuzzyWuzzy sequence similarity scores for all pairs of sequences in a FASTA file.
        Histogram shows the distribution of FuzzyWuzzy sequence similarity scores for all pairs of sequences in your FASTA file.
        The x-axis represents the range of possible FuzzyWuzzy scores, which range from 0 (completely dissimilar sequences) to 100 (identical sequences).
        The y-axis represents the frequency of sequence pairs with a given FuzzyWuzzy score. Each bar in the histogram represents a range of FuzzyWuzzy scores,
        and the height of the bar indicates the frequency of sequence pairs with a score in that range.
        '''
        font = {'family': 'serif', 'color': 'black', 'size': 17}
        plt.style.use('seaborn-v0_8')
        fig, ax = plt.subplots(figsize=(15, 8))  # Increase the figure size

        # Create a dictionary to hold the sequence similarity scores for each header
        scores_dict = {}
        for header in self.headerList:
            scores_dict[header] = []

        # Compute sequence similarity scores for all pairs of sequences
        for i in range(len(self.sequenceList)):
            for j in range(i + 1, len(self.sequenceList)):
                score = fuzz.ratio(self.sequenceList[i], self.sequenceList[j])
                scores_dict[self.headerList[i]].append(score)
                scores_dict[self.headerList[j]].append(score)

        # Plot the histogram for each header's sequence similarity scores
        for header in self.headerList:
            scores = scores_dict[header]
            plt.hist(scores, alpha=0.5, label=header, bins=range(0, 101, 5))

        plt.title('FuzzyWuzzy Sequence Similarity Scores Histogram', font, fontweight='bold')
        plt.xlabel('FuzzyWuzzy Score', font, fontweight='bold')
        plt.ylabel('Frequency', font, fontweight='bold')
        # plt.legend(loc='center left', fontsize=12)  # Add a legend for headers
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=14)  # Increase legend font size
        plt.tight_layout()  # Automatically adjust subplot parameters to give specified padding
        plt.savefig('FuzzyWuzzySequenceSimilarityScoresHistogram.png', dpi=600)

    def find_orfs(self, seq):
        '''
        Finds all Open Reading Frames (ORFs) in the given sequence.
        Returns a list of tuples, where each tuple contains the start
        and end positions of an ORF.
        '''
        orfs = []
        for orf in Orfs.findOrfs(seq, minLength=100):
            orfs.append(orf.getOrfIndex)
        return orfs

    def calculate_similarity_score(self, seq1, seq2):
        '''
        Calculates the FuzzyWuzzy similarity score between two sequences.
        '''
        return fuzz.ratio(seq1, seq2)
        # return fuzz.partial_ratio(seq1, seq2)

    def generate_orf_similarity_plot(self):
        '''
        Generates a heatmap plot that displays the similarity scores between the open reading frames (ORFs) of two DNA sequences.
        The x-axis and y-axis represent the headers of the DNA sequences. The color of each square in the heatmap represents
        the similarity score between the ORFs of the corresponding pair of sequences. The higher the score, the redder the
        square, and the lower the score, the bluer the square.
        '''
        fig, ax = plt.subplots(figsize=(10, 10))  # Reduce the figure size

        # Calculate ORF similarity scores
        scores = []
        for i, (header_i, seq_i) in enumerate(self.combineList):
            row = []
            for j, (header_j, seq_j) in enumerate(self.combineList):
                orfs_i = self.find_orfs(seq_i)
                orfs_j = self.find_orfs(seq_j)
                if orfs_i and orfs_j:
                    orf_similarities = []
                    for orf_i in orfs_i:
                        for orf_j in orfs_j:
                            # Compare only the longest ORFs
                            orf_i_len = orf_i[1] - orf_i[0]
                            orf_j_len = orf_j[1] - orf_j[0]
                            if orf_i_len >= orf_j_len:
                                score = self.calculate_similarity_score(seq_i[orf_i[0]:orf_i[1]], seq_j[orf_j[0]:orf_j[1]])
                            else:
                                score = self.calculate_similarity_score(seq_i[orf_i[0]:orf_i[1]], seq_j[orf_j[0]:orf_j[1]])
                            orf_similarities.append(score)
                    if orf_similarities:
                        row.append(sum(orf_similarities) / len(orf_similarities))
                    else:
                        row.append(0)
                else:
                    row.append(0)
            scores.append(row)

        # Plot ORF similarities
        plt.imshow(scores, cmap='coolwarm', interpolation='nearest')
        plt.colorbar()
        plt.xticks(range(len(self.headerList)), self.headerList, rotation=90, fontsize=12)  # Reduce the font size
        plt.yticks(range(len(self.headerList)), self.headerList, fontsize=12)  # Reduce the font size
        plt.title('ORF Similarity Scores Heatmap', fontdict={'fontsize': 20, 'fontweight': 'bold'})
        plt.xlabel('Sequences', fontdict={'fontsize': 16, 'fontweight': 'bold'})
        plt.ylabel('Sequences', fontdict={'fontsize': 16, 'fontweight': 'bold'})
        plt.tight_layout()  # Automatically adjust subplot parameters to give specified padding
        plt.savefig('ORFSimilarityScores.png', dpi=600)

    def generate_orf_similarity_plot_with_dendogram(self):
        '''
        Generates a heatmap plot that displays the similarity scores between the open reading frames (ORFs) of two DNA sequences with a dendogram.
        The x-axis and y-axis represent the headers of the DNA sequences. The color of each square in the heatmap represents
        the similarity score between the ORFs of the corresponding pair of sequences. A dendrogram is included to show the hierarchical
        clustering of the sequences based on their ORF similarity. In a dendrogram, the colors represent different clusters
        or groups of data points. The dendrogram is constructed by recursively merging the data points or clusters that are
        the most similar to each other based on a chosen distance metric. The y-axis of the dendrogram represents the
        distance or dissimilarity between the clusters being merged. The height at which two branches merge corresponds
        to the distance between the clusters at the time of the merge. Therefore, the numbers on the y-axis represent
        the distance between the clusters. The clusters that are more similar to each other are merged at lower heights,
        while the clusters that are less similar are merged at higher heights. The distance values can be normalized or
        scaled to a specific range for easier interpretation.
        '''
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10), gridspec_kw={'width_ratios': [1, 0.1]})

        # Calculate ORF similarity scores
        scores = []
        for i, (header_i, seq_i) in enumerate(self.combineList):
            row = []
            for j, (header_j, seq_j) in enumerate(self.combineList):
                orfs_i = self.find_orfs(seq_i)
                orfs_j = self.find_orfs(seq_j)
                if orfs_i and orfs_j:
                    orf_similarities = []
                    for orf_i in orfs_i:
                        for orf_j in orfs_j:
                            # Compare only the longest ORFs
                            orf_i_len = orf_i[1] - orf_i[0]
                            orf_j_len = orf_j[1] - orf_j[0]
                            if orf_i_len >= orf_j_len:
                                score = self.calculate_similarity_score(seq_i[orf_i[0]:orf_i[1]], seq_j[orf_j[0]:orf_j[1]])
                            else:
                                score = self.calculate_similarity_score(seq_i[orf_i[0]:orf_i[1]], seq_j[orf_j[0]:orf_j[1]])
                            orf_similarities.append(score)
                    if orf_similarities:
                        row.append(sum(orf_similarities) / len(orf_similarities))
                    else:
                        row.append(0)
                else:
                    row.append(0)
            scores.append(row)

        # Convert scores matrix into a condensed distance matrix
        distances = pdist(scores)

        # Generate linkage matrix from condensed distance matrix
        linkage_matrix = hierarchy.linkage(distances, method='complete')
        # Calculate dendrogram
        # linkage_matrix = hierarchy.linkage(np.asarray(scores), method='complete')
        dendro = hierarchy.dendrogram(linkage_matrix, ax=ax2, orientation='right')
        ax2.set_title('Dendrogram', fontdict={'fontsize': 20, 'fontweight': 'bold'})
        ax2.set_xticks([])

        # Plot ORF similarities
        im = ax1.imshow(scores, cmap='coolwarm', interpolation='nearest')
        ax1.set_xticks(range(len(self.headerList)))
        ax1.set_xticklabels(self.headerList, rotation=90, fontsize=12)
        ax1.set_yticks(range(len(self.headerList)))
        ax1.set_yticklabels(self.headerList, fontsize=12)
        ax1.set_title('ORF Similarity Scores Heatmap with Dendrogram', fontdict={'fontsize': 20, 'fontweight': 'bold'})
        ax1.set_xlabel('Sequences', fontdict={'fontsize': 16, 'fontweight': 'bold'})
        ax1.set_ylabel('Sequences', fontdict={'fontsize': 16, 'fontweight': 'bold'})
        fig.colorbar(im, ax=ax1)
        plt.tight_layout()
        plt.savefig('ORFSimilarityScores_with_dendrogram.png', dpi=600)

    def generate_orf_similarity_histogram(self):
        '''
        Generates a histogram of the similarity scores between all ORFs in the sequences.
        Saves the histogram as a PNG file.
        '''
        scores = []
        for i, (header_i, seq_i) in enumerate(self.combineList):
            for j, (header_j, seq_j) in enumerate(self.combineList):
                orfs_i = self.find_orfs(seq_i)
                orfs_j = self.find_orfs(seq_j)
                if orfs_i and orfs_j:
                    for orf_i in orfs_i:
                        for orf_j in orfs_j:
                            # Compare only the longest ORFs
                            orf_i_len = orf_i[1] - orf_i[0]
                            orf_j_len = orf_j[1] - orf_j[0]
                            if orf_i_len >= orf_j_len:
                                score = self.calculate_similarity_score(seq_i[orf_i[0]:orf_i[1]],
                                                                        seq_j[orf_j[0]:orf_j[1]])
                            else:
                                score = self.calculate_similarity_score(seq_i[orf_i[0]:orf_i[1]],
                                                                        seq_j[orf_j[0]:orf_j[1]])
                            scores.append(score)

        # Plot histogram of ORF similarity scores
        fig, ax = plt.subplots(figsize=(10, 10))
        n, bins, patches = plt.hist(scores, bins=20, color='blue')
        plt.title('Histogram of ORF Similarity Scores', fontdict={'fontsize': 20, 'fontweight': 'bold'})
        plt.xlabel('Similarity Score', fontdict={'fontsize': 16, 'fontweight': 'bold'})
        plt.ylabel('Frequency', fontdict={'fontsize': 16, 'fontweight': 'bold'})

        # Add a vertical line indicating the mean similarity score
        mean_score = np.mean(scores)
        ax.axvline(mean_score, color='red', linestyle='dashed', linewidth=2)
        ax.text(mean_score + 0.1, ax.get_ylim()[1] * 0.9, f"Mean = {mean_score:.2f}", fontsize=12, color='red')

        plt.savefig('ORFSimilarityHistogram.png', dpi=600)

class CommandLine:
    '''
    Command line class that takes a filename argument from the user and passes it to a SequenceSimilarity class to generate different kinds of plots.
    '''
    def __init__(self):
        '''
        Sets-up the command line arguments to run the script
        Usage: python3 <filename.py> <fastafile.fa>
        Example: python3 SequenceComparison.py A23EEV21.results.fa
        '''
        self.parser = argparse.ArgumentParser(description='Calculate sequence similarity using FuzzyWuzzy')
        self.parser.add_argument('filename', type=str, help='FASTA file containing sequences to compare')
        self.args = self.parser.parse_args()

    def run(self):
        '''
        Calls methods from the SequenceSimilarity class and runs methods to generate plots
        '''
        ss = SequenceSimilarity(self.args.filename)
        ss.generate_scatter_plot()
        ss.generate_violin_plots()
        ss.generate_histogram_plot()
        ss.generate_orf_similarity_plot()
        ss.generate_orf_similarity_plot_with_dendogram()
        ss.generate_orf_similarity_histogram()

if __name__ == '__main__':
    cmd = CommandLine()
    cmd.run()