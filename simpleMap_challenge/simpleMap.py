import array
import sys
import numpy
import pysam
import argparse
import logging 
import time
logger = logging.getLogger()

"""See the comments below to see the code you need to complete.
"""

class MinimizerIndexer(object):
    """ Simple minimizer based substring-indexer. 
    
    Please read: https://doi.org/10.1093/bioinformatics/bth408
    
    Related to idea of min-hash index and other "sketch" methods.
    """
    def __init__(self, targetString, w, k, t):
        """ The target string is a string/array of form "[ACGT]*".
        
        Stores the lexicographically smallest k-mer in each window of length w, such that w >= k positions. This
        smallest k-mer is termed a minmer. 
        
        If a minmer occurs in the target sequence more than t times as a minmer then it is omitted from the index, i.e. if the given minmer (kmer) is a minmer
        in more than t different locations in the target string. Note, a minmer may be the minmer for more than t distinct windows
        and not be pruned, we remove minmers only if they have more than t distinct occurrences as minmers in the sequence.
        """
        
        self.targetString = targetString
        self.w = w
        self.k = k
        self.t = t # If a minmer occurs more than t times then its entry is removed from the index
        # This is a heuristic to remove repetitive minmers that would create many spurious alignments between
        # repeats
        
        # Hash of minmers to query locations, stored as a map whose keys
        # are minmers and whose values are lists of the start indexes of
        # occurrences of the corresponding minmer in the targetString, 
        # sorted in ascending order of index in the targetString.
        #
        # For example if k = 2 and w = 4 and targetString = "GATTACATTT"
        #
        # GATTACATTT
        # GATT (AT)
        #  ATTA (AT)
        #   TTAC (AC)
        #    TACA (AC)
        #     ACAT (AC)
        #      CATT (AT)
        #       ATTT (AT)
        #
        # then self.minimizerMap = { "AT":(1,6), "AC":(4,) }
        self.minimizerMap = {}

        # Code to complete to build index - you are free to define additional functions
        corrected_range = len(self.targetString) - self.w + 1
        repeats = 0

        temp_dict = {}
        for i in range(0, corrected_range):
            curr_w = self.targetString[i:i + self.w]
            k_size = len(curr_w) - self.k
            k_collection = set()

            for i_k in range(0, k_size + 1):
                curr_k = curr_w[i_k:i_k + self.k]
                k_collection.add((curr_k, i + i_k))
            
            min_kmer = min(k_collection, key=lambda x:(x[0], x[1]))

            if min_kmer[0] in temp_dict and len(temp_dict[min_kmer[0]]) >= self.t:
                continue

            if min_kmer[0] not in temp_dict:
                temp_dict[min_kmer[0]] = set([min_kmer[1]])
            else:
                temp_dict[min_kmer[0]].add(min_kmer[1])
        
        self.minimizerMap = temp_dict
        # print(self.targetString, self.minimizerMap)
        
        # unit-test

    def getMatches(self, searchString):
        """ Iterates through search string finding minmers in searchString and
        yields their list of minmer occurrences in targetString, each as a pair of (x, (y,)*N),
        where x is the index in searchString and y is an occurrence in targetString.
        
        For example if k = 2 and w = 4 and targetString = "GATTACATTT" and searchString = "GATTTAC"
        then self.minimizerMap = { "AT":(1,6), "AC":(4,) }
        and getMatches will yield the following sequence:
        (1, (1,6)), (5, (4,))
        
        You will need to use the "yield" keyword
        """
        # Code to complete - you are free to define additional functions
        for i in range(0, len(searchString) - self.k + 1):
            curr_str = searchString[i:i + self.k]
            if curr_str in self.minimizerMap:
                yield i, self.minimizerMap[curr_str]

t_str = "GATTACATTT"
new_min = MinimizerIndexer(t_str, 4, 2, 1000)
result = [x for x in new_min.getMatches("GATTTAC")]
result.sort(key=lambda x:(x[0],x[1]))
#print(result)

class SeedCluster:
    """ Represents a set of seeds between two strings.
    """
    def __init__(self, seeds):
        """ Seeds is a list of pairs [ (x_1, y_1), (x_2, y_2), ..., ], each is an instance of a seed 
        (see static cluster seeds method below: static methods: https://realpython.com/blog/python/instance-class-and-static-methods-demystified/)
        """
        seeds = list(seeds)
        seeds.sort()
        self.seeds = seeds
        #print(f'seeds: {self.seeds}')
        # Gather the minimum and maximum x and y coordinates
        self.minX = seeds[0][0]
        self.maxX = seeds[-1][0]
        ys = [y for x,y in seeds] # python3: map(lambda...) changed to list comprehension
        self.minY = min(ys)
        self.maxY = max(ys)

    @staticmethod
    def clusterSeeds(seeds, l):
        """ Cluster seeds (k-mer instances) in two strings. This is a static constructor method that creates a set
        of SeedCluster instances.
        
        Here seeds is a list of tuples, each tuple has the form (x, (y_1, y_2, ... )), where x is the coordinate
        in the first string and y_1, y_2, ... are coordinates in the second string. Each pair of x and y_i
        is an occurence of a shared k-mer in both strings, termed a *seed*, such that the k-mer 
        occurrence starts at position x in the first string and starts at position y_i in the second string.
        The input seeds list contains no duplicates and is sorted in ascending order, 
        first by x coordinate (so each successive tuple will have a greater  
        x coordinate), and then each in tuple the y coordinates are sorted in ascending order.
        
        Two seeds (x_1, y_1), (x_2, y_2) are *close* if the absolute distances | x_2 - x_1 | and | y_2 - y_1 |
        are both less than or equal to l.   
        
        Consider a *seed graph* in which the nodes are the seeds, and there is an edge between two seeds if they
        are close. clusterSeeds returns the connected components of this graph
        (https://en.wikipedia.org/wiki/Connected_component_(graph_theory)).
        
        The return value is a Python set of SeedCluster object, each representing a connected component of seeds in the 
        seed graph.
        
        (QUESTION 1): The clustering of seeds is very simplistic. Can you suggest alternative strategies by
        which the seeds could be clustered, and what the potential benefits such alternative strategies could
        have? Consider the types of information you could use.  
        """ 
        
        # Code to complete - you are free to define other functions as you like
        
        #After much struggle, Roni inspired me to use BFS

        pairs = []
        if len(seeds) > 0:
            for seed in seeds:
                if len(seed[1]) > 0:
                    for ys in seed[1]:
                        pairs.append((seed[0], ys))
    
        final_cluster = set()

        def BFS(pair, pairs):
            queue = [pair]
            visited = set()
            visited.add(pair)
            cluster = set()
            cluster.add(pair)

            curr_pairs = set(pairs)

            while len(queue) > 0:
                x1,y1 = queue.pop(0)
                
                found_pairs = set()
                for x2,y2 in curr_pairs:
                    if abs(x2-x1) <= l and abs(y2-y1) <=l and (x2,y2) not in cluster:
                        if (x2,y2) not in visited:
                            queue.append((x2,y2))
                            visited.add((x2,y2))
                        cluster.add((x2,y2))
                        found_pairs.add((x2,y2))

                curr_pairs = curr_pairs.difference(found_pairs)
                #print(curr_pairs)
            return cluster

        for pair in pairs:
            cluster = BFS(pair,pairs)
            final_cluster.add(SeedCluster(cluster))
        
        return final_cluster

#result = [(1, (6,)), (2, ()), (3, (1, 2)), (4, (2,)), (5, ())] # correct = exp: [(1, 6), (3, 1), (3, 2), (4, 2)], testing:(3, 1), result: [(3, 1), (3, 2), (4, 2)] l = 4
#result = [(0, ()), (1, ()), (2, (2,)), (3, (2,)), (5, ()), (6, (1,))] # correct = exp: [(2, 2)], testing:(2, 2), result: [(2, 2), (3, 2)]
#result = [(1, (8,)), (2, (2, 4, 5))] # correct = exp: [(2, 4), (2, 5)], testing:(2, 4), result: [(2, 4)] l = 1
#sc = SeedCluster.clusterSeeds(result, 1) 


class SmithWaterman(object):
    def __init__(self, string1, string2, gapScore=-2, matchScore=3, mismatchScore=-3):
        """ Finds an optimal local alignment of two strings.
        
        Implements the Smith-Waterman algorithm: 
        https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
        
        (QUESTION 2): The Smith-Waterman algorithm finds the globally optimal local alignment between to 
        strings, but requires O(|string1| * |string2|) time. Suggest alternative strategies you could implement
        to accelerate the finding of reasonable local alignments. What drawbacks might such alternatives have?
        """
        # Code to complete to compute the edit matrix

        self.string1 = string1
        self.string2 = string2
        self.gapScore = gapScore
        self.matchScore = matchScore
        self.mismatchScore = mismatchScore

        self.H = numpy.zeros(shape=[len(string1) + 1, len(string2) + 1], dtype=int)

        self.highestScore = (-1 * float("inf"), (0,0)) #highest score and its coords
        self.scoring = lambda i,j: self.matchScore if self.string1[i-1] == self.string2[j-1] else self.mismatchScore
        self.scoreTuple = lambda i,j: (self.H[i-1,j-1] + self.scoring(i,j), \
            self.H[i-1,j] + self.gapScore, \
            self.H[i,j-1] + self.gapScore \
            )
        
        for i in range(1, len(self.H) ):
            for j in range(1, len(self.H[0])):
                self.H[i,j] = max(self.scoreTuple(i,j))
                if self.H[i,j] > self.highestScore[0]:
                    self.highestScore = (self.H[i,j], (i,j))
        #print(self.H, self.highestScore)

                
    def getAlignment(self):
        """ Returns an optimal local alignment of two strings. Alignment
        is returned as an ordered list of aligned pairs.
        
        e.g. For the two strings GATTACA and CTACC an optimal local alignment
        is (GAT)TAC(A)
             (C)TAC(C)
        where the characters in brackets are unaligned. This alignment would be returned as
        [ (3, 1), (4, 2), (5, 3) ] 
        
        In cases where there is a tie between optimal sub-alignments use the following rule:
        Let (i, j) be a point in the edit matrix, if there is a tie between possible sub-alignments
        (e.g. you could chooose equally between different possibilities), choose the (i, j) to (i-1, j-1)
        (match) in preference, then the (i, j) to (i-1, j) (insert in string1) in preference and
        then (i, j) to (i, j-1) (insert in string2).
        """
        # Code to complete - generated by traceback through matrix to generate aligned pairs

        #inspired by NB 4
        score = self.highestScore
        alignment = []
        
        i, j = score[1][0], score[1][1]

        tracebackScore = lambda i,j: ( self.H[i-1, j-1] + self.scoring(i, j), \
            self.H[i-1, j] + self.gapScore, \
            self.H[i, j-1] + self.gapScore
        ) # 0: match/mismatch, 1: vertical, 2: horizontal

        while self.H[i,j] != 0:
            curr_score = tracebackScore(i,j)
            probable_predecessor = max(curr_score)

            if probable_predecessor == curr_score[0]:
                i -= 1
                j -= 1
                alignment.append((i,j))
            elif probable_predecessor == curr_score[1]:
                i -= 1
            elif probable_predecessor == curr_score[2]:
                j -= 1
        
        alignment.reverse()
        return alignment


    def getMaxAlignmentScore(self):
        """ Returns the maximum alignment score
        """
        # Code to complete
        return self.highestScore[0]
    
# str1 = 'GATTACA'
# str2 = 'CTACC'
# sw = SmithWaterman(str1, str2)

#print(sw.getAlignment())

def simpleMap(targetString, minimizerIndex, queryString, config):
    """ Function takes a target string with precomputed minimizer index and a query string
    and returns the best alignment it finds between target and query, using the given options specified in config.
    
    Maps the string in both its forward and reverse complement orientations.
    
    (QUESTION 3): The code below is functional, but very slow. Suggest ways you could potentially accelerate it, 
    and note any drawbacks this might have.
    """
    bestAlignment = [None]
    
    def mapForwards(queryString):
        """ Maps the query string forwards
        """
        # Find seed matches, aka "aligned kmers"
        seeds = list(minimizerIndex.getMatches(queryString))
        
        # For each cluster of seeds
        for seedCluster in SeedCluster.clusterSeeds(list(seeds), l=config.l):
            
            # Get substring of query and target to align
            queryStringStart = max(0, seedCluster.minX - config.c) # Inclusive coordinate
            queryStringEnd = min(len(queryString), seedCluster.maxX + config.k + config.c) # Exclusive coordinate
            querySubstring = queryString[queryStringStart:queryStringEnd]
            
            targetStringStart = max(0, seedCluster.minY - config.c) # Inclusive coordinate
            targetStringEnd = min(len(targetString), seedCluster.maxY + config.k + config.c) # Exclusive coordinate
            targetSubstring = targetString[targetStringStart:targetStringEnd]
            
            #print "target_aligning", targetStringStart, targetStringEnd, targetSubstring
            #print "query_aligning", queryStringStart, queryStringEnd, querySubstring
            
            # Align the genome and read substring
            alignment = SmithWaterman(targetSubstring, querySubstring, 
                                      gapScore=config.gapScore, 
                                      matchScore=config.matchScore,
                                      mismatchScore=config.mismatchScore)
            
            # Update best alignment if needed
            if bestAlignment[0] == None or alignment.getMaxAlignmentScore() > bestAlignment[0].getMaxAlignmentScore():
                bestAlignment[0] = alignment
        
        return bestAlignment
    
    def reverseComplement(string):
        """Computes the reverse complement of a string
        """
        rMap = { "A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
        return "".join(rMap[i] for i in string[::-1])
                
    # Run mapping forwards and reverse
    mapForwards(queryString)
    mapForwards(reverseComplement(queryString))
    
    return bestAlignment[0]

class Config():
    """ Minimal configuration class for handing around parameters
    """
    def __init__(self):
        self.w = 30
        self.k = 20
        self.t = 10
        self.l = 30
        self.c = 100
        self.gapScore=-2
        self.matchScore=3
        self.mismatchScore=-3
        self.logLevel = "INFO"
        
def main():
    # Read parameters
    config = Config()
    
    #Parse the inputs args/options
    parser = argparse.ArgumentParser(usage="target_fasta query_fastq [options]") # , version="%prog 0.1")

    parser.add_argument("target_fasta", type=str,
                        help="The target genome fasta file.")
    parser.add_argument("query_fastq", type=str,
                        help="The query sequences.")
    
    parser.add_argument("--w", dest="w", help="Length of minimizer window. Default=%s" % config.w, default=config.w)
    parser.add_argument("--k", dest="k", help="Length of k-mer. Default=%s" % config.k, default=config.k)
    parser.add_argument("--t", dest="t", help="Discard minmers that occur more frequently " 
                                            "in the target than t. Default=%s" % config.w, default=config.w)
    parser.add_argument("--l", dest="l", help="Cluster two minmers into the same cluster if within l bases of"
                                            " each other in both target and query. Default=%s" % config.l, default=config.l)
    parser.add_argument("--c", dest="c", help="Add this many bases to the prefix and suffix of a seed cluster in the"
                                            " target and query sequence. Default=%s" % config.c, default=config.c)
    parser.add_argument("--gapScore", dest="gapScore", help="Smith-Waterman gap-score. Default=%s" % 
                      config.gapScore, default=config.gapScore)
    parser.add_argument("--matchScore", dest="matchScore", help="Smith-Waterman match-score. Default=%s" % 
                      config.gapScore, default=config.gapScore)
    parser.add_argument("--mismatchScore", dest="mismatchScore", help="Smith-Waterman mismatch-score. Default=%s" % 
                      config.mismatchScore, default=config.mismatchScore)
    parser.add_argument("--log", dest="logLevel", help="Logging level. Default=%s" % 
                      config.logLevel, default=config.logLevel)
    
    options = parser.parse_args()
    
    # Parse the log level
    numeric_level = getattr(logging, options.logLevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % options.logLevel)
    
    # Setup a logger
    logger.setLevel(numeric_level)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(numeric_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.debug("Established logger")
    
    startTime = time.time()
    
    # Parse the target sequence and read the first sequence
    with pysam.FastaFile(options.target_fasta) as targetFasta:
        targetString = targetFasta.fetch(targetFasta.references[0])
    logger.info("Parsed target string. Length: %s" % len(targetString))
    
    # Build minimizer index
    minimizerIndex = MinimizerIndexer(targetString.upper(), w=options.w, k=options.k, t=options.t)
    minmerInstances = sum(map(len, minimizerIndex.minimizerMap.values()))
    logger.info("Built minimizer index in %s seconds. #minmers: %s, #minmer instances: %s" %
                 ((time.time()-startTime), len(minimizerIndex.minimizerMap), minmerInstances))
    
    # Open the query files
    alignmentScores = [] # Array storing the alignment scores found
    with pysam.FastqFile(options.query_fastq) as queryFastq:
        # For each query string build alignment
        for query, queryIndex in zip(queryFastq, range(sys.maxsize)):  # python3: xrange to range, maxint to maxsize
            print (queryIndex) # python3
            alignment = simpleMap(targetString, minimizerIndex, query.sequence.upper(), config)
            alignmentScore = 0 if alignment is None else alignment.getMaxAlignmentScore()
            alignmentScores.append(alignmentScore)
            logger.debug("Mapped query sequence #%i, length: %s alignment_found?: %s "
                         "max_alignment_score: %s" % 
                         (queryIndex, len(query.sequence), alignment is not None, alignmentScore)) 
            # Comment this out to test on a subset
            #if queryIndex > 100:
            #    break
            
    # Print some stats
    logger.critical("Finished alignments in %s total seconds, average alignment score: %s" % 
                    (time.time()-startTime, float(sum(alignmentScores))/len(alignmentScores)))
    
if __name__ == '__main__':
    main()