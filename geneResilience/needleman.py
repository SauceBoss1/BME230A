import numpy
from itertools import product
from collections import deque

class NeedlemanWunsch(object):
    def __init__(self, string1, string2, gapScore=-2, matchScore=3, 
                 mismatchScore=-3):        
        self.editMatrix = numpy.zeros(shape=[len(string1)+1, len(string2)+1], dtype=int) # Numpy matrix representing edit matrix
        self.string1 = string1
        self.string2 = string2
        self.gapScore = gapScore
        self.matchScore = matchScore
        self.mismatchScore = mismatchScore

        
        self.instantiateMatrixLines()
        self.fillEditMatrix()
        #print('NW:\n',self.editMatrix)
        
        self.comparedList = self.getAlignment()
        self.score = self.editMatrix[-1][-1]
        

    def instantiateMatrixLines(self):
        for dist in range(len(self.editMatrix[0])):
            self.editMatrix[0][dist] = (dist*self.gapScore)
            
        for dist in range(len(self.editMatrix)):
            self.editMatrix[dist][0] = (dist*self.gapScore)

                
    def fillEditMatrix(self):
        for string1char in range(1,len(self.string1)+1):
            for string2char in range(1,len(self.string2)+1):
                
                s1charString = self.string1[string1char-1]
                s2charString = self.string2[string2char-1]
                
                topPlus = self.editMatrix[string1char-1][string2char]+self.gapScore
                leftPlus = self.editMatrix[string1char][string2char-1]+self.gapScore
                
                if s1charString != s2charString:
                    diag = self.editMatrix[string1char-1][string2char-1] + self.mismatchScore
                else:
                    diag = self.editMatrix[string1char-1][string2char-1] + self.matchScore
                
                if (diag>=leftPlus)and(diag>=topPlus):
                    newFill = diag
                elif(leftPlus>=topPlus):
                    newFill = leftPlus
                else:
                    newFill = topPlus
                self.editMatrix[string1char][string2char] = newFill
                                   
                
    def getAlignment(self):
        """
        method based on of "Parallel Needleman-Wunsch Algorithm for Grid" and John Lekberg's blog
        
        
        Returns an optimal global alignment of two strings. Aligned
        is returned as an ordered list of aligned pairs.
        
        e.g. For the two strings GATTACA and TACA an global alignment is
        is GATTACA
           ---TACA
        This alignment would be returned as:
        
        [(3, 0), (4, 1), (5, 2), (6, 3)]
        """
        string1=self.string1
        string2=self.string2
        string1Len, string2Len = len(string1), len(string2)
        subStrs = lambda a, b: int(a == b)
        F = {}
        Ptr = {}
    
        F[-1, -1] = 0
        for i in range(string1Len):
            F[i, -1] = -i
        for j in range(string2Len):
            F[-1, j] = -j
    
        DIAG = -1, -1
        LEFT = -1, 0
        UP = 0, -1
        option_Ptr = DIAG, LEFT, UP
        for i, j in product(range(string1Len), range(string2Len)):
            option_F = (
                F[i - 1, j - 1] + subStrs(string1[i], string2[j]),
                F[i - 1, j] - 1,
                F[i, j - 1] - 1,
            )
            F[i, j], Ptr[i, j] = max(zip(option_F, option_Ptr))
    
        alignment = deque()
        i, j = string1Len - 1, string2Len - 1
        while i >= 0 and j >= 0:
            direction = Ptr[i, j]
            if direction == DIAG:
                element = i, j
            elif direction == LEFT:
                element = i, None
            elif direction == UP:
                element = None, j
            alignment.appendleft(element)
            di, dj = direction
            i, j = i + di, j + dj
        while i >= 0:
            alignment.appendleft((i, None))
            i -= 1
        while j >= 0:
            alignment.appendleft((None, j))
            j -= 1
        return list(alignment)
    
    
    def show_stringA(self):
        new_string = ''
        for tup in self.comparedList:
            if tup[0] != None:
                new_string += self.string1[tup[0]]
            else:
                new_string+='_'
        return new_string
                
    def show_stringB(self):
        new_string = ''
        for tup in self.comparedList:
            if tup[1] != None:
                new_string += self.string2[tup[1]]
            else:
                new_string+='_'
        return new_string
      
#string1 = "TTTT"
#string2 =    "TTATT"

#needlemanWunsch = NeedlemanWunsch(string1, string2)
#needlemanWunsch.show_stringA()
#needlemanWunsch.show_stringB()