"""
Needs the fasta inFile name -- inFile = "" -- located at bottom

outputs a fasta file called "alignment_output.txt"


#alignOb = alignTest(stringList, idList)
#alignOb.returnAlignments()

"""
from needleman import NeedlemanWunsch
from Bio import SeqIO
class alignTest:
    def __init__(self, sList, IDs, pipeLine = True):
        self.pipLine = pipeLine
        self.sList, self.IDs = sList, IDs
        self.setScoreMeans = []
        self.loopMains()
        self.selectBestSet()
        
        
    def compareSet(self, mainSeqIndex, needsReturn):
        thisSet = self.sList
        initMain = thisSet[mainSeqIndex]
        initLast = thisSet[mainSeqIndex]
        done=False
        loopCount = 0
        while not done:
            setIDs = self.IDs.copy()
            setScores = []
            
            for i in range(len(self.sList)):
                if i == 0:
                    initMain = initLast
                    aligned = [initMain]
                if i != mainSeqIndex:
                    thisNM = NeedlemanWunsch(initLast,thisSet[i])
                    revNM = NeedlemanWunsch(initLast,self.revComp(thisSet[i]))
                    
                    if revNM.score > thisNM.score:
                        thisNM = revNM
                        setIDs[i]+='-rev'
                    
                    setScores.append(thisNM.score)
                    initLast = thisNM.show_stringA()
                    print(initLast, self.IDs[mainSeqIndex])
                    adjString2 = thisNM.show_stringB()
                    aligned.append(adjString2)
                    print(adjString2, self.IDs[i])
                    #print('\n')    
                    
                        
            if initMain == initLast:
                done = True
            #print('\n\n')
            loopCount+=1
        #print('loopCount: '+str(loopCount))
        self.setScoreMeans.append(float(sum(setScores))/float(len(setScores)))
        if needsReturn:
            return aligned, setIDs
        
    def loopMains(self):
        for i in range(len(self.sList)):
            self.compareSet(i, False)
            
    def selectBestSet(self):
        bestScore = 0
        bestIndex = -1
        for i in range(len(self.setScoreMeans)):
            score = self.setScoreMeans[i]
            if score > bestScore:
                bestScore = score
                bestIndex = i
        self.bestIndex = bestIndex
        print('Main = '+self.IDs[bestIndex], bestScore)
        
    def returnAlignments(self):
        alignments, IDs = self.compareSet(self.bestIndex, True)
        if not self.pipLine:
            with open('alignment_output.txt', 'w') as f:
                for i in range(len(self.IDs)):
                    print(alignments[i],IDs[i])
                    f.write('>'+IDs[i]+'\n')
                    f.write(alignments[i]+'\n\n')
        else:
            return alignments
                
            
    def revComp(self, seq):
        newString = ''
        for i in range(len(seq)-1, 0, -1):
            base = seq[i]
            if base == '_':
                newString +='_'
            if base == 'A':
                newString+='T'
            if base == 'C':
                newString+='G'
            if base == 'G':
                newString+='C'
            if base == 'T':
                newString+='A'
        return newString
        
def main(inFile):
    idList = []
    stringList = []
    lineCount = 0
    #inFile = open(inFile, 'r')
    
    # for line in inFile:
    #     line.strip()
    #     if line[0] == '>':
    #         end = line.index('\n')
    #         idList.append(line[1:end])
    #         lineCount+=1
    #     elif lineCount == 1:
    #         stringList.append(line.strip())
    #         lineCount = 0
    # inFile.close()

    records = list(SeqIO.parse(inFile, 'fasta'))
    stringList = [r.seq for r in records]
    idList = [r.id for r in records]
    
    alignOb = alignTest(stringList, idList, pipeLine=False)
    alignOb.returnAlignments()
    

if __name__ == '__main__':
  inFile = 'final/A23EEV21.results.fa'
  main(inFile)
    