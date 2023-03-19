'''
geneFinder.py

by Derfel Terciano version 0.1

This library takes in sequences and finds open reading frames using regular expressions.
'''
from collections import OrderedDict
import re
import sys


class Orfs:
  '''
  Orfs

  Finds Open Reading Frames (ORFs) and stores it this object

  Inputs:
    - orfs: tuple -> this tuple stores the start and end index of where the ORF was found
    - seq: str -> stores the orf's corresponding sequence
    - length: int -> stores the length of the found orf

  public functions:
    - getOrfIndex (property) -> returns the indices where the ORF was found
    - getSeq (property) -> returns the sequence of the ORF
    - getLength (property) -> returns the length of the ORF
    - findOrfs (static method) -> finds all ORFs of a given seq and minimum length and
                                  returns a list of Orf objects
  '''
  def __init__(self, orfs: tuple, seq: str, length: int):
    self.orfs = orfs
    self.seq = seq
    self.length = length

  @property
  def getOrfIndex(self):
    '''Return the index of where the Orf was found'''
    return self.orfs

  @property
  def getSeq(self):
    '''Return the ORF sequence'''
    return self.seq
  
  @property
  def getLength(self):
    '''Return the length of the sequence'''
    return self.length

  @staticmethod
  def findOrfs(dnaseq: str, minLength=100):
    '''
    Find all ORFs given a dna sequence

    inputs:
      - danaseq: str -> sequence to find ORFs from
      - minLength -> minimum length size of an ORF (default: 100)
    
    output:
      - Returns a list of ORF objects
    '''
    orfs = [] # orfs will contain an object of ORFs which contain index, length and sequence
    #dnaseq = dnaseq.replace('_', '')
    start_codons = 'A.?TG|A.?_TG|GTG|TTG|ATG' #regex expression that takes into account some mut. starts
    #start_codons = 'ATG' #regex expression that takes into account some mut. starts
    stop_codons = 'TAA|TAG|TGA' #regex expression that takes into account stop codons (non mutated)

    for frame in range(3):
      start_pos = [(m.start(), m.end()) for m in re.finditer(start_codons, dnaseq[frame:])] # find the indices of ALL start codons

      for starts in start_pos:
        curr_start_frame = dnaseq[starts[0]+frame:]
        for stops in re.finditer(stop_codons, curr_start_frame):
          #print(stops)
          stop_pos = stops.start() + starts[0] + 3 + frame

          if stop_pos > starts[1] and stop_pos > (starts[0] + frame) and (stop_pos - (starts[0]+frame)) % 3 == 0 :

            if (stop_pos - starts[0] + frame) >= minLength:
              orfs.append(Orfs((starts[0]+frame, stop_pos), curr_start_frame[:stops.end()], (stop_pos - starts[0]+frame)))
            break

    return orfs

class AlignAllOrfs:
  '''
  AlignOrfs

  Given a list of multiple sequences, find all orfs within those sequences and
  align them.

  inputs:
    - seq: list -> list of sequences to find Orfs
    - minLength -> the minimum ORF length to find (default 100)
  
  public functions:
    - getOrfs (property) -> returns a list of list of Orf objects that are unaligned
    - printAlignedOrfs -> prints all Orfs found but in an aligned format

  '''
  def __init__(self, data: list, minLength = 100):
    self.seqs = data
    self.minLength = minLength

    self.foundOrfs = []
    for seq in self.seqs: # for each sequences, find all ORFs within that sequence
      self.foundOrfs.append((seq[0], Orfs.findOrfs(seq[1], minLength=self.minLength)))
    #print(self.seqs)
  
  @property
  def getOrfs(self):
    '''Return unaligned orfs'''
    return self.foundOrfs

  def printAlignedOrfs(self, printStrings = True, outFile = sys.stdout, spacer = '-', ultraAlign = True, strLim = 500):
    '''
    Prints all aligned orfs

    inputs:
      printStrings -> optional parameter to enable printing (default -> True)
      outFile -> what file to print alignments out (default -> sys.stdout)
      spacer -> parameter to manipulate spacing chars (default -> '-')

    '''
    finalAlignment = []
    orfsOnly = [orf[1] for orf in self.foundOrfs]

    for orfs in orfsOnly:

      finalString = ''
      lastDist = 0
      ele_num = 0

      for orf in orfs:
        currDist = orf.getOrfIndex

        if ele_num == 0: # if this is the first element, then no spacing alignment needs to be done
          finalString += spacer * (currDist[0])
          ele_num += 1
        else:
          finalString += spacer * (currDist[0] - lastDist[1]) # calculate how much space is needed
                                                              # this is done by taking into account
                                                              # the current starting postion
                                                              # and the last end position found
        
        lastDist = currDist

        if ultraAlign:
          try:
            if finalString[-1] == spacer:
              finalString += orf.getSeq
            else:
              finalString += '\n'+ (spacer * currDist[0]) + orf.getSeq
          except:
            finalString += orf.getSeq
        else:
          finalString += orf.getSeq

      if finalString != '':
        finalAlignment.append(finalString[:strLim])

      if printStrings and finalString != '':
        print(finalString[:strLim] + '\n', file=outFile)
    return finalAlignment

def main():
  testString = ['A_TGAA__AA_CCCTAA_CC__GG_GATTTT__TTTAG__', \
                'A_TGAA__AA_CCTGAA_CC__GG_GATGTT__TT_AGT_', \
                'A_TGAA__AA_CC_TAAGCC__GG_GATGTT__TT_AGT_', \
                'A_TGAA__AA_CC_TAA_CC__GG_GATGTTG_T__AGTT', \
                'ACTGAA__AA_CC_TAA_CC__GG_GATGTTG_TTTAG__', \
                'A_TGACA_AA_CC_TAA_CC__GG_GATGTT__TT_AGT_', \
                'A_TGACA_AA_CC_TAA_CC__GG_GATGTTGGTTTAA__', \
                'A_TGAAA_AC_CC_TAA_CC__GG_GATGTT__TTTAG__', \
                'A_TGA__AAA_CC_TACACC__GG_GATGTT__TT_AGT_', \
                'A_TGAA__AA_CC_TA_ACC__GGCGATGTT__TTTAG__'
               ]
  
  test2 = ['ATGAATT_A_TGAGTATAAATT_A_TGAGTATAAAAA_TAAA_TGAGTATAAAAA_TAAAA_G_ACCGCAGAAAGCTTA_CG___AG_AA___G_TTAA']
  
  #out = open('out.out', 'w')
  out = sys.stdout
  final = AlignAllOrfs(testString, minLength=10)
  final.printAlignedOrfs(outFile=out, ultraAlign=True)
  #unal = [[y.getSeq for y in x] for x in final.foundOrfs]
  #print(unal)
  

  out.close()

if __name__ == '__main__':
  main()