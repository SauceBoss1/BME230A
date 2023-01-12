# Problem 1 - implement a naive string matching algorithm:

pattern = "CCTTTTGC"
#pattern = 'ACTTA'
text =    "CGTGCCTACTTACTTACTTACCTTTTGCCTTTTGCACGCGAA"

def naive(p, t):
  characterComparisons = 0
  occurrences = []

  # Code to complete - do not use find or the "in" operator
  for i in range(0, len(t) - len(p) + 1):
    matched = True
    for j in range(0, len(p)):
      characterComparisons += 1
      if t[i+j] != p[j]:
        matched = False
        break
    if matched:
      occurrences.append(i)
  return occurrences, characterComparisons

#naive(pattern, text) == ([20, 27], 60)


# Problem 2 - implement the bad character rule:

dnaAlphabet="ACGT"

def makeBadCharacterRuleLookupTable(p, alphabet):
  
  # Make len(alphabet) x len(p) lookup table
  badCharacterRuleLookupTable = [ [1]*len(alphabet) for i in range(len(p)) ]
  
  # Code to complete

  for i in range(1,len(p)):
    currChar = p[i-1]
    currIt = badCharacterRuleLookupTable[i]
    ind = alphabet.index(currChar)
    for j in range(0,len(currIt)):
      badCharacterRuleLookupTable[i][j] = badCharacterRuleLookupTable[i-1][j] + 1
    badCharacterRuleLookupTable[i][ind] = 1
      
  
  return badCharacterRuleLookupTable
      
          
# makeBadCharacterRuleLookupTable(pattern, dnaAlphabet) == [[1, 1, 1, 1],
#  [2, 1, 2, 2],
#  [3, 1, 3, 3],
#  [4, 2, 4, 1],
#  [5, 3, 5, 1],
#  [6, 4, 6, 1],
#  [7, 5, 7, 1],
#  [8, 6, 1, 2]]

 # Problem 3 - use the bad character rule lookup table to reduce the total number of character comparisons

def naivePlusBadCharacter(p, t, bclut, alphabet):
   characterComparisons = 0
   occurrences = []
  
   # Code to complete
   lastIndex = len(p)-1

   characterComparisons += 1
   while lastIndex <= len(t) - 1:
      if lastIndex > len(t) - 1:
         break
      
      currT = t[lastIndex - (len(p) - 1): lastIndex + 1]
      match = True
      for i in range(len(p) - 1, -1, -1):
         characterComparisons += 1
         indexCounter = lastIndex - ((len(p) - 1) - i) #mathematically decrements index from string
         if t[indexCounter] != p[i]:
            lastIndex = bclut[i][alphabet.index(t[indexCounter])] + lastIndex
            #print(bclut[i][alphabet.index(t[indexCounter])])
            match = False
            break

      if match:
         characterComparisons += 1
         occurrences.append(lastIndex - (len(p) - 1))
         lastIndex += (len(p) - 1)

   return occurrences, characterComparisons

bclut = makeBadCharacterRuleLookupTable(pattern, dnaAlphabet)
naivePlusBadCharacter(pattern, text, bclut, dnaAlphabet) == ([20, 27], 24)
print(naivePlusBadCharacter(pattern, text, bclut, dnaAlphabet))

# pattern = "                    CCTTTTGC"
# pattern = 'ACTTA'
# text =    "CGTGCCTACTTACTTACTTACCTTTTGCCTTTTGCACGCGAA"