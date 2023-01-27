# Problem 2: Build a simple suffix array
text = 'ACTTGGAGATCTTTGAGGCTAGGTATTCGGGATCGAAGCTCATTTCGGGGATCGATTACGATATGGTGGGTATTCGGGA'
pattern = 'GGTATTCGGGA'
K = 3


class SuffixArray(object):
    def __init__(self, t):
      ''' Create suffix array representing suffixes in t '''
      
      self.td = t + "$"
      self.index = [] ## Array of integers representing lexicographically sorted suffixes of t
      # e.g. for t$ = ATA$
      # have suffixes
      # 0 = ATA$
      # 1 = TA$
      # 2 = A$
      # 3 = $
      # such that self.index == [ 3, 2, 0, 1 ]
      
      # Code to complete - finish building self.index for t
      tempArr = [(self.td, 0)]
      for i in range(1, len(self.td)):
        tempArr.append((self.td[i:], i))
      #print(sorted(tempArr,key=lambda x:x[0]))
      self.index = [i[1] for i in sorted(tempArr,key=lambda x:x[0])]
      
    def query(self, p):
      ''' Return occurrences of pattern p in t'''
      
      # Code to complete - find all occurrences of p in t by writing binary search
      # function on self.index
      occurrences = set()
      l = 0
      r = len(self.index) - 1
      while l < r:
         mid = int((l+r)/2)

         if p == self.td[self.index[mid]:self.index[mid]+len(p)]:
            occurrences.add(self.index[mid])

         if p > self.td[self.index[mid]:]:
            l = mid + 1
         else:
            r = mid
      
      s = l
      r = len(self.index) - 1
      while l < r:
         mid = int((s+r)/2)
         if p == self.td[self.index[mid]:self.index[mid]+len(p)]:
            occurrences.add(self.index[mid])
         if p < self.td[self.index[mid]:]:
            r = mid
         else:
            s = mid + 1
      
      return occurrences

      # find all suffixes fo character, move onto next char of p

# Test suffix array construction
sa = SuffixArray("ATA")
print(sa.index == [ 3, 2, 0, 1 ])

# # Test suffix array search
sa = SuffixArray(text)
#print(text[50:])
print(sa.query(pattern))
print(sorted(sa.query(pattern)) == [21, 68])