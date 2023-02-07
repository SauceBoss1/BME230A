import sys
import numpy
import math

"""The following uses Python to challenge you to create an algorithm for finding
matches between a set of aligned strings. Minimal familiarity with Python is 
necessary, notably list and Numpy array slicing. 

Coding by: Derfel Terciano
Collaborators: Roni Altshuler
"""

"""Problem 1.

Let X be a list of M binary strings (over the alphabet { 0, 1 }) each of length 
N. 

For integer 0<=i<=N we define an ith prefix sort as a lexicographic sort 
(here 0 precedes 1) of the set of ith prefixes: { x[:i] | x in X }.
Similarly an ith reverse prefix sort is a lexicographic sort of the set of
ith prefixes after each prefix is reversed.

Let A be an Mx(N+1) matrix such that for all 0<=i<M, 0<=j<=N, A[i,j] is the 
index in X of the ith string ordered by jth reverse prefix. To break ties 
(equal prefixes) the ordering of the strings in X is used. 

Complete code for the following function that computes A for a given X.

Here X is a Python list of Python strings. 
To represent A we use a 2D Numpy integer array.

Example:

>>> X = getRandomX() #This is in the challenge1UnitTest.py file
>>> X
['110', '000', '001', '010', '100', '001', '100'] #Binary strings, M=7 and N=3
>>> A = constructReversePrefixSortMatrix(X)
>>> A
array([[0, 1, 1, 1],
       [1, 2, 2, 4],
       [2, 3, 5, 6],
       [3, 5, 4, 3],
       [4, 0, 6, 0],
       [5, 4, 3, 2],
       [6, 6, 0, 5]])
>>> 

Hint:
Column j (0 < j <= N) of the matrix can be constructed from column j-1 and the 
symbol in each sequence at index j-1.  

Question 1: In terms of M and N what is the asymptotic cost of your algorithm?

A: O(MN) this is because you need to iterate for each column and row
"""

def constructReversePrefixSortMatrix(X):
    #Creates the Mx(N+1) matrix
   A = numpy.empty(shape=[len(X), 1 if len(X) == 0 else len(X[0])+1 ], dtype=int) 
   
   #Code to write - you're free to define extra functions 
   #(inline or outside of this function) if you like.

   #inspired by Alg 1 in Durbin's PBWT paper 

   A[0,0] = 0
   for i in range(1, len(X)): #set col 0 to increasing order
      A[i,0] = A[i-1,0]+1

   M = len(X)
   N = len(X[0])

   a_k = [i for i in range(M)] #permutation of numbers 0,...,M-1

   for k in range(1, N+1): #for each char 'k'
      a=[]
      b=[]

      for i in a_k: # for each index i based on the ordering of a_k
         if X[i][:k][::-1].startswith('0'): # 0s go in list a
            a.append(i)
         else: # 1s go in list b
            b.append(i) 
      
      a_k = a+b #combine the sorted a and b lists

      for i in range(M): #update each column in row k
         A[i,k] = a_k[i]


   return A 
   
"""Problem 2: 

Following on from the previous problem, let Y be the MxN matrix such that for 
all 0 <= i < M, 0 <= j < N, Y[i,j] = X[A[i,j]][j].

Complete the following to construct Y for X. 

Hint: You can either use your solution to constructReversePrefixSortMatrix() 
or adapt the code from that algorithm to create Y without using 
constructReversePrefixSortMatrix().

Question 2: In terms of M and N what is the asymptotic cost of your algorithm?

A: O(MN) This is because we also need to iterate through each column and row
"""
def constructYFromX(X):
   #Creates the MxN matrix
   Y = numpy.empty(shape=[len(X), 0 if len(X) == 0 else len(X[0]) ], dtype=int)
   
   #Code to write - you're free to define extra functions
   #(inline or outside of this function) if you like.

   A = constructReversePrefixSortMatrix(X) #create Matrix A

   M = len(X)
   N = len(X[0])

   for i in range(0,M): #iterate through each column and row
      for j in range(0,N):
         Y[i,j] = X[A[i,j]][j] #procedure describe to us in the above description of problem 2
   

   return Y


"""Problem 3.

Y is a transformation of X. Complete the following to construct X from Y, 
returning X as a list of strings as defined in problem 1.
Hint: This is the inverse of X to Y, but the code may look very similar.

Question 3a: In terms of M and N what is the asymptotic cost of your algorithm?
A: O(MN) This is because the algorithm itself needs to iterate through each index (i) and it needs to
   iterate through each string char k.

Question 3b: What could you use the transformation of Y for? 
Hint: consider the BWT.
A: O(MN) This transformation can be used like a BWT transformation. With this algorithm, we
   can see similatries as to how we can approach this problem to when we did the bwt notebook.
   Here, we basically move the 'subtrings' around until we got our original strings just like
   in the bwt 

Question 3c: Can you come up with a more efficient data structure for storing Y?
A: Since it was mentioned in Durbins positional BWT algorithm paper, it is best if we use
   Huffman encoding since Huffman encoding uses primarily 0s and 1s to decide the order
   of the trees. 
   With that being said, it seems like the most efficient data structure for storing Y
   would be using a series heaps and trees.
"""
def constructXFromY(Y):
   #Creates the MxN matrix
   X = numpy.empty(shape=[len(Y), 0 if len(Y) == 0 else len(Y[0]) ], dtype=int)
   
   #Code to write - you're free to define extra functions
   #(inline or outside of this function) if you like.

   #inspired by Roni's thought process and by the BWT algorithm
   
   M = len(Y)
   N = len(Y[0])

   a_k = [([], i) for i in range(M)] # just like in problem 1, a modified
                                     # permutation of 0,...,M-1
                                     # but this time each number is associated with
                                     # a list of 0s and 1s that gets updated as we 
                                     # go through the algorithm

                                     # in the end each list should contain the ordered
                                     # string that is an element of X

   for k in range(N): #go through each row k
      a=[]
      b=[]
      for i in range(M): #go through each col i
         a_k[i][0].append(Y[i,k]) #append the char from Y[i,k] from each item in a_k

         if a_k[i][0][-1] == 0: #like in problem 1, use the pbwt sorting algorithm
            a.append(a_k[i])
         else:
            b.append(a_k[i])
         
      a_k = a + b #combine lists like in problem 1

   for i in range(N):
      for j in a_k:
         X[j[1],i] = j[0][i] #update matrix X (update is O(MN))
   
   return list(map(lambda i : "".join(map(str, i)), X)) #Convert back to a list of strings

"""Problem 4.

Define the common suffix of two strings to be the maximum length suffix shared 
by both strings, e.g. for "10110" and "10010" the common suffix is "10" because 
both end with "10" but not both "110" or both "010". 

Let D be a Mx(N+1) Numpy integer array such that for all 1<=i<M, 1<=j<=N, 
D[i,j] is the length of the common suffix between the substrings X[A[i,j]][:j] 
and X[A[i-1,j]][:j].  

Complete code for the following function that computes D for a given A.

Example:

>>> X = getRandomX()
>>> X
['110', '000', '001', '010', '100', '001', '100']
>>> A = constructReversePrefixSortMatrix(X)
>>> A
array([[0, 1, 1, 1],
       [1, 2, 2, 4],
       [2, 3, 5, 6],
       [3, 5, 4, 3],
       [4, 0, 6, 0],
       [5, 4, 3, 2],
       [6, 6, 0, 5]])
>>> D = constructCommonSuffixMatrix(A, X)
>>> D
array([[0, 0, 0, 0],
       [0, 1, 2, 2],
       [0, 1, 2, 3],
       [0, 1, 1, 1],
       [0, 0, 2, 2],
       [0, 1, 0, 0],
       [0, 1, 1, 3]])

Hints: 

As before, column j (0 < j <= N) of the matrix can be constructed from column j-1 
and thesymbol in each sequence at index j-1.

For an efficient algorithm consider that the length of the common suffix 
between X[A[i,j]][:j] and X[A[i-k,j]][:j], for all 0<k<=i is 
min(D[i-k+1,j], D[i-k+2,j], ..., D[i,j]).

Question 4: In terms of M and N what is the asymptotic cost of your algorithm?
A: O(MN) This is because we need to iterate through both the columns and rows
"""
#A4: O(MN)

def constructCommonSuffixMatrix(A, X):
   D = numpy.zeros(shape=A.shape, dtype=int) #Creates the Mx(N+1) D matrix 

   #Code to write - you're free to define extra functions 
   #(inline or outside of this function) if you like.

   M = len(X)
   N = len(X[0])
   
   #note row 0 and col 0 must be all 0

   for j in range(1,N+1): #each row j
      for i in range(1,M): #each col i
         str_1, str_2 = X[A[i,j]][:j][::-1], X[A[i-1,j]][:j][::-1]

         # just like descrbed in the problem above, on string
         # contains: X[A[i,j]][:j], while another contains: X[A[i-1,j]][:j][::-1]
         # NOTE: these strings are reversed so it is easier to iterate and compare through each
         #       string from the beginning. Once the strings are not equal, then we can simply break
         #       this avoids me from calculating the length of the string and iterating backwards


         str_len = 0 #intialize count
         assert len(str_1) == len(str_2) #make sure that both strings are of equal length
         while str_len < len(str_1) and str_1[str_len] == str_2[str_len]: # keep looping until the index of each string does not match
            str_len += 1 #idea inspired by grading script

         D[i,j] = str_len

   return D

"""Problem 5.
    
For a pair of strings X[x], X[y], a long match ending at j is a common substring
of X[x] and X[y] that ends at j (so that X[x][j] != X[y][j] or j == N) that is longer
than a threshold 'minLength'. E.g. for strings "0010100" and "1110111" and length
threshold 2 (or 3) there is a long match "101" ending at 5.
    
The following algorithm enumerates for all long matches between all substrings of
X, except for simplicity those long matches that are not terminated at
the end of the strings.
    
Question 5a: What is the asymptotic cost of the algorithm in terms of M, N and the
number of long matches?
A: Based on Durbin's position BWT, the algorithm is examnined to be of 
   O(max(NM, number of matches)). This means that the runtime is 
   proportional to the maximum number of matches for each row and column.
   However, when return the results of the matches, the runtime of that would
   be O(NM)

Question 5b: Can you see any major time efficiencies that could be gained by
refactoring?
A: In this algorithm, I can see that if we combine the loops that iterate through 
   list b and c into 2 for loops, then it would make the algorithm run so much farther.
   In addition, if we used a different approach to where we would only need one matrix
   instead of dealing with 2 extra for loops, we can then iterate through the entirety
   of the matrix once.
    
Question 5c: Can you see any major space efficiencies that could be gained by
refactoring?
A: A space efficiency I found would be to find a way to combine matrices A and
   D together into one single operation. This is because we need to create a 
   matrix of O(MN) and then create another matrix of O(MN) and this would take up
   extra space. 
    
Question 5d: Can you imagine alternative algorithms to compute such matches?,
if so, what would be the asymptotic cost and space usage?
A: Based on Durbin's PBWT paper, we can see that in algorithm 4, it has a time complexity of 
   O(MN) while having a space usage of O(M). Algorithm 4 computes the set maximal matches in
   string k whereas getLongMatches() gets only the matches based on a specific threshold.
"""
def getLongMatches(X, minLength):
   assert minLength > 0
   
   A = constructReversePrefixSortMatrix(X)
   D = constructCommonSuffixMatrix(A, X)
   
   #For each column, in ascending order of column index
   for j in range(1, 0 if len(X) == 0 else len(X[0])):
      #Working arrays used to store indices of strings containing long matches
      #b is an array of strings that have a '0' at position j
      #c is an array of strings that have a '1' at position j
      #When reporting long matches we'll report all pairs of indices in b X c,
      #as these are the long matches that end at j.
      b, c = [], []
      
      #Iterate over the aligned symbols in column j in reverse prefix order
      for i in range(len(X)):
         #For each string in the order check if there is a long match between
         #it and the previous string.
         #If there isn't a long match then this implies that there can
         #be no long matches ending at j between sequences indices in A[:i,j]
         #and sequence indices in A[i:,j], thus we report all long matches
         #found so far and empty the arrays storing long matches.
         if D[i,j] < minLength:
               for x in b:
                  for y in c:
                     #The yield keyword converts the function into a
                     #generator - alternatively we could just to append to
                     #a list and return the list
                     
                     #We return the match as tuple of two sequence
                     #indices (ordered by order in X) and coordinate at which
                     #the match ends
                     yield (x, y, j) if x < y else (y, x, j)
               b, c = [], []
         
         #Partition the sequences by if they have '0' or '1' at position j.
         if X[A[i,j]][j] == '0':
               b.append(A[i,j])
         else:
               c.append(A[i,j])
      
      #Report any leftover long matches for the column
      for x in b:
         for y in c:
               yield (x, y, j) if x < y else (y, x, j)