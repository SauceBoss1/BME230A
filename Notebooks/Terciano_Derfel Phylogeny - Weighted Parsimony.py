# Let's do weighted small parsimony 

## Problem 1

import sys

class N: 
  """ Class to represent internal node, leaves are strings """
  def __init__(self, left, right, leftBranchLength=1.0, rightBranchLength=1.0):
    self.left, self.right, self.leftBranchLength, self.rightBranchLength = left, right, leftBranchLength, rightBranchLength
    
tree = N(N("A", N("T", "T")), N(N(N("T", "T"), "A"), "A")) # Example tree

def subCost(ancestorChar, descendantChar, branchLength):
  """ Substitution cost function """
  return 0 if ancestorChar == descendantChar else 1

# positive infinity
p_inf = float("inf")
def parsimonyCost(t, alphabet="ACGT", subCostFn=subCost):
  """ Calculates the cost of substitutions for the given tree t node of a tree, 
  returns dictionary of alphabet characters to costs"""
  # Code to write - hint use isinstance function to determine if node is internal or leaf

  if isinstance(t.left, N):
    left = parsimonyCost(t.left, alphabet, subCostFn)
  else:
    left = {a:(0 if a == t.left else p_inf) for a in alphabet} 

  if isinstance(t.right, N):
    right = parsimonyCost(t.right, alphabet, subCostFn)
  else:
    right = {a:(0 if a == t.right else p_inf) for a in alphabet}

  #print(f'left = {left} right = {right}')
  
  cost = {}
  for a in alphabet:
    cost[a] = min([x + subCostFn(a, b, t.leftBranchLength) for b, x in left.items()]) + min([x + subCostFn(a, b, t.rightBranchLength) for b, x in right.items()])
  
  return cost

                                        
print(parsimonyCost(tree)) # Should print {'A': 2, 'C': 4, 'G': 4, 'T': 3}