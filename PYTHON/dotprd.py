#-------------------------------------------------------------------------
# dot product.  n_==-1 is a flag to exclude parameter T.
#-------------------------------------------------------------------------
def dotprd(n_, a, b):
  n1 = 1
  if n_ == -1:
    n1 = 2
  d = sum(a[n1:] * b[n1:])
  return d


