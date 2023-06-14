#-------------------------------------------------------------------------
#  function to be minimised
#-------------------------------------------------------------------------
import settings
def getrhs(n_,x):
  y_ = steporbit(settings.ndts,x) 
  y  = y_ - x 				# difference
  y[0] = 0E0 					# constraints, rhs=0
  return y


