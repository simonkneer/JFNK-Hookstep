#-------------------------------------------------------------------------
#  Global Variables
#-------------------------------------------------------------------------
def init():
   global new_x	# Current best x
   global epsJ		# epsilon used in Jacobian approximation
   global ndts		# Number of timesteps taken in period T
   global fixT		# Fix T for equilibrium, rather than PO solution
   global p		# Parameters of dynamical system
   global dt # timestep