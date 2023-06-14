#-------------------------------------------------------------------------
#  advance x by ndts_ timesteps
#-------------------------------------------------------------------------
from Lorenz_f import Lorenz_f
import settings

def steporbit(ndts_,x):
   if ndts_ != 1:		# Set timestep size dt=T/ndts_
      settings.dt = x[0] / ndts_ 	# If only doing one step to calc \dot{x},
   				# then use previously set dt.

   a = x[1:] 
         # second-order predictor-corrector method
   for n in range(ndts_):
      fa  = Lorenz_f(a) 
      a1  = a + settings.dt*fa 
      fa1 = Lorenz_f(a1) 
      a  = a + 0.5E0*settings.dt*(fa+fa1) 

   y = x.copy()
   y[1:] = a 
   return y
 

