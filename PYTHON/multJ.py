#-------------------------------------------------------------------------
#  Action of Jacobian on update dx, and evaluation of constraint on update.
#  Use approximation    dF(x_n)/dx . dx = (F(x_n+eps.dx)-F(x_n))/eps
#-------------------------------------------------------------------------
import settings
from dotprd import dotprd
from getrhs import getrhs

def multJ(n_,dx):
   new_x = settings.new_x
   new_fx = settings.new_fx
   epsJ = settings.epsJ
   ndts = settings.ndts
   fixT = settings.fixT
   				# (F(x0+eps.x)-F(x0))/eps
   eps = np.sqrt(dotprd(1,dx,dx)) 
   eps = epsJ * np.sqrt(dotprd(1,new_x,new_x)) / eps 
   y = new_x + eps*dx 
   s = getrhs(n_,y) 
   y = (s - new_fx) / eps 
   if fixT:			# no extra constraint if T fixed
      y[0] = 0E0 
   else:				# constraint: dx . \dot{x} = 0				# no update in trajectory direction
      s = steporbit(1,new_x) 
      dt = new_x[0]/ndts 
      s = (s - new_x) / dt 
      y[0] = dotprd(-1,s,dx) 
   return y
 


