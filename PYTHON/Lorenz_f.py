#-------------------------------------------------------------------------
#  Lorenz system: dx/dt = f(x)
#            state:  x[0]==X, x[1]==Y, x[2]==Z.
#       parameters:  p[0]==s, p[1]==r, p[2]==b  (classic: s=10,r=28,b=8/3)
#-------------------------------------------------------------------------
import numpy as np
import settings 

def Lorenz_f(x):
   p = settings.p
   dx = np.zeros(x.shape) 
   dx[0] = p[0]*(-x[0]+x[1]) 
   dx[1] = x[0]*(-x[2]+p[1]) - x[1] 
   dx[2] = x[0]*x[1] - p[2]*x[2] 
   return dx


