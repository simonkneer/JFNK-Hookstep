#----------------------------------------------------------------------b()
# https://github.com/apwillis1/JFNK-Hookstep
#
# Please cite either of
# https://doi.org/10.48550/arXiv.1908.06730
# https://doi.org/10.1016/j.softx.2017.05.003 
#                                      Thanks in advance! Ashley 2023.
#----------------------------------------------------------------------
# solve A x = b for x   
# minimise |Ax-b| subject to constraint |x| < delta .
# requires lapack routines dgelsy, dgesvd.
#----------------------------------------------------------------------
# m	  gmres dimension
# n 	  dimension of x
# x	  on input:  guess for x, can be 0
#         on exit:  solution x, subject to constraint if del>0
# b	  input b
# matvec  performs y = A x:  y = matvec(N,x)
# psolve  preconditioner, solve M y = x:  y = psolve(N,x)
# dotprd  dot product, d = dotprd(n,a,b)
# h       Hessian matrix,  size (m+1)*m
# v       Krylov subspace, size n*(m+1)
# res	  on input: |Ax-b|/|b|<res 
#         on exit:  residual reached
# del     on input: if(del>0) then the x returned is the hookstep
#         on exit:  norm of next b predicted by hook
# its	  on input: max num its 
#         on exit:  number of its taken
# info	  on input: if(info==1) print* residuals
#                   if(info==2) recalc hookstep with new del
# 	  on exit:  0 sucessful, 1 method breakdown, 2 max its
#----------------------------------------------------------------------
from dotprd import dotprd

def GMRESm(m,n,x,b,matvec,psolve,dotprd,res,del,its,info):
   persistent h
   persistent v
   persistent beta_
   persistent j_
   """
   if info == 2  
      [y,del] = GMREShook(j_,h,m,beta_,del) 
      z = v(:,1:j_)*y(1:j_) 
      x = psolve(n,z) 
      info = 0 
      return
   end	 
   """
   tol = res
   imx = its
   its = 0 
   v   = zeros(n,m+1) 

   while True: #(restart)
      res_ = 1E99 
      stgn = 1E0 - 1E-14 

      beta_ = np.sqrt(dotprd(n,x,x)) 
      if beta_ == 0E0: 
         w = 0E0 
      else:
         w = matvec(n,x) 
       
      w = b - w 
      beta_ = np.sqrt(dotprd(n,w,w)) 
      v(:,1) = w / beta_ 
      
      h = np.zeros([m+1,m]) 
      for j in range(m):
         j_ = j 
         its = its + 1 
         z = v[:,j]  
         z = psolve(n,z) 
         w = matvec(n,z) 
         for i in range(j):
            h[i,j] = dotprd(n,w,v[:,i]) 
            w = w - h[i,j]*v[:,i] 
         
         h[j+1,j] = np.sqrt(dotprd(n,w,w)) 
         v[:,j+1] = w / h[j+1,j] 
            
         p = np.zeros([j+1,1]) 
         p[0] = beta_ 
         h_[:j+1,:j] = h[:j+1,:j] 
         y = np.linalg.pinv(h_)*p
         
         p = - h[:j+1,:j]*y 
         p[0] = p[0] + beta_ 
         res  = np.sqrt(np.sum(p*p)) 
         if info==1: 
            print('gmresm: it={}  res={}\n'.format( its, res))  
         
         
         done = ((res <= tol) or (its == imx) or (res > res_)) 
         if done or j==m:
            if del > 0E0:  
               [y,del] = GMREShook(j,h,m,beta_,del) 
            
            z = v[:,:j]*y[:j] 
            z = psolve(n,z) 
            x = x + z 
            if its==imx:
               info = 2 
            
            if res>res_:
               info = 1 
            
            if res<=tol:
               info = 0 
             
            if done: 
               return 
            
            if del>0E0:   
               print('gmres: WARNING: m too small. restart affects hookstep.\n')
            
         
         res_ = res*stgn 

      
   return [x,res,del,its,info]



