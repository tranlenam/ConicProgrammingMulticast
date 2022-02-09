import numpy as np
import sys
from mosek.fusion import *
#from scipy.sparse import coo_matrix
import matplotlib.pyplot as plt
# Define a stream printer to grab output from MOSEK
def complex2realchannel(complexchannel):
    channelreal =  np.zeros([nUsers,2,2*nTx])
    for iUser in range(nUsers):
        b = np.concatenate((complexchannel[iUser,:].real,-complexchannel[iUser,:].imag),axis=0)
        c= np.concatenate((complexchannel[iUser,:].imag,complexchannel[iUser,:].real),axis=0)
        #u = np.vstack((b,c))
        channelreal[iUser,:,:] = np.vstack((b,c))
    return channelreal 
# build the convex model for (6)
def initializeModel(nTx, nUsers):
    M = Model()
    t = M.variable('t', 1, Domain.unbounded())
    w = M.variable('w', 2*nTx, Domain.unbounded())
    M.objective("obj", ObjectiveSense.Minimize, t)
    z = Var.vstack(t, w)
    M.constraint("qc", z, Domain.inQCone())
    b = M.parameter("b",nUsers)
    A = M.parameter("A",nUsers,2*nTx)
    M.constraint(Expr.sub(Expr.mul(A,w),b), Domain.greaterThan(0.0))
    # Return the ready model
    return M


np.random.seed(0) 
nTx = 100 # number of Tx antennas

nUsers = 20 # Number of users

nIterations = 1  # Maximum number of iterations

nWarmups = 1 # number of warmups

totaltime = 0 # total run time of the SCA method

nRuns = 10 # total number of runs
inf = 0.0
nIters = 100; # maximum number of iterations   
errtol =1e-2   
# Make a MOSEK environment


M = initializeModel(nTx, nUsers)

  
b = M.getParameter("b")
A = M.getParameter("A")
w = M.getVariable("w")

for iRun in range(nRuns):
    objseq =[0]*nIters # objective sequence
    channel = 0.5 ** 0.5 * (np.random.randn(nUsers,nTx)+ 1j * np.random.randn(nUsers,nTx))
    wnextwarmup = np.random.randn(nTx,nWarmups) + 1j * np.random.randn(nTx,nWarmups)
    wnextwarmup = wnextwarmup/np.reshape(np.absolute(np.matmul(channel,wnextwarmup)).min(axis=0),(1,nWarmups))
    idx = (np.linalg.norm(wnextwarmup,axis=0)).argmin(axis=0)
    wnext = wnextwarmup[:,idx,None]            
    wnextreal = np.concatenate((wnext.real,wnext.imag),axis=0)
    channelreal = complex2realchannel(channel)

    for iIter in range(nIters):
        a = np.matmul(np.matmul(wnextreal.T,channelreal.transpose(0,2,1)),channelreal)
        a = np.squeeze(a)
     
        b.setValue((np.matmul(a,wnextreal) + 1).squeeze())
        A.setValue(2*a)
        M.solve()
        wnextreal = np.array(w.level()).reshape(2*nTx,1)
        objseq[iIter]=M.primalObjValue()
        if (iIter>9) and (abs(objseq[iIter]-objseq[iIter-9])<errtol):
            objseq = np.delete(objseq,slice(min(iIter+1,nIters),nIters))
            break
            
M.dispose() # clear the model
# plot the convergence of the objective sequence for the last run            
plt.plot(range(len(objseq)), objseq)
plt.show()
    
    
