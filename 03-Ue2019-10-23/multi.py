from multiprocessing import Pool
import time
import numpy as np
import os


def DotProduct(a,b):
    if len(a)==len(b):
        sum=0
        for i in range(0,len(a)):
            sum=sum+a[i]*b[i]
    else:
        print("Vectors are not of the same length")
    return sum


def Matrixmultiplication(A,B):
    work=[]
    mdim, dummy= np.shape(A)
    dummy,ndim= np.shape(B)
    for m in range(0,mdim):
        for n in range(0,ndim):
            a=A[m,:]
            b=B[:,n]
            Vectors=[m,n,a,b]
            work.append(Vectors)
    return work


def do_work(work):
    Result=[work[0],work[1],DotProduct(work[2],work[3])]
    return Result

def pool_handler():
    p = Pool(1)
    Calcs=p.map(do_work, work)
    return Calcs


A=np.random.rand(500,500)
B=np.random.rand(500,500)
m=np.shape(A)[0]
n=np.shape(B)[1]
C=np.empty([m,n])
work=Matrixmultiplication(A,B)
Big_Result=[]

start=time.clock()

#if __name__ == '__main__':
hope=pool_handler()
for i in hope:
    C[i[0],i[1]]=i[2]
    
end=time.clock()

print("Numpy:",abs(np.subtract(C,np.matmul(A,B))) > 1e-8 )
print("Duration is ", end-start)
