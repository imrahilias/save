import sys, string, socket
import math as ma
import numpy as np

host = "localhost"
port = 56700
message = ""

typs = [0, 1, 2, 3, 4, 5, 6, 71, 10, 11, 12, 13, 80, 20, 21, 22, 23, 100, 90, 70];
batch = 1024;
sample = {}
m = {}
sigma = {}
ps = {}
cnt = {}
eps = 1e-8
es = {}
emax = ma.log(6,2)
ediff = {}
em = 0
emt = {}

for typ in typs:
   
   # open socket
   try:
      s = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
   except socket.error, message:
      s = None
   try:
      s.connect((host,port))
   except socket.error, message:
      s.close()
      s = None
   if s is None:
      print 'Can not open socket ',message
      sys.exit(1)

   # send data
   message = "throw" + " " + str(batch) + " " + str(typ) + "\n"
   s.send( message )
   #print( "typ" + " " + str(typ) + " " + str(batch) + " bytes sent")
   
   # read answer, no buffer, we now size & decode
   block = s.recv( batch ).decode()
   #print( "typ" + " " + str(typ) + " " +  str(batch) + " bytes recieved")
   s.close()

   # make a nice array of the blockdata
   vec = []
   n = [0]*6
   p = [0]*6
   for i in range(len(block)):
        a = int(block[i])
        vec.append(a)
        n[a-1] += 1 # count faces, careful python is code not math!
  
   # do lowlevel statistics
   m[typ] = np.mean(vec)
   sigma[typ] = np.std(vec)
   cnt[typ] = n # save cause why not

   # its all about entropy
   e = 0
   for i in range(len(p)):
      p[i] = float(n[i]) / float(batch) # probability := count / throws per face
      if abs(p[i]) > eps:
         e -= p[i]*ma.log(p[i],2)
      
   ediff[typ] =  emax - e 
   if ediff[typ] > em:
      em = ediff[typ]
      emt = typ
   
   es[typ] = e
   ps[typ] = p

print( "and the winner of entropy maximum is: typ " + str(emt) )
print( "with delta entropy:  " + str(ediff[emt]) )
print( "with mean: " + str(m[emt]) + " +/- " + str(sigma[emt]) )

