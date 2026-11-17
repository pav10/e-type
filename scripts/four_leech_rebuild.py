#!/usr/bin/env python3
"""Rebuild the 96x96 Gram matrix G of the SRFP311T1 'four-Leech graph-code' manuscript.
Data and construction copied verbatim from their verify_publication_core.py (supplied by the
user); everything after the construction (their certificate checks) is omitted.
Writes external/four_leech_G96.txt."""
import os
#!/usr/bin/env python3
"""Exact core certificate for the SRFP311T1 96-dimensional graph-code manuscript.

Independent of external matrix/data files. Requires numpy, sympy.
Checks: embedded Gram construction, Leech identification via exhaustive exact
rational norm-2 enumeration, exact graph-code identities and determinant,
mod-7 anti-isometry and module structure, original + Hadamard-dual eigenvalue
obstructions, full 96-generator Hadamard congruence, and eight new shell-28 projective-line determinant tests.

Does NOT certify remaining shells 28,42,56, zero-component hash exclusions,
or minimal-companion CVP claims. Does NOT prove extremality.
"""
import time, math, hashlib
from fractions import Fraction
import numpy as np
import sympy as sp

AP_TEXT = """84 -1 -16 -3 13 -28 -50 52 19 -34 -7 8
-34 84 -1 -55 -3 13 -7 -50 52 74 -34 -7
-33 -34 84 60 -55 -3 27 -7 -50 -5 74 -34
3 -57 17 84 -1 -16 -3 13 -28 -50 52 19
22 3 -57 -34 84 -1 -55 -3 13 -7 -50 52
-49 22 3 -33 -34 84 60 -55 -3 27 -7 -50
-34 -7 8 3 -57 17 84 -1 -16 -3 13 -28
74 -34 -7 22 3 -57 -34 84 -1 -55 -3 13
-5 74 -34 -49 22 3 -33 -34 84 60 -55 -3
-50 52 19 -34 -7 8 3 -57 17 84 -1 -16
-7 -50 52 74 -34 -7 22 3 -57 -34 84 -1
27 -7 -50 -5 74 -34 -49 22 3 -33 -34 84"""
BP_TEXT = """0 -33 -17 6 9 -21 16 22 -24 -16 0 19
33 0 -33 -2 6 9 0 16 22 -22 -16 0
17 33 0 -43 -2 6 -19 0 16 24 -22 -16
-6 2 43 0 -33 -17 6 9 -21 16 22 -24
-9 -6 2 33 0 -33 -2 6 9 0 16 22
21 -9 -6 17 33 0 -43 -2 6 -19 0 16
-16 0 19 -6 2 43 0 -33 -17 6 9 -21
-22 -16 0 -9 -6 2 33 0 -33 -2 6 9
24 -22 -16 21 -9 -6 17 33 0 -43 -2 6
16 22 -24 -16 0 19 -6 2 43 0 -33 -17
0 16 22 -22 -16 0 -9 -6 2 33 0 -33
-19 0 16 24 -22 -16 21 -9 -6 17 33 0"""
R_TEXT = """632123604334354443530236
041344364566545016610102
501625560535255036565213
364442263334255153051325
543321241204120412162650
245434632225061425425321
230140024145660244442025
430323624301331046254125
213222634664405123132034
144063030601163422421613
020233502262001455255250
246253204026031446555032
645013026555134121516113
605665564342345666353254
164315560305062130363065
060644560060144366366504
316461115335451520522645
551242011331616551034613
366052534022041665504226
115350654661054546563306
253521563553651203416330
264554153225222521036203
321330256305526400023560
142346254121220440344332"""
# Frozen 24x24 unimodular basis transformation. This need not be LLL-verified:
# only its exact determinant and its congruence P_red=U P U^T matter.
U_TEXT = """-2 -1 0 -1 1 -1 -4 1 -3 -1 1 -3 -2 1 0 0 -1 -1 0 1 0 0 0 0
2 0 2 0 1 2 4 -1 3 2 1 3 0 0 1 0 2 1 -1 0 0 0 0 0
-1 -1 -1 1 -1 -3 -4 3 -6 -1 -2 -5 0 1 -1 1 -6 0 3 1 -4 3 -3 0
-1 2 -3 2 -2 0 -1 4 -5 3 -3 0 1 2 0 2 -3 2 2 1 -2 1 0 0
2 -1 3 0 0 2 3 -1 1 1 1 1 0 -1 2 -1 1 1 0 -1 0 0 0 0
-4 1 -3 -1 2 -3 -6 2 -2 -1 0 -3 -4 4 -4 1 -2 -3 -1 3 -1 0 -2 0
2 -1 3 0 1 2 3 -3 6 0 3 2 0 -1 2 -1 4 -1 -1 -1 3 -1 0 0
4 0 1 2 -3 3 7 -3 4 2 -2 5 4 -4 2 -1 4 2 1 -4 3 0 2 0
-1 0 -3 -1 1 -2 -3 -1 0 -3 0 -1 -1 0 -4 0 0 -2 -1 0 0 -1 1 0
-5 2 -4 -4 5 -4 -6 0 2 -4 1 -2 -4 4 -7 0 -1 -4 -4 3 -1 -3 0 0
-3 -1 2 -4 4 0 -4 -1 0 -3 4 -2 -4 0 0 -2 2 -2 -5 1 1 -4 3 -1
-4 2 -5 -1 1 -3 -7 3 -3 -2 -1 -2 -1 4 -5 3 -4 -1 -1 3 -1 -1 0 1
2 0 -3 2 -2 -1 -1 -1 -1 -2 1 0 0 -1 -3 2 -1 -1 -1 1 -3 1 0 -1
-5 0 1 -6 6 -2 -3 -1 3 -3 2 -1 -3 2 -2 -2 2 -2 -4 1 3 -4 2 1
1 2 -6 3 -2 -3 -3 3 -4 -2 0 -2 3 1 -3 4 -5 0 1 3 -6 1 0 -1
4 -1 2 2 -2 2 5 -1 0 2 0 2 2 -3 3 0 1 2 1 -1 -2 2 0 -1
2 -3 8 0 0 4 5 -3 5 3 1 2 -1 -2 6 -3 5 0 1 -3 6 1 -2 1
5 1 -1 4 -3 2 6 -1 1 2 0 3 3 -3 2 0 2 1 1 -1 -3 1 1 -2
9 -5 8 5 -4 6 8 -4 1 4 3 2 1 -6 8 -1 4 2 2 -3 -1 5 -2 -2
2 4 -11 5 -4 -3 0 3 -6 0 -3 0 5 -1 -4 3 -3 1 2 1 -8 1 3 -3
1 -2 10 -3 0 7 9 -4 7 5 -2 7 1 -3 7 -4 7 4 0 -6 11 -2 2 3
3 0 4 2 -4 4 7 -1 4 3 -2 5 6 -3 6 0 2 4 3 -4 6 0 1 2
2 0 2 2 -4 3 5 0 1 4 -4 5 3 -1 2 2 0 4 2 -2 3 2 -1 2
-4 0 9 -10 6 5 4 -6 12 -1 0 7 -4 0 -1 -6 8 0 -8 -4 13 -9 4 4"""
# The fixed U matrix is checked by exact determinant and congruence identities.

def parse_mat(s): return np.array([[int(x) for x in line.split()] for line in s.strip().splitlines()],dtype=np.int64)
def rank7(A):
 A=np.array(A,dtype=np.int64)%7; n,m=A.shape;r=0
 for col in range(m):
  if r>=n:break
  piv=next((i for i in range(r,n) if A[i,col]!=0),None)
  if piv is None:continue
  if piv!=r:A[[piv,r]]=A[[r,piv]]
  A[r]=(A[r]*pow(int(A[r,col]),-1,7))%7
  for i in range(n):
   if i!=r and A[i,col]:A[i]=(A[i]-A[i,col]*A[r])%7
  r+=1
 return r

def inv7(A):
 n=A.shape[0];X=np.column_stack((np.array(A,dtype=np.int64)%7,np.eye(n,dtype=np.int64)))%7
 for j in range(n):
  piv=next((i for i in range(j,n) if X[i,j]),None)
  assert piv is not None
  X[[j,piv]]=X[[piv,j]]
  X[j]=X[j]*pow(int(X[j,j]),-1,7)%7
  for i in range(n):
   if i!=j and X[i,j]:X[i]=(X[i]-X[i,j]*X[j])%7
 return X[:,n:]

def ceilq(f):return -((-f.numerator)//f.denominator)
def floorq(f):return f.numerator//f.denominator

def exact_root_enumerate(Pr):
 """Complete Fincke-Pohst at bound 2 using only Fraction and integer sqrt.
 Bounds deliberately widened by 1; no floating-point pruning.
 """
 M=sp.Matrix(Pr.tolist());L,D=M.LDLdecomposition(hermitian=False)
 assert L*D*L.T==M
 d=[Fraction(int(D[i,i].p),int(D[i,i].q)) for i in range(24)]
 ell=[[Fraction(int(L[j,i].p),int(L[j,i].q)) for j in range(24)] for i in range(24)]
 assert min(d)>0
 nodes=0;roots=0;z=[0]*24
 def rec(i,rem):
  nonlocal nodes,roots
  nodes+=1
  if i<0:
   if any(z):roots+=1
   return
  shift=sum(ell[i][j]*z[j] for j in range(i+1,24))
  B=rem/d[i]
  if B<0:return
  radius=math.isqrt(ceilq(B))+1
  lo=ceilq(-shift)-radius
  hi=floorq(-shift)+radius
  for k in range(lo,hi+1):
   cost=d[i]*(k+shift)**2
   if cost<=rem:
    z[i]=k;rec(i-1,rem-cost)
  z[i]=0
 rec(23,Fraction(2))
 return nodes,roots,min(d)

def main():
 start=time.time()
 AP,BP=parse_mat(AP_TEXT),parse_mat(BP_TEXT)
 R=np.array([[int(c) for c in row] for row in R_TEXT.splitlines()],dtype=np.int64)
 Ured=parse_mat(U_TEXT)
 assert R.shape==(24,24) and Ured.shape==(24,24)
 Pnum=np.zeros((24,24),dtype=np.int64)
 Pnum[::2,::2]=2*AP+BP;Pnum[::2,1::2]=AP+4*BP
 Pnum[1::2,::2]=AP-3*BP;Pnum[1::2,1::2]=4*AP+2*BP
 assert np.all(Pnum%7==0)
 P=Pnum//7;H=np.kron(np.eye(2,dtype=np.int64),P)
 assert np.array_equal(P,P.T) and np.all(np.diag(P)%2==0)
 Ps=sp.Matrix(P.tolist()); assert Ps.det(method='domain-ge')==1
 assert sp.det(sp.Matrix(Ured.tolist())) in (-1,1)
 Pr=Ured@P@Ured.T
 assert np.all(np.diag(Pr)==4)
 nodes,roots,dmin=exact_root_enumerate(Pr)
 assert roots==0
 print(f'PASS: P has determinant 1, a norm-4 vector and NO roots (exact FP: {nodes} nodes, {roots} roots, D_min={dmin}).',flush=True)

 AE=np.array([[7,0,-3,-2],[0,7,-5,-3],[3,5,35,21],[2,3,21,14]],dtype=np.int64)
 BE=np.array([[0,0,6,4],[0,0,10,6],[-6,-10,0,0],[-4,-6,0,0]],dtype=np.int64)
 Uh=np.kron(AE,AP)-2*np.kron(BE,BP)
 Vh=np.kron(AE,BP)+np.kron(BE,AP)+np.kron(BE,BP)
 Nh=np.zeros((96,96),dtype=np.int64)
 Nh[::2,::2]=2*Uh+Vh;Nh[::2,1::2]=Uh+4*Vh
 Nh[1::2,::2]=Uh-3*Vh;Nh[1::2,1::2]=4*Uh+2*Vh
 C=np.column_stack((np.eye(24,dtype=np.int64),R))
 AA=np.zeros((48,48),dtype=np.int64);BB=np.zeros((48,48),dtype=np.int64);AA[:24]=C
 for j in range(24):AA[24+j,24+j]=-1;BB[24+j,24+j]=2
 T=np.zeros((96,96),dtype=np.int64)
 for i in range(48):
  for j in range(48):
   aa,bb=int(AA[i,j]),int(BB[i,j])
   T[2*i:2*i+2,2*j:2*j+2]=[[aa,bb],[-2*bb,aa+bb]]
 NG=T@Nh@T.T
 assert np.all(NG%49==0)
 G=NG//49
 assert np.array_equal(G,G.T) and np.all(np.diag(G)%2==0)
 Gs=sp.Matrix(G.tolist());assert Gs.det(method='domain-ge')==1
 A,B,C=Gs[:48,:48],Gs[:48,48:],Gs[48:,48:]
 T0=sp.Matrix([[5,3],[3,2]]);K0=sp.Matrix([[5,-3],[-3,2]])
 U=sp.Matrix([[1,-1],[1,-2]]);W=sp.Matrix([[1,1],[1,2]])
 assert C==7*sp.kronecker_product(T0,Ps)
 assert U*T0*U.T==sp.eye(2) and W*K0*W.T==sp.eye(2)
 D=sp.kronecker_product(T0.inv(),Ps.inv())*B.T
 assert all(e.q==1 for e in D)
 assert 7*A-B*D==sp.kronecker_product(K0,Ps)
 X=sp.kronecker_product(W.T,sp.eye(24));Y=sp.kronecker_product(U.T,sp.eye(24))
 Fsp=Y.inv()*D*X
 assert all(e.q==1 for e in Fsp)
 F=np.array(Fsp.tolist(),dtype=np.int64)%7
 assert rank7(F)==48
 assert np.all((F.T@(H%7)@F+H)%7==0)
 NN=(F@F+np.eye(48,dtype=np.int64))%7
 assert np.all((NN@NN)%7==0) and rank7(NN)==24
 out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "external", "four_leech_G96.txt")
 with open(out, "w") as f:
  f.write("[" + "; ".join(", ".join(str(int(v)) for v in r) for r in G) + "]\n")
 print("written", out)
if __name__ == "__main__":
    main()
