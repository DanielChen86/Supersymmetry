import numpy as np
import sympy as sym
import pandas as pd
from Grassmannian import*
from numpy import sqrt

'---------------------- Examples ----------------------'

'''

Eta([]) , Eta(['111']) , Eta(['101','202','01'])
Eta.etalist , Eta.sign

SUSY(1.,[]) , SUSY(3.+2.j,['111','112'])
SUSY.number , SUSY.etalist , SUSY.order
SUSY.copy() , SUSY.vanish(1e-6)

addSUSY([SUSY,SUSY,SUSY,...])
addSUSY.arr , addSUSY.alleta
addSUSY.order(3) , addSUSY.get(['101','01']) , addSUSY.estimate(1e-6) , addSUSY.copy() , addSUSY.vanish(1e-6)
addSUSY_zeros()

diff(3.,'111') , diff(SUSY,'111')

Diff(3.,'111')
Diff.number , Diff.differentials
Diff.diff(SUSY) , Diff.diff(addSUSY) , Diff.copy()

addDiff([Diff,Diff,Diff,...])
addDiff.arr , addDiff.alldiff
addDiff.diff(SUSY) , addDiff.diff(addSUSY) , addDiff.copy()

Operator(multi,diff)    # multi=SUSY/addSUSY , diff=Diff/addDiff
Operator.operate(SUSY) , Operator.operate(addSUSY)

'''



'---------------------- basic ingredients ----------------------'


data=pd.read_csv('spinor helicity data (4pt)/spinor helicity data 1.csv')
data=data['data']


def get(matrix,precision=5):
    return sym.Matrix(matrix.round(precision))


unit=np.matrix([[1.,0.],[0.,1.]])

epsilon_U=np.matrix([[0.,1.],[-1.,0.]])
epsilon_L=np.matrix([[0.,-1.],[1.,0.]])


def epsilonL(num1,num2,err=True):
    if num1==num2:
        return 0.
    elif num2%2==0 and num1==num2-1:
        return complex(epsilon_L[0,1])
    elif num1%2==0 and num1==num2+1:
        return complex(epsilon_L[1,0])
    else:
        if err:
            print ("error: epsilonL")
        return 0.

def epsilonU(num1,num2,err=True):
    if num1==num2:
        return 0.
    elif num2%2==0 and num1==num2-1:
        return complex(epsilon_U[0,1])
    elif num1%2==0 and num1==num2+1:
        return complex(epsilon_U[1,0])
    else:
        if err:
            print ("error: epsilonU")
        return 0.



sigma={}
sigma[0]=np.matrix([[1.,0.],[0.,1.]])
sigma[1]=np.matrix([[0.,1.],[1.,0.]])
sigma[2]=np.matrix([[0.,-1.j],[1.j,0.]])
sigma[3]=np.matrix([[1.,0.],[0.,-1.]])

metric=np.matrix(np.diag([1.,-1.,-1.,-1.]))


'---------------------- spinor helicity brackets ----------------------'


RA={}
RA[0]=np.transpose(np.matrix([data[0]+data[1]*1.j,data[2]+data[3]*1.j]))
RA[1]=np.matrix([[data[4]+data[5]*1.j,data[6]+data[7]*1.j],[data[8]+data[9]*1.j,data[10]+data[11]*1.j]])
RA[2]=np.matrix([[data[12]+data[13]*1.j,data[14]+data[15]*1.j],[data[16]+data[17]*1.j,data[18]+data[19]*1.j]])
RA[3]=np.matrix([[data[20]+data[21]*1.j,data[22]+data[23]*1.j],[data[24]+data[25]*1.j,data[26]+data[27]*1.j]])
RA[4]=np.matrix([[data[28]+data[29]*1.j,data[30]+data[31]*1.j],[data[32]+data[33]*1.j,data[34]+data[35]*1.j]])



LB={}
LB[0]=np.matrix([data[36]+data[37]*1.j,data[38]+data[39]*1.j])
LB[1]=np.matrix([[data[40]+data[41]*1.j,data[42]+data[43]*1.j],[data[44]+data[45]*1.j,data[46]+data[47]*1.j]])
LB[2]=np.matrix([[data[48]+data[49]*1.j,data[50]+data[51]*1.j],[data[52]+data[53]*1.j,data[54]+data[55]*1.j]])
LB[3]=np.matrix([[data[56]+data[57]*1.j,data[58]+data[59]*1.j],[data[60]+data[61]*1.j,data[62]+data[63]*1.j]])
LB[4]=np.matrix([[data[64]+data[65]*1.j,data[66]+data[67]*1.j],[data[68]+data[69]*1.j,data[70]+data[71]*1.j]])


m=data[72]


for i in range(5):
    RA[i]=RA[i]/sqrt(m)
    LB[i]=LB[i]/sqrt(m)

m=1


LA={}
RB={}
LA[0]=np.transpose(epsilon_U*RA[0])
RB[0]=epsilon_U*np.transpose(LB[0])
for i in range(1,5):
    LA[i]=np.transpose(epsilon_U*RA[i])
    RB[i]=epsilon_U*np.transpose(LB[i])


pL={}
pU={}
pL[0]=RA[0]*LB[0]
pU[0]=RB[0]*LA[0]
for i in range(1,5):
    pL[i]=RA[i]*epsilon_L*LB[i]
    pU[i]=-RB[i]*epsilon_L*LA[i]


Pvec={}
for i in range(5):
    Pvec[i]=np.transpose(np.matrix((0.+0.j,0.+0.j,0.+0.j,0.+0.j)))
    for j in range(4):
        Pvec[i][j,0]=np.trace(pU[i]*sigma[j]/2)


xL=(pL[1]*RB[0])[0,0]/m/(RA[0][0,0])
xR=-(pL[3]*RB[0])[0,0]/m/(RA[0][0,0])


'---------------------- useful tools ----------------------'

def Num(matrix):
    if matrix.shape==(1,1):
        return complex(matrix[0,0])
    else:
        print ("error")
        return None



sig={1:1.,2:-1.,3:1.,4:-1.}


ref=np.matrix(np.random.randn(1,2)+np.ones((1,2)))



'---------------------- verify ----------------------'

if __name__=='__main__':
    print (get(LA[1]*RA[1])==-get(LB[1]*RB[1]))
    print (abs(RB[0][0,0]-((pU[1]*RA[0])[0,0]*xL/m))<1e-8)
    print (abs(RB[0][0,0]+((pU[3]*RA[0])[0,0]*xR/m))<1e-8)


