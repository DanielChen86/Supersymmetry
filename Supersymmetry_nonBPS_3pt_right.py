import numpy as np
import sympy as sym
import pandas as pd
from Grassmannian import*



'---------------------- Examples ----------------------'

'''

Generators:
FZ , GZ
Da_upper[A] , Db_upper[A] , Da_lower[A] , Db_lower[A]
Da[A] , Db[A] , QA[A] , QB[A]

building blocks:
delta[A] , delta12
etaP[A] , etas(i,A)
Qa[A] , Qb[A] , R[A]

Super amplitude:
solution[#]

'''






'---------------------- basic ingredients ----------------------'


data=pd.read_csv('spinor helicity data (3pt) (right)/spinor helicity data 1 (right).csv')
data=data['data']


def get(matrix,precision=5):
    return sym.Matrix(matrix.round(precision))


unit=np.matrix([[1.,0.],[0.,1.]])

epsilon_U=np.matrix([[0.,1.],[-1.,0.]])
epsilon_L=np.matrix([[0.,-1.],[1.,0.]])



def epsilonL(num1,num2):
    if (num1 in (1,2)) and (num2 in (1,2)):
        return complex(epsilon_L[num1-1,num2-1])
    else:
        print ("error")
        return 0.


def epsilonU(num1,num2):
    if (num1 in (1,2)) and (num2 in (1,2)):
        return complex(epsilon_U[num1-1,num2-1])
    else:
        print ("error")
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


LB={}
LB[0]=np.matrix([data[20]+data[21]*1.j,data[22]+data[23]*1.j])
LB[1]=np.matrix([[data[24]+data[25]*1.j,data[26]+data[27]*1.j],[data[28]+data[29]*1.j,data[30]+data[31]*1.j]])
LB[2]=np.matrix([[data[32]+data[33]*1.j,data[34]+data[35]*1.j],[data[36]+data[37]*1.j,data[38]+data[39]*1.j]])


m=data[40]


for i in range(3):
    RA[i]=RA[i]/np.sqrt(m)
    LB[i]=LB[i]/np.sqrt(m)

m=1


LA={}
RB={}
LA[0]=np.transpose(epsilon_U*RA[0])
RB[0]=epsilon_U*np.transpose(LB[0])
for i in range(1,3):
    LA[i]=np.transpose(epsilon_U*RA[i])
    RB[i]=epsilon_U*np.transpose(LB[i])


pL={}
pU={}
pL[0]=RA[0]*LB[0]
pU[0]=RB[0]*LA[0]
for i in range(1,3):
    pL[i]=RA[i]*epsilon_L*LB[i]
    pU[i]=-RB[i]*epsilon_L*LA[i]


Pvec={}
for i in range(3):
    Pvec[i]=np.transpose(np.matrix((0.+0.j,0.+0.j,0.+0.j,0.+0.j)))
    for j in range(4):
        Pvec[i][j,0]=np.trace(pU[i]*sigma[j]/2)


x12=-(pL[1]*RB[0])[0,0]/m/(RA[0][0,0])


'---------------------- useful tools ----------------------'

def Num(matrix):
    if matrix.shape==(1,1):
        return complex(matrix[0,0])
    else:
        print ("error")
        return None



sig={1:1.,2:-1.}


ref=np.matrix(np.random.randn(1,2)+np.ones((1,2)))






'---------------------- Generators ----------------------'


FZ=1/np.sqrt(2)
GZ=float(np.sqrt(1-FZ**2))



Da_upper={}
for A in range(1,3):
    Da_upper[A]=SUSY_zeros()
    for i in range(1,3):
        for I in range(2):
            Da_upper[A]+=SUSY(FZ*Num(ref*RA[i][:,I]),[str(i)+str(I)+str(A)])

Db_upper={}
for A in range(1,3):
    Db_upper[A]=SUSY_zeros()
    for B in range(1,3):
        for i in range(1,3):
            for I in range(2):
                Db_upper[A]+=SUSY((-GZ)*Num(ref*RB[i][:,I])*sig[i]*epsilonL(A,B),[str(i)+str(I)+str(B)])

Da_lower={}
for A in range(1,3):
    Da_lower[A]=Diff_zeros()
    for i in range(1,3):
        for I in range(2):
            for J in range(2):
                for B in range(1,3):
                    Da_lower[A]+=Diff(GZ*sig[i]*Num(ref*RA[i][:,I])*epsilonU(A,B)*epsilon_L[I,J],str(i)+str(J)+str(B))

Db_lower={}
for A in range(1,3):
    Db_lower[A]=Diff_zeros()
    for i in range(1,3):
        for I in range(2):
            for J in range(2):
                Db_lower[A]+=Diff(-FZ*Num(ref*RB[i][:,I])*epsilon_L[I,J],str(i)+str(J)+str(A))


Da={}
Db={}
QA={}
QB={}
for A in range(1,3):
    Da[A]=Da_upper[A]+Da_lower[A]
    Db[A]=Db_upper[A]+Db_lower[A]
    QA[A]=Da[A]-Num(ref*RA[0])*uSUSY('0'+str(A))
    QB[A]=Db[A]+Num(ref*RB[0])*uDiff('0'+str(A))



'---------------------- building blocks ----------------------'


delta={}
for A in range(1,3):
    delta[A]=SUSY_zeros()
    for i in range(1,3):
        for I in range(2):
            delta[A]+=SUSY(Num(LA[0]*RA[i][:,I]),[str(i)+str(I)+str(A)])
delta12=delta[1]*delta[2]



etaP={}
etaP[1]=uSUSY('01')
etaP[2]=uSUSY('02')


def etas(particle,Rcharge):
    result=SUSY_zeros()
    for I in range(2):
        for J in range(2):
            result+=(-0.5)*epsilon_U[I,J]*SUSY(1.,[str(particle)+str(I)+str(Rcharge),str(particle)+str(J)+str(Rcharge)])
    return result


Qa={}
for A in range(1,3):
    Qa[A]=SUSY_zeros()
    for I in range(2):
        for J in range(2):
            Qa[A]+=(LA[1]*RA[2])[I,J]*uSUSY(str(1)+str(I)+str(A))*uSUSY(str(2)+str(J)+str(A))
    Qa[A]+=(m*etas(1,A)+m*etas(2,A))



Qb={}
for A in range(1,3):
    Qb[A]=SUSY_zeros()
    for I in range(2):
        for J in range(2):
            Qb[A]+=(LB[1]*RB[2])[I,J]*uSUSY(str(1)+str(I)+str(A))*uSUSY(str(2)+str(J)+str(A))
    Qb[A]+=(m*etas(1,A)+m*etas(2,A))



R={}
for A in range(1,3):
    R[A]=SUSY_zeros()
    for I in range(2):
        R[A]+=SUSY(Num(LB[0]*RB[1][:,I]),[str(1)+str(I)+str(A)])*etas(2,A)+SUSY(Num(LB[0]*RB[2][:,I]),[str(2)+str(I)+str(A)])*etas(1,A)






'---------------------- Super amplitude ----------------------'



solution={}
solution[1]=\
etaP[1]*delta[1]*etaP[2]*delta[2]\
+(GZ/FZ/x12)*(Qa[1]*etaP[1]*Qb[2]*etaP[2]+Qb[1]*etaP[1]*Qa[2]*etaP[2])\
+(GZ**2/FZ**2/x12**2)*R[1]*etaP[1]*R[2]*etaP[2]\
-(1/FZ)*(etaP[1]*delta[1]*Qa[2]+Qa[1]*etaP[2]*delta[2])\
-(GZ/FZ**2/x12)*(Qa[1]*etaP[1]*R[2]+R[1]*Qa[2]*etaP[2])\
+(1/FZ**2)*Qa[1]*Qa[2]

solution[2]=\
(-1/x12**2)*Qb[1]*etaP[1]*Qb[2]*etaP[2]\
+(FZ/x12**2)*(Qb[1]*etaP[1]*R[2]+R[1]*Qb[2]*etaP[2])\
+(GZ/x12)*(etaP[1]*delta[1]*Qb[2]+Qb[1]*etaP[2]*delta[2])\
-(GZ**2)*delta12\
-(FZ*GZ/x12)*(Qa[1]*Qb[2]+Qb[1]*Qa[2])\
-(FZ**2/x12**2)*R[1]*R[2]





'---------------------- verify ----------------------'

if __name__=='__main__':

    print((
    len(QA[1].operate(solution[1]).estimate().arr)\
    +len(QA[2].operate(solution[1]).estimate().arr)\
    +len(QB[1].operate(solution[1]).estimate().arr)\
    +len(QB[2].operate(solution[1]).estimate().arr)\
    +len(QA[1].operate(solution[2]).estimate().arr)\
    +len(QA[2].operate(solution[2]).estimate().arr)\
    +len(QB[1].operate(solution[2]).estimate().arr)\
    +len(QB[2].operate(solution[2]).estimate().arr))\
    ==0)




