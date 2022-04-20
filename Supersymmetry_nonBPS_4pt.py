import numpy as np
from Grassmannian import*
from SpinorBrackets4pt import*
from itertools import product
from numpy import sqrt

'---------------------- Central charges ----------------------'


FL12=0.2
FL34=0.1
FR12=0.1
FR34=0.1
GL12=float(sqrt(1-FL12**2))
GL34=float(sqrt(1-FL34**2))
GR12=float(sqrt(1-FR12**2))
GR34=float(sqrt(1-FR34**2))


'---------------------- Generators (left) ----------------------'


DaL_upper={}
for A in range(1,3):
    DaL_upper[A]=SUSY_zeros()
    for i in range(1,3):
        for I in range(2):
            DaL_upper[A]+=SUSY(FL12*Num(ref*RA[i][:,I]),[str(i)+str(I)+str(A)])
for A in range(3,5):
    DaL_upper[A]=SUSY_zeros()
    for i in range(1,3):
        for I in range(2):
            DaL_upper[A]+=SUSY(FL34*Num(ref*RA[i][:,I]),[str(i)+str(I)+str(A)])

DbL_upper={}
for A in range(1,3):
    DbL_upper[A]=SUSY_zeros()
    for B in range(1,3):
        for i in range(1,3):
            for I in range(2):
                DbL_upper[A]+=SUSY((-GL12)*Num(ref*RB[i][:,I])*sig[i]*epsilonL(A,B),[str(i)+str(I)+str(B)])
for A in range(3,5):
    DbL_upper[A]=SUSY_zeros()
    for B in range(3,5):
        for i in range(1,3):
            for I in range(2):
                DbL_upper[A]+=SUSY((-GL34)*Num(ref*RB[i][:,I])*sig[i]*epsilonL(A,B),[str(i)+str(I)+str(B)])


DaL_lower={}
for A in range(1,3):
    DaL_lower[A]=Diff_zeros()
    for i in range(1,3):
        for I in range(2):
            for J in range(2):
                for B in range(1,3):
                    DaL_lower[A]+=Diff(GL12*sig[i]*Num(ref*RA[i][:,I])*epsilonU(A,B)*epsilon_L[I,J],str(i)+str(J)+str(B))
for A in range(3,5):
    DaL_lower[A]=Diff_zeros()
    for i in range(1,3):
        for I in range(2):
            for J in range(2):
                for B in range(3,5):
                    DaL_lower[A]+=Diff(GL34*sig[i]*Num(ref*RA[i][:,I])*epsilonU(A,B)*epsilon_L[I,J],str(i)+str(J)+str(B))


DbL_lower={}
for A in range(1,3):
    DbL_lower[A]=Diff_zeros()
    for i in range(1,3):
        for I in range(2):
            for J in range(2):
                DbL_lower[A]+=Diff(-FL12*Num(ref*RB[i][:,I])*epsilon_L[I,J],str(i)+str(J)+str(A))
for A in range(3,5):
    DbL_lower[A]=Diff_zeros()
    for i in range(1,3):
        for I in range(2):
            for J in range(2):
                DbL_lower[A]+=Diff(-FL34*Num(ref*RB[i][:,I])*epsilon_L[I,J],str(i)+str(J)+str(A))


DaL={}
DbL={}
QAL={}
QBL={}
for A in range(1,5):
    DaL[A]=DaL_upper[A]+DaL_lower[A]
    DbL[A]=DbL_upper[A]+DbL_lower[A]
    QAL[A]=DaL[A]+Num(ref*RA[0])*uSUSY('0'+str(A))
    QBL[A]=DbL[A]+Num(ref*RB[0])*uDiff('0'+str(A))




'---------------------- Generators (right) ----------------------'

DaR_upper={}
for A in range(1,3):
    DaR_upper[A]=SUSY_zeros()
    for i in range(3,5):
        for I in range(2):
            DaR_upper[A]+=SUSY(FR12*Num(ref*RA[i][:,I]),[str(i)+str(I)+str(A)])
for A in range(3,5):
    DaR_upper[A]=SUSY_zeros()
    for i in range(3,5):
        for I in range(2):
            DaR_upper[A]+=SUSY(FR34*Num(ref*RA[i][:,I]),[str(i)+str(I)+str(A)])



DbR_upper={}
for A in range(1,3):
    DbR_upper[A]=SUSY_zeros()
    for B in range(1,3):
        for i in range(3,5):
            for I in range(2):
                DbR_upper[A]+=SUSY((-GR12)*Num(ref*RB[i][:,I])*sig[i]*epsilonL(A,B),[str(i)+str(I)+str(B)])
for A in range(3,5):
    DbR_upper[A]=SUSY_zeros()
    for B in range(3,5):
        for i in range(3,5):
            for I in range(2):
                DbR_upper[A]+=SUSY((-GR34)*Num(ref*RB[i][:,I])*sig[i]*epsilonL(A,B),[str(i)+str(I)+str(B)])


DaR_lower={}
for A in range(1,3):
    DaR_lower[A]=Diff_zeros()
    for i in range(3,5):
        for I in range(2):
            for J in range(2):
                for B in range(1,3):
                    DaR_lower[A]+=Diff(GR12*sig[i]*Num(ref*RA[i][:,I])*epsilonU(A,B)*epsilon_L[I,J],str(i)+str(J)+str(B))
for A in range(3,5):
    DaR_lower[A]=Diff_zeros()
    for i in range(3,5):
        for I in range(2):
            for J in range(2):
                for B in range(3,5):
                    DaR_lower[A]+=Diff(GR34*sig[i]*Num(ref*RA[i][:,I])*epsilonU(A,B)*epsilon_L[I,J],str(i)+str(J)+str(B))



DbR_lower={}
for A in range(1,3):
    DbR_lower[A]=Diff_zeros()
    for i in range(3,5):
        for I in range(2):
            for J in range(2):
                DbR_lower[A]+=Diff(-FR12*Num(ref*RB[i][:,I])*epsilon_L[I,J],str(i)+str(J)+str(A))
for A in range(3,5):
    DbR_lower[A]=Diff_zeros()
    for i in range(3,5):
        for I in range(2):
            for J in range(2):
                DbR_lower[A]+=Diff(-FR34*Num(ref*RB[i][:,I])*epsilon_L[I,J],str(i)+str(J)+str(A))



DaR={}
DbR={}
QAR={}
QBR={}
for A in range(1,5):
    DaR[A]=DaR_upper[A]+DaR_lower[A]
    DbR[A]=DbR_upper[A]+DbR_lower[A]
    QAR[A]=DaR[A]-Num(ref*RA[0])*uSUSY('0'+str(A))
    QBR[A]=DbR[A]+Num(ref*RB[0])*uDiff('0'+str(A))







'---------------------- Generators (left+right) ----------------------'


Da_upper={}
for A in range(1,5):
    Da_upper[A]=SUSY_zeros()

for A in range(1,3):
    for i in range(1,3):
        for I in range(2):
            Da_upper[A]+=SUSY(FL12*Num(ref*RA[i][:,I]),[str(i)+str(I)+str(A)])
for A in range(3,5):
    for i in range(1,3):
        for I in range(2):
            Da_upper[A]+=SUSY(FL34*Num(ref*RA[i][:,I]),[str(i)+str(I)+str(A)])
for A in range(1,3):
    for i in range(3,5):
        for I in range(2):
            Da_upper[A]+=SUSY(FR12*Num(ref*RA[i][:,I]),[str(i)+str(I)+str(A)])
for A in range(3,5):
    for i in range(3,5):
        for I in range(2):
            Da_upper[A]+=SUSY(FR34*Num(ref*RA[i][:,I]),[str(i)+str(I)+str(A)])


Da_lower={}
for A in range(1,5):
    Da_lower[A]=Diff_zeros()

for A in range(1,3):
    for i in range(1,3):
        for I in range(2):
            for J in range(2):
                for B in range(1,3):
                    Da_lower[A]+=Diff(GL12*sig[i]*Num(ref*RA[i][:,I])*epsilonU(A,B)*epsilon_L[I,J],str(i)+str(J)+str(B))
for A in range(3,5):
    for i in range(1,3):
        for I in range(2):
            for J in range(2):
                for B in range(3,5):
                    Da_lower[A]+=Diff(GL34*sig[i]*Num(ref*RA[i][:,I])*epsilonU(A,B)*epsilon_L[I,J],str(i)+str(J)+str(B))
for A in range(1,3):
    for i in range(3,5):
        for I in range(2):
            for J in range(2):
                for B in range(1,3):
                    Da_lower[A]+=Diff(GR12*sig[i]*Num(ref*RA[i][:,I])*epsilonU(A,B)*epsilon_L[I,J],str(i)+str(J)+str(B))
for A in range(3,5):
    for i in range(3,5):
        for I in range(2):
            for J in range(2):
                for B in range(3,5):
                    Da_lower[A]+=Diff(GR34*sig[i]*Num(ref*RA[i][:,I])*epsilonU(A,B)*epsilon_L[I,J],str(i)+str(J)+str(B))




Db_upper={}
for A in range(1,5):
    Db_upper[A]=SUSY_zeros()

for A in range(1,3):
    for B in range(1,3):
        for i in range(1,3):
            for I in range(2):
                Db_upper[A]+=SUSY((-GL12)*Num(ref*RB[i][:,I])*sig[i]*epsilonL(A,B),[str(i)+str(I)+str(B)])
for A in range(3,5):
    for B in range(3,5):
        for i in range(1,3):
            for I in range(2):
                Db_upper[A]+=SUSY((-GL34)*Num(ref*RB[i][:,I])*sig[i]*epsilonL(A,B),[str(i)+str(I)+str(B)])
for A in range(1,3):
    for B in range(1,3):
        for i in range(3,5):
            for I in range(2):
                Db_upper[A]+=SUSY((-GR12)*Num(ref*RB[i][:,I])*sig[i]*epsilonL(A,B),[str(i)+str(I)+str(B)])
for A in range(3,5):
    for B in range(3,5):
        for i in range(3,5):
            for I in range(2):
                Db_upper[A]+=SUSY((-GR34)*Num(ref*RB[i][:,I])*sig[i]*epsilonL(A,B),[str(i)+str(I)+str(B)])


Db_lower={}
for A in range(1,5):
    Db_lower[A]=Diff_zeros()

for A in range(1,3):
    for i in range(1,3):
        for I in range(2):
            for J in range(2):
                Db_lower[A]+=Diff(-FL12*Num(ref*RB[i][:,I])*epsilon_L[I,J],str(i)+str(J)+str(A))
for A in range(3,5):
    for i in range(1,3):
        for I in range(2):
            for J in range(2):
                Db_lower[A]+=Diff(-FL34*Num(ref*RB[i][:,I])*epsilon_L[I,J],str(i)+str(J)+str(A))
for A in range(1,3):
    for i in range(3,5):
        for I in range(2):
            for J in range(2):
                Db_lower[A]+=Diff(-FR12*Num(ref*RB[i][:,I])*epsilon_L[I,J],str(i)+str(J)+str(A))
for A in range(3,5):
    for i in range(3,5):
        for I in range(2):
            for J in range(2):
                Db_lower[A]+=Diff(-FR34*Num(ref*RB[i][:,I])*epsilon_L[I,J],str(i)+str(J)+str(A))



QA={}
QB={}
for A in range(1,5):
    QA[A]=Da_upper[A]+Da_lower[A]
    QB[A]=Db_upper[A]+Db_lower[A]




'---------------------- building blocks ----------------------'



deltaL={}
for A in range(1,5):
    deltaL[A]=SUSY_zeros()
    for i in range(1,3):
        for I in range(2):
            deltaL[A]+=SUSY(Num(LA[0]*RA[i][:,I]),[str(i)+str(I)+str(A)])
deltaL12=deltaL[1]*deltaL[2]
deltaL34=deltaL[3]*deltaL[4]

deltaR={}
for A in range(1,5):
    deltaR[A]=SUSY_zeros()
    for i in range(3,5):
        for I in range(2):
            deltaR[A]+=SUSY(Num(LA[0]*RA[i][:,I]),[str(i)+str(I)+str(A)])
deltaR12=deltaR[1]*deltaR[2]
deltaR34=deltaR[3]*deltaR[4]


etaP={}
etaP[1]=uSUSY('01')
etaP[2]=uSUSY('02')
etaP[3]=uSUSY('03')
etaP[4]=uSUSY('04')
etaPs12=0
for (A,B) in product(*([range(1,3)]*2)):
    etaPs12+=(1/2)*epsilonL(A,B)*etaP[A]*etaP[B]
etaPs34=0
for (A,B) in product(*([range(3,5)]*2)):
    etaPs34+=(1/2)*epsilonL(A,B)*etaP[A]*etaP[B]






zeta=np.matrix(np.random.randn(1,2)+np.ones((1,2)))
xi=np.matrix(np.random.randn(1,2)+np.ones((1,2)))



HaL={}
for A in range(1,5):
    HaL[A]=SUSY_zeros()
    for i in range(1,3):
        for I in range(2):
            HaL[A]+=SUSY(Num(zeta*RA[i][:,I])/Num(zeta*RA[0]),[str(i)+str(I)+str(A)])


HaR={}
for A in range(1,5):
    HaR[A]=SUSY_zeros()
    for i in range(3,5):
        for I in range(2):
            HaR[A]+=SUSY(Num(zeta*RA[i][:,I])/Num(zeta*RA[0]),[str(i)+str(I)+str(A)])



HbL={}
for A in range(1,5):
    HbL[A]=SUSY_zeros()
    for B in range(1,5):
        for i in range(1,3):
            for I in range(2):
                HbL[A]+=SUSY(Num(xi*RB[i][:,I])*sig[i]*epsilonL(A,B,err=False)/Num(xi*RB[0]),[str(i)+str(I)+str(B)])


HbR={}
for A in range(1,5):
    HbR[A]=SUSY_zeros()
    for B in range(1,5):
        for i in range(3,5):
            for I in range(2):
                HbR[A]+=SUSY(Num(xi*RB[i][:,I])*sig[i]*epsilonL(A,B,err=False)/Num(xi*RB[0]),[str(i)+str(I)+str(B)])


HHaL12=-HaL[1]*HaL[2]
HHaL34=-HaL[3]*HaL[4]
HHaR12=-HaR[1]*HaR[2]
HHaR34=-HaR[3]*HaR[4]
HHbL12=HbL[1]*HbL[2]
HHbL34=HbL[3]*HbL[4]
HHbR12=HbR[1]*HbR[2]
HHbR34=HbR[3]*HbR[4]
HHabL12=(1/2)*(HaL[1]*HbL[1]+HaL[2]*HbL[2])
HHabL34=(1/2)*(HaL[3]*HbL[3]+HaL[4]*HbL[4])
HHabR12=(1/2)*(HaR[1]*HbR[1]+HaR[2]*HbR[2])
HHabR34=(1/2)*(HaR[3]*HbR[3]+HaR[4]*HbR[4])



'---------------------- Super amplitude ----------------------'



solutionL12={}

solutionL12[1]=0

solutionL12[1]+=deltaL12*etaPs12*(1+(-2*GL12/FL12)*HHabL12+(GL12**2/FL12**2)*HHaL12*HHbL12)
for (A,B) in product(*([range(1,3)]*2)):
    solutionL12[1]+=deltaL12*etaP[A]*(1/FL12)*epsilonL(A,B)*HaL[B]
for A in range(1,3):
    solutionL12[1]+=deltaL12*etaP[A]*(GL12/(FL12**2))*HHaL12*HbL[A]
solutionL12[1]+=deltaL12*(1/(FL12**2))*HHaL12

solutionL12[1]*=FL12**2


solutionL12[2]=0

solutionL12[2]+=deltaL12*etaPs12*HHbL12
for (A,B) in product(*([range(1,3)]*2)):
    solutionL12[2]+=deltaL12*etaP[A]*FL12*epsilonL(A,B)*HaL[B]*HHbL12
for A in range(1,3):
    solutionL12[2]+=deltaL12*etaP[A]*(GL12)*HbL[A]
solutionL12[2]+=deltaL12*((GL12**2)+(2*FL12*GL12)*HHabL12+(FL12**2)*HHaL12*HHbL12)

solutionL12[2]*=-xL







solutionL34={}

solutionL34[1]=0

solutionL34[1]+=deltaL34*etaPs34*(1+(-2*GL34/FL34)*HHabL34+(GL34**2/FL34**2)*HHaL34*HHbL34)
for (A,B) in product(*([range(3,5)]*2)):
    solutionL34[1]+=deltaL34*etaP[A]*(1/FL34)*epsilonL(A,B)*HaL[B]
for A in range(3,5):
    solutionL34[1]+=deltaL34*etaP[A]*(GL34/(FL34**2))*HHaL34*HbL[A]
solutionL34[1]+=deltaL34*(1/(FL34**2))*HHaL34

solutionL34[1]*=FL34**2


solutionL34[2]=0

solutionL34[2]+=deltaL34*etaPs34*HHbL34
for (A,B) in product(*([range(3,5)]*2)):
    solutionL34[2]+=deltaL34*etaP[A]*FL34*epsilonL(A,B)*HaL[B]*HHbL34
for A in range(3,5):
    solutionL34[2]+=deltaL34*etaP[A]*(GL34)*HbL[A]
solutionL34[2]+=deltaL34*((GL34**2)+(2*FL34*GL34)*HHabL34+(FL34**2)*HHaL34*HHbL34)

solutionL34[2]*=-xL






solutionR12={}

solutionR12[1]=0

solutionR12[1]+=deltaR12*etaPs12*(1+(2*GR12/FR12)*HHabR12+(GR12**2/FR12**2)*HHaR12*HHbR12)
for (A,B) in product(*([range(1,3)]*2)):
    solutionR12[1]-=deltaR12*etaP[A]*(1/FR12)*epsilonL(A,B)*HaR[B]
for A in range(1,3):
    solutionR12[1]+=deltaR12*etaP[A]*(GR12/(FR12**2))*HHaR12*HbR[A]
solutionR12[1]+=deltaR12*(1/(FR12**2))*HHaR12

solutionR12[1]*=FR12**2


solutionR12[2]=0

solutionR12[2]+=deltaR12*etaPs12*HHbR12
for (A,B) in product(*([range(1,3)]*2)):
    solutionR12[2]-=deltaR12*etaP[A]*FR12*epsilonL(A,B)*HaR[B]*HHbR12
for A in range(1,3):
    solutionR12[2]+=deltaR12*etaP[A]*(GR12)*HbR[A]
solutionR12[2]+=deltaR12*((GR12**2)-(2*FR12*GR12)*HHabR12+(FR12**2)*HHaR12*HHbR12)

solutionR12[2]*=-xR






solutionR34={}

solutionR34[1]=0

solutionR34[1]+=deltaR34*etaPs34*(1+(2*GR34/FR34)*HHabR34+(GR34**2/FR34**2)*HHaR34*HHbR34)
for (A,B) in product(*([range(3,5)]*2)):
    solutionR34[1]-=deltaR34*etaP[A]*(1/FR34)*epsilonL(A,B)*HaR[B]
for A in range(3,5):
    solutionR34[1]+=deltaR34*etaP[A]*(GR34/(FR34**2))*HHaR34*HbR[A]
solutionR34[1]+=deltaR34*(1/(FR34**2))*HHaR34

solutionR34[1]*=FR34**2


solutionR34[2]=0

solutionR34[2]+=deltaR34*etaPs34*HHbR34
for (A,B) in product(*([range(3,5)]*2)):
    solutionR34[2]-=deltaR34*etaP[A]*FR34*epsilonL(A,B)*HaR[B]*HHbR34
for A in range(3,5):
    solutionR34[2]+=deltaR34*etaP[A]*(GR34)*HbR[A]
solutionR34[2]+=deltaR34*((GR34**2)-(2*FR34*GR34)*HHabR34+(FR34**2)*HHaR34*HHbR34)

solutionR34[2]*=-xR








'---------------------- verify ----------------------'


if __name__=='__main__':

    accuarcy=1e-8

    print((
    len(QAL[1].operate(solutionL12[1]).estimate(accuarcy).arr)
    +len(QAL[1].operate(solutionL12[2]).estimate(accuarcy).arr)
    +len(QAL[2].operate(solutionL12[1]).estimate(accuarcy).arr)
    +len(QAL[2].operate(solutionL12[2]).estimate(accuarcy).arr)
    +len(QAL[3].operate(solutionL34[1]).estimate(accuarcy).arr)
    +len(QAL[3].operate(solutionL34[2]).estimate(accuarcy).arr)
    +len(QAL[4].operate(solutionL34[1]).estimate(accuarcy).arr)
    +len(QAL[4].operate(solutionL34[2]).estimate(accuarcy).arr)
    +len(QAR[1].operate(solutionR12[1]).estimate(accuarcy).arr)
    +len(QAR[1].operate(solutionR12[2]).estimate(accuarcy).arr)
    +len(QAR[2].operate(solutionR12[1]).estimate(accuarcy).arr)
    +len(QAR[2].operate(solutionR12[2]).estimate(accuarcy).arr)
    +len(QAR[3].operate(solutionR34[1]).estimate(accuarcy).arr)
    +len(QAR[3].operate(solutionR34[2]).estimate(accuarcy).arr)
    +len(QAR[4].operate(solutionR34[1]).estimate(accuarcy).arr)
    +len(QAR[4].operate(solutionR34[2]).estimate(accuarcy).arr)
    +len(QBL[1].operate(solutionL12[1]).estimate(accuarcy).arr)
    +len(QBL[1].operate(solutionL12[2]).estimate(accuarcy).arr)
    +len(QBL[2].operate(solutionL12[1]).estimate(accuarcy).arr)
    +len(QBL[2].operate(solutionL12[2]).estimate(accuarcy).arr)
    +len(QBL[3].operate(solutionL34[1]).estimate(accuarcy).arr)
    +len(QBL[3].operate(solutionL34[2]).estimate(accuarcy).arr)
    +len(QBL[4].operate(solutionL34[1]).estimate(accuarcy).arr)
    +len(QBL[4].operate(solutionL34[2]).estimate(accuarcy).arr)
    +len(QBR[1].operate(solutionR12[1]).estimate(accuarcy).arr)
    +len(QBR[1].operate(solutionR12[2]).estimate(accuarcy).arr)
    +len(QBR[2].operate(solutionR12[1]).estimate(accuarcy).arr)
    +len(QBR[2].operate(solutionR12[2]).estimate(accuarcy).arr)
    +len(QBR[3].operate(solutionR34[1]).estimate(accuarcy).arr)
    +len(QBR[3].operate(solutionR34[2]).estimate(accuarcy).arr)
    +len(QBR[4].operate(solutionR34[1]).estimate(accuarcy).arr)
    +len(QBR[4].operate(solutionR34[2]).estimate(accuarcy).arr))
    ==0)



