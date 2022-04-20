from Grassmannian import*
from SpinorBrackets3pt import*
from itertools import product


'---------------------- Examples ----------------------'

'''

Generators:
Qa1[A] , Qb1[A] , Qa2[A] , Qb2[A]

building blocks (Lorentz covariant):
etas(i,A)
delta[A] , delta12
Qa[A] , Qb[A] , R[A]

Super amplitude:
solution[#]

'''



'---------------------- Generators ----------------------'


# first set

QA1={}
QB1={}

QA1[1]=Operator_zeros()
for i in range(1,3):
    for I in range(2):
        QA1[1]+=Num(ref*RA[i][:,I])*uSUSY(str(i)+str(I)+'1')
QA1[1]+=SUSY(Num(ref*RA[0]),['01'])

QA1[2]=Operator_zeros()
for i in range(1,3):
    for I in range(2):
        for J in range(2):
            QA1[2]+=(-1)*sig[i]*Num(ref*RA[i][:,I])*epsilon_L[I,J]*uDiff(str(i)+str(J)+'1')
QA1[2]+=SUSY(Num(ref*RA[0]),['02'])

QB1[1]=Operator_zeros()
for i in range(1,3):
    for I in range(2):
        for J in range(2):
            QB1[1]+=(-1)*Num(ref*RB[i][:,I])*epsilon_L[I,J]*uDiff(str(i)+str(J)+'1')
QB1[1]+=Diff(Num(ref*RB[0]),'01')

QB1[2]=SUSY_zeros()
for i in range(1,3):
    for I in range(2):
        QB1[2]+=(-1)*sig[i]*Num(ref*RB[i][:,I])*uSUSY(str(i)+str(I)+'1')
QB1[2]+=Diff(Num(ref*RB[0]),'02')


# second set

QA2={}
QB2={}

QA2[1]=Operator_zeros()
for i in range(1,3):
    for I in range(2):
        for J in range(2):
            QA2[1]+=sig[i]*Num(ref*RA[i][:,I])*epsilon_L[I,J]*uDiff(str(i)+str(J)+'2')
QA2[1]+=SUSY(Num(ref*RA[0]),['01'])

QA2[2]=Operator_zeros()
for i in range(1,3):
    for I in range(2):
        QA2[2]+=Num(ref*RA[i][:,I])*uSUSY(str(i)+str(I)+'2')
QA2[2]+=SUSY(Num(ref*RA[0]),['02'])

QB2[1]=Operator_zeros()
for i in range(1,3):
    for I in range(2):
            QB2[1]+=sig[i]*Num(ref*RB[i][:,I])*uSUSY(str(i)+str(I)+'2')
QB2[1]+=Diff(Num(ref*RB[0]),'01')

QB2[2]=SUSY_zeros()
for i in range(1,3):
    for I in range(2):
        for J in range(2):
            QB2[2]+=(-1)*Num(ref*RB[i][:,I])*epsilon_L[I,J]*uDiff(str(i)+str(J)+'2')
QB2[2]+=Diff(Num(ref*RB[0]),'02')







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


etaPs=0
for (A,B) in product(*([range(1,3)]*2)):
    etaPs+=(1/2)*epsilonL(A,B)*etaP[A]*etaP[B]



zeta=np.matrix(np.random.randn(1,2)+np.ones((1,2)))
xi=np.matrix(np.random.randn(1,2)+np.ones((1,2)))



Ha={}
for A in range(1,3):
    Ha[A]=SUSY_zeros()
    for i in range(1,3):
        for I in range(2):
            Ha[A]+=SUSY(Num(zeta*RA[i][:,I])/Num(zeta*RA[0]),[str(i)+str(I)+str(A)])


Hb={}
for A in range(1,3):
    Hb[A]=SUSY_zeros()
    for B in range(1,3):
        for i in range(1,3):
            for I in range(2):
                Hb[A]+=SUSY(Num(xi*RB[i][:,I])*sig[i]*epsilonL(A,B)/Num(xi*RB[0]),[str(i)+str(I)+str(B)])


HHa=-Ha[1]*Ha[2]
HHb=Hb[1]*Hb[2]
HHab=(1/2)*(Ha[1]*Hb[1]+Ha[2]*Hb[2])




'---------------------- Super amplitude ----------------------'



solution={}
solution[1]=etaP[1]*delta[1]-delta[1]*Ha[1]-delta[1]*Hb[2]*etaP[1]*etaP[2]+delta[1]*Ha[1]*Hb[2]*etaP[2]
solution[2]=etaP[2]*delta[2]-delta[2]*Ha[2]+delta[2]*Hb[1]*etaP[1]*etaP[2]+delta[2]*Ha[2]*Hb[1]*etaP[1]




'---------------------- verify ----------------------'


if __name__=='__main__':

    print ((
    len(QA1[1].operate(solution[1]).estimate().arr)
    +len(QA1[2].operate(solution[1]).estimate().arr)
    +len(QB1[1].operate(solution[1]).estimate().arr)
    +len(QB1[2].operate(solution[1]).estimate().arr)
    +len(QA2[1].operate(solution[2]).estimate().arr)
    +len(QA2[2].operate(solution[2]).estimate().arr)
    +len(QB2[1].operate(solution[2]).estimate().arr)
    +len(QB2[2].operate(solution[2]).estimate().arr))
    ==0)




