from Grassmannian import*
from SpinorBrackets3pt import*



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





'---------------------- building blocks (Lorentz covariant) ----------------------'


def etas(particle,Rcharge):
    result=SUSY_zeros()
    for I in range(2):
        for J in range(2):
            result+=(-0.5)*epsilon_U[I,J]*SUSY(1.,[str(particle)+str(I)+str(Rcharge),str(particle)+str(J)+str(Rcharge)])
    return result


delta={}
for A in range(1,3):
    delta[A]=SUSY_zeros()
    for i in range(1,3):
        for I in range(2):
            delta[A]+=SUSY(Num(LA[0]*RA[i][:,I]),[str(i)+str(I)+str(A)])
delta12=delta[1]*delta[2]



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
solution[1]=Qa[1]+(1/x12)*Qb[1]*uSUSY('01')*uSUSY('02')+uSUSY('01')*delta[1]+(1/x12)*uSUSY('02')*R[1]
solution[2]=Qa[2]+(1/x12)*Qb[2]*uSUSY('01')*uSUSY('02')+uSUSY('02')*delta[2]-(1/x12)*uSUSY('01')*R[2]



'---------------------- verify ----------------------'



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



