import numpy as np
from Grassmannian2 import*
from SpinorBrackets3pt import*
from itertools import product
from numpy import sqrt

'---------------------- Central charges ----------------------'


FZ12=0.2
FZ34=0.1
FZ56=0.3
FZ78=0.4
GZ12=float(sqrt(1-FZ12**2))
GZ34=float(sqrt(1-FZ34**2))
GZ56=float(sqrt(1-FZ56**2))
GZ78=float(sqrt(1-FZ78**2))


'---------------------- Generators (left) ----------------------'


Da_upper={}
for A in range(1,3):
    Da_upper[A]=SUSY_zeros()
    for i in range(1,3):
        for I in range(2):
            Da_upper[A]+=SUSY(FZ12*Num(ref*RA[i][:,I]),[str(i)+str(I)+str(A)])
for A in range(3,5):
    Da_upper[A]=SUSY_zeros()
    for i in range(1,3):
        for I in range(2):
            Da_upper[A]+=SUSY(FZ34*Num(ref*RA[i][:,I]),[str(i)+str(I)+str(A)])
for A in range(5,7):
    Da_upper[A]=SUSY_zeros()
    for i in range(1,3):
        for I in range(2):
            Da_upper[A]+=SUSY(FZ56*Num(ref*RA[i][:,I]),[str(i)+str(I)+str(A)])
for A in range(7,9):
    Da_upper[A]=SUSY_zeros()
    for i in range(1,3):
        for I in range(2):
            Da_upper[A]+=SUSY(FZ78*Num(ref*RA[i][:,I]),[str(i)+str(I)+str(A)])



Db_upper={}
for A in range(1,3):
    Db_upper[A]=SUSY_zeros()
    for B in range(1,3):
        for i in range(1,3):
            for I in range(2):
                Db_upper[A]+=SUSY((-GZ12)*Num(ref*RB[i][:,I])*sig[i]*epsilonL(A,B),[str(i)+str(I)+str(B)])
for A in range(3,5):
    Db_upper[A]=SUSY_zeros()
    for B in range(3,5):
        for i in range(1,3):
            for I in range(2):
                Db_upper[A]+=SUSY((-GZ34)*Num(ref*RB[i][:,I])*sig[i]*epsilonL(A,B),[str(i)+str(I)+str(B)])
for A in range(5,7):
    Db_upper[A]=SUSY_zeros()
    for B in range(5,7):
        for i in range(1,3):
            for I in range(2):
                Db_upper[A]+=SUSY((-GZ56)*Num(ref*RB[i][:,I])*sig[i]*epsilonL(A,B),[str(i)+str(I)+str(B)])
for A in range(7,9):
    Db_upper[A]=SUSY_zeros()
    for B in range(7,9):
        for i in range(1,3):
            for I in range(2):
                Db_upper[A]+=SUSY((-GZ78)*Num(ref*RB[i][:,I])*sig[i]*epsilonL(A,B),[str(i)+str(I)+str(B)])


Da_lower={}
for A in range(1,3):
    Da_lower[A]=Diff_zeros()
    for i in range(1,3):
        for I in range(2):
            for J in range(2):
                for B in range(1,3):
                    Da_lower[A]+=Diff(GZ12*sig[i]*Num(ref*RA[i][:,I])*epsilonU(A,B)*epsilon_L[I,J],str(i)+str(J)+str(B))
for A in range(3,5):
    Da_lower[A]=Diff_zeros()
    for i in range(1,3):
        for I in range(2):
            for J in range(2):
                for B in range(3,5):
                    Da_lower[A]+=Diff(GZ34*sig[i]*Num(ref*RA[i][:,I])*epsilonU(A,B)*epsilon_L[I,J],str(i)+str(J)+str(B))
for A in range(5,7):
    Da_lower[A]=Diff_zeros()
    for i in range(1,3):
        for I in range(2):
            for J in range(2):
                for B in range(5,7):
                    Da_lower[A]+=Diff(GZ56*sig[i]*Num(ref*RA[i][:,I])*epsilonU(A,B)*epsilon_L[I,J],str(i)+str(J)+str(B))
for A in range(7,9):
    Da_lower[A]=Diff_zeros()
    for i in range(1,3):
        for I in range(2):
            for J in range(2):
                for B in range(7,9):
                    Da_lower[A]+=Diff(GZ78*sig[i]*Num(ref*RA[i][:,I])*epsilonU(A,B)*epsilon_L[I,J],str(i)+str(J)+str(B))


Db_lower={}
for A in range(1,3):
    Db_lower[A]=Diff_zeros()
    for i in range(1,3):
        for I in range(2):
            for J in range(2):
                Db_lower[A]+=Diff(-FZ12*Num(ref*RB[i][:,I])*epsilon_L[I,J],str(i)+str(J)+str(A))
for A in range(3,5):
    Db_lower[A]=Diff_zeros()
    for i in range(1,3):
        for I in range(2):
            for J in range(2):
                Db_lower[A]+=Diff(-FZ34*Num(ref*RB[i][:,I])*epsilon_L[I,J],str(i)+str(J)+str(A))
for A in range(5,7):
    Db_lower[A]=Diff_zeros()
    for i in range(1,3):
        for I in range(2):
            for J in range(2):
                Db_lower[A]+=Diff(-FZ56*Num(ref*RB[i][:,I])*epsilon_L[I,J],str(i)+str(J)+str(A))
for A in range(7,9):
    Db_lower[A]=Diff_zeros()
    for i in range(1,3):
        for I in range(2):
            for J in range(2):
                Db_lower[A]+=Diff(-FZ78*Num(ref*RB[i][:,I])*epsilon_L[I,J],str(i)+str(J)+str(A))


Da={}
Db={}
QA={}
QB={}
for A in range(1,9):
    Da[A]=Da_upper[A]+Da_lower[A]
    Db[A]=Db_upper[A]+Db_lower[A]
    QA[A]=Da[A]+Num(ref*RA[0])*uSUSY('0'+str(A))
    QB[A]=Db[A]+Num(ref*RB[0])*uDiff('0'+str(A))





'---------------------- building blocks ----------------------'



delta={}
for A in range(1,9):
    delta[A]=SUSY_zeros()
    for i in range(1,3):
        for I in range(2):
            delta[A]+=SUSY(Num(LA[0]*RA[i][:,I]),[str(i)+str(I)+str(A)])
delta12=delta[1]*delta[2]
delta34=delta[3]*delta[4]
delta56=delta[5]*delta[6]
delta78=delta[7]*delta[8]




etaP={}
etaP[1]=uSUSY('01')
etaP[2]=uSUSY('02')
etaP[3]=uSUSY('03')
etaP[4]=uSUSY('04')
etaP[5]=uSUSY('05')
etaP[6]=uSUSY('06')
etaP[7]=uSUSY('07')
etaP[8]=uSUSY('08')
etaPs12=0
for (A,B) in product(*([range(1,3)]*2)):
    etaPs12+=(1/2)*epsilonL(A,B)*etaP[A]*etaP[B]
etaPs34=0
for (A,B) in product(*([range(3,5)]*2)):
    etaPs34+=(1/2)*epsilonL(A,B)*etaP[A]*etaP[B]
etaPs56=0
for (A,B) in product(*([range(5,7)]*2)):
    etaPs56+=(1/2)*epsilonL(A,B)*etaP[A]*etaP[B]
etaPs78=0
for (A,B) in product(*([range(7,9)]*2)):
    etaPs78+=(1/2)*epsilonL(A,B)*etaP[A]*etaP[B]




zeta=np.matrix(np.random.randn(1,2)+np.ones((1,2)))
xi=np.matrix(np.random.randn(1,2)+np.ones((1,2)))



Ha={}
for A in range(1,9):
    Ha[A]=SUSY_zeros()
    for i in range(1,3):
        for I in range(2):
            Ha[A]+=SUSY(Num(zeta*RA[i][:,I])/Num(zeta*RA[0]),[str(i)+str(I)+str(A)])



Hb={}
for A in range(1,9):
    Hb[A]=SUSY_zeros()
    for B in range(1,9):
        for i in range(1,3):
            for I in range(2):
                Hb[A]+=SUSY(Num(xi*RB[i][:,I])*sig[i]*epsilonL(A,B,False)/Num(xi*RB[0]),[str(i)+str(I)+str(B)])




HHa12=-Ha[1]*Ha[2]
HHa34=-Ha[3]*Ha[4]
HHa56=-Ha[5]*Ha[6]
HHa78=-Ha[7]*Ha[8]
HHb12=Hb[1]*Hb[2]
HHb34=Hb[3]*Hb[4]
HHb56=Hb[5]*Hb[6]
HHb78=Hb[7]*Hb[8]
HHab12=(1/2)*(Ha[1]*Hb[1]+Ha[2]*Hb[2])
HHab34=(1/2)*(Ha[3]*Hb[3]+Ha[4]*Hb[4])
HHab56=(1/2)*(Ha[5]*Hb[5]+Ha[6]*Hb[6])
HHab78=(1/2)*(Ha[7]*Hb[7]+Ha[8]*Hb[8])


'---------------------- Super amplitude ----------------------'



solution12={}

solution12[1]=0

solution12[1]+=delta12*etaPs12*(1+(-2*GZ12/FZ12)*HHab12+(GZ12**2/FZ12**2)*HHa12*HHb12)
for (A,B) in product(*([range(1,3)]*2)):
    solution12[1]+=delta12*etaP[A]*(1/FZ12)*epsilonL(A,B)*Ha[B]
for A in range(1,3):
    solution12[1]+=delta12*etaP[A]*(GZ12/(FZ12**2))*HHa12*Hb[A]
solution12[1]+=delta12*(1/(FZ12**2))*HHa12

solution12[1]*=FZ12**2


solution12[2]=0

solution12[2]+=delta12*etaPs12*HHb12
for (A,B) in product(*([range(1,3)]*2)):
    solution12[2]+=delta12*etaP[A]*FZ12*epsilonL(A,B)*Ha[B]*HHb12
for A in range(1,3):
    solution12[2]+=delta12*etaP[A]*(GZ12)*Hb[A]
solution12[2]+=delta12*((GZ12**2)+(2*FZ12*GZ12)*HHab12+(FZ12**2)*HHa12*HHb12)

solution12[2]*=-x12







solution34={}

solution34[1]=0

solution34[1]+=delta34*etaPs34*(1+(-2*GZ34/FZ34)*HHab34+(GZ34**2/FZ34**2)*HHa34*HHb34)
for (A,B) in product(*([range(3,5)]*2)):
    solution34[1]+=delta34*etaP[A]*(1/FZ34)*epsilonL(A,B)*Ha[B]
for A in range(3,5):
    solution34[1]+=delta34*etaP[A]*(GZ34/(FZ34**2))*HHa34*Hb[A]
solution34[1]+=delta34*(1/(FZ34**2))*HHa34

solution34[1]*=FZ34**2


solution34[2]=0

solution34[2]+=delta34*etaPs34*HHb34
for (A,B) in product(*([range(3,5)]*2)):
    solution34[2]+=delta34*etaP[A]*FZ34*epsilonL(A,B)*Ha[B]*HHb34
for A in range(3,5):
    solution34[2]+=delta34*etaP[A]*(GZ34)*Hb[A]
solution34[2]+=delta34*((GZ34**2)+(2*FZ34*GZ34)*HHab34+(FZ34**2)*HHa34*HHb34)

solution34[2]*=-x12



solution56={}

solution56[1]=0

solution56[1]+=delta56*etaPs56*(1+(-2*GZ56/FZ56)*HHab56+(GZ56**2/FZ56**2)*HHa56*HHb56)
for (A,B) in product(*([range(5,7)]*2)):
    solution56[1]+=delta56*etaP[A]*(1/FZ56)*epsilonL(A,B)*Ha[B]
for A in range(5,7):
    solution56[1]+=delta56*etaP[A]*(GZ56/(FZ56**2))*HHa56*Hb[A]
solution56[1]+=delta56*(1/(FZ56**2))*HHa56

solution56[1]*=FZ56**2


solution56[2]=0

solution56[2]+=delta56*etaPs56*HHb56
for (A,B) in product(*([range(5,7)]*2)):
    solution56[2]+=delta56*etaP[A]*FZ56*epsilonL(A,B)*Ha[B]*HHb56
for A in range(5,7):
    solution56[2]+=delta56*etaP[A]*(GZ56)*Hb[A]
solution56[2]+=delta56*((GZ56**2)+(2*FZ56*GZ56)*HHab56+(FZ56**2)*HHa56*HHb56)

solution56[2]*=-x12




solution78={}

solution78[1]=0

solution78[1]+=delta78*etaPs78*(1+(-2*GZ78/FZ78)*HHab78+(GZ78**2/FZ78**2)*HHa78*HHb78)
for (A,B) in product(*([range(7,9)]*2)):
    solution78[1]+=delta78*etaP[A]*(1/FZ78)*epsilonL(A,B)*Ha[B]
for A in range(7,9):
    solution78[1]+=delta78*etaP[A]*(GZ78/(FZ78**2))*HHa78*Hb[A]
solution78[1]+=delta78*(1/(FZ78**2))*HHa78

solution78[1]*=FZ78**2


solution78[2]=0

solution78[2]+=delta78*etaPs78*HHb78
for (A,B) in product(*([range(7,9)]*2)):
    solution78[2]+=delta78*etaP[A]*FZ78*epsilonL(A,B)*Ha[B]*HHb78
for A in range(7,9):
    solution78[2]+=delta78*etaP[A]*(GZ78)*Hb[A]
solution78[2]+=delta78*((GZ78**2)+(2*FZ78*GZ78)*HHab78+(FZ78**2)*HHa78*HHb78)

solution78[2]*=-x12









'---------------------- verify ----------------------'




accuarcy=1e-8



print((
len(QA[1].operate(solution12[1]).estimate(accuarcy).arr)
+len(QA[1].operate(solution12[2]).estimate(accuarcy).arr)
+len(QA[2].operate(solution12[1]).estimate(accuarcy).arr)
+len(QA[2].operate(solution12[2]).estimate(accuarcy).arr)
+len(QA[3].operate(solution34[1]).estimate(accuarcy).arr)
+len(QA[3].operate(solution34[2]).estimate(accuarcy).arr)
+len(QA[4].operate(solution34[1]).estimate(accuarcy).arr)
+len(QA[4].operate(solution34[2]).estimate(accuarcy).arr)
+len(QB[1].operate(solution12[1]).estimate(accuarcy).arr)
+len(QB[1].operate(solution12[2]).estimate(accuarcy).arr)
+len(QB[2].operate(solution12[1]).estimate(accuarcy).arr)
+len(QB[2].operate(solution12[2]).estimate(accuarcy).arr)
+len(QB[3].operate(solution34[1]).estimate(accuarcy).arr)
+len(QB[3].operate(solution34[2]).estimate(accuarcy).arr)
+len(QB[4].operate(solution34[1]).estimate(accuarcy).arr)
+len(QB[4].operate(solution34[2]).estimate(accuarcy).arr))
==0)




'---------------------- previous ----------------------'


