import numpy as np
from Grassmannian import*
from SpinorBrackets3pt import*
from itertools import product
from numpy import sqrt

'---------------------- Examples ----------------------'

'''

Generators:
FZ , GZ
Da_upper[A] , Db_upper[A] , Da_lower[A] , Db_lower[A]
Da[A] , Db[A] , QA[A] , QB[A]

building blocks:
delta[A] , delta12
Ha[A] , Hb[A] , HHa , HHb , HHab

building blocks (Lorentz covariant):
etaP[A] , etas(i,A)
Qa[A] , Qb[A] , R[A]

Super amplitude:
solution_HaHb[#] , solution[#]

'''


'---------------------- Generators ----------------------'


# FZ=1
# GZ=float(sqrt(1-FZ**2))


FZ = 1/sqrt(2)
GZ = float(sqrt(1-FZ**2))


# FZ=0.2
# GZ=float(sqrt(1-FZ**2))


Da_upper = {}
for A in range(1, 3):
    Da_upper[A] = SUSY_zeros()
    for i in range(1, 3):
        for I in range(2):
            Da_upper[A] += SUSY(FZ*Num(ref*RA[i][:, I]),
                                [str(i)+str(I)+str(A)])

Db_upper = {}
for A in range(1, 3):
    Db_upper[A] = SUSY_zeros()
    for B in range(1, 3):
        for i in range(1, 3):
            for I in range(2):
                Db_upper[A] += SUSY((-GZ)*Num(ref*RB[i][:, I])
                                    * sig[i]*epsilonL(A, B), [str(i)+str(I)+str(B)])

Da_lower = {}
for A in range(1, 3):
    Da_lower[A] = Diff_zeros()
    for i in range(1, 3):
        for I in range(2):
            for J in range(2):
                for B in range(1, 3):
                    Da_lower[A] += Diff(GZ*sig[i]*Num(ref*RA[i][:, I]) *
                                        epsilonU(A, B)*epsilon_L[I, J], str(i)+str(J)+str(B))

Db_lower = {}
for A in range(1, 3):
    Db_lower[A] = Diff_zeros()
    for i in range(1, 3):
        for I in range(2):
            for J in range(2):
                Db_lower[A] += Diff(-FZ*Num(ref*RB[i][:, I])
                                    * epsilon_L[I, J], str(i)+str(J)+str(A))


Da = {}
Db = {}
QA = {}
QB = {}
for A in range(1, 3):
    Da[A] = Da_upper[A]+Da_lower[A]
    Db[A] = Db_upper[A]+Db_lower[A]
    QA[A] = Da[A]+Num(ref*RA[0])*uSUSY('0'+str(A))
    QB[A] = Db[A]+Num(ref*RB[0])*uDiff('0'+str(A))


'---------------------- building blocks ----------------------'


delta = {}
for A in range(1, 3):
    delta[A] = SUSY_zeros()
    for i in range(1, 3):
        for I in range(2):
            delta[A] += SUSY(Num(LA[0]*RA[i][:, I]), [str(i)+str(I)+str(A)])
delta12 = delta[1]*delta[2]


etaP = {}
etaP[1] = uSUSY('01')
etaP[2] = uSUSY('02')


etaPs = 0
for (A, B) in product(*([range(1, 3)]*2)):
    etaPs += (1/2)*epsilonL(A, B)*etaP[A]*etaP[B]


zeta = np.matrix(np.random.randn(1, 2)+np.ones((1, 2)))
xi = np.matrix(np.random.randn(1, 2)+np.ones((1, 2)))


Ha = {}
for A in range(1, 3):
    Ha[A] = SUSY_zeros()
    for i in range(1, 3):
        for I in range(2):
            Ha[A] += SUSY(Num(zeta*RA[i][:, I])/Num(zeta*RA[0]),
                          [str(i)+str(I)+str(A)])


Hb = {}
for A in range(1, 3):
    Hb[A] = SUSY_zeros()
    for B in range(1, 3):
        for i in range(1, 3):
            for I in range(2):
                Hb[A] += SUSY(Num(xi*RB[i][:, I])*sig[i]*epsilonL(A,
                              B)/Num(xi*RB[0]), [str(i)+str(I)+str(B)])


HHa = -Ha[1]*Ha[2]
HHb = Hb[1]*Hb[2]
HHab = (1/2)*(Ha[1]*Hb[1]+Ha[2]*Hb[2])


'---------------------- Super amplitude ----------------------'


solution = {}


solution[1] = 0

solution[1] += delta12*(etaPs)*(1+(-2*GZ/FZ)*HHab+(GZ**2/FZ**2)*HHa*HHb)
for (A, B) in product(*([range(1, 3)]*2)):
    solution[1] += delta12*etaP[A]*(1/FZ)*epsilonL(A, B)*Ha[B]
for A in range(1, 3):
    solution[1] += delta12*etaP[A]*(GZ/(FZ**2))*HHa*Hb[A]
solution[1] += delta12*(1/(FZ**2))*HHa

solution[1] *= FZ**2


solution[2] = 0

solution[2] += delta12*(etaPs)*HHb
for (A, B) in product(*([range(1, 3)]*2)):
    solution[2] += delta12*etaP[A]*(FZ)*epsilonL(A, B)*Ha[B]*HHb
for A in range(1, 3):
    solution[2] += delta12*etaP[A]*(GZ)*Hb[A]
solution[2] += delta12*((GZ**2)+(2*FZ*GZ)*HHab+(FZ**2)*HHa*HHb)


solution[2] *= -x12


'---------------------- verify ----------------------'

if __name__ == '__main__':
    precision = 1e-7

    print((
        len(QA[1].operate(solution[1]).estimate(precision).arr)
        + len(QA[2].operate(solution[1]).estimate(precision).arr)
        + len(QB[1].operate(solution[1]).estimate(precision).arr)
        + len(QB[2].operate(solution[1]).estimate(precision).arr)
        + len(QA[1].operate(solution[2]).estimate(precision).arr)
        + len(QA[2].operate(solution[2]).estimate(precision).arr)
        + len(QB[1].operate(solution[2]).estimate(precision).arr)
        + len(QB[2].operate(solution[2]).estimate(precision).arr))
        == 0)
