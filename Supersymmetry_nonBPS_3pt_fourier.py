import numpy as np
from Grassmannian import*
from SpinorBrackets3pt import*
from itertools import product
from Supersymmetry_nonBPS_3pt import FZ, GZ, solution


'---------------------- building blocks ----------------------'


deltaf = {}
for A in range(1, 3):
    deltaf[A] = SUSY_zeros()
    for i in range(1, 3):
        for I in range(2):
            deltaf[A] += SUSY(Num(LB[0]*RB[i][:, I]),
                              [str(i)+str(I)+str(A)+'f'])
deltaf12 = deltaf[1]*deltaf[2]


etaPf = {}
etaPf[1] = uSUSY('01f')
etaPf[2] = uSUSY('02f')


etaPs = 0
for (A, B) in product(*([range(1, 3)]*2)):
    etaPs += (1/2)*epsilonL(A, B)*etaPf[A]*etaPf[B]


zeta = np.matrix(np.random.randn(1, 2)+np.ones((1, 2)))
xi = np.matrix(np.random.randn(1, 2)+np.ones((1, 2)))


Haf = {}
for A in range(1, 3):
    Haf[A] = SUSY_zeros()
    for i in range(1, 3):
        for I in range(2):
            Haf[A] += SUSY(Num(zeta*RB[i][:, I])/Num(zeta*RB[0]),
                           [str(i)+str(I)+str(A)+'f'])


Hbf = {}
for A in range(1, 3):
    Hbf[A] = SUSY_zeros()
    for B in range(1, 3):
        for i in range(1, 3):
            for I in range(2):
                Hbf[A] -= SUSY(Num(xi*RA[i][:, I])*sig[i]*epsilonL(A,
                               B)/Num(xi*RA[0]), [str(i)+str(I)+str(B)+'f'])


HHaf = -Haf[1]*Haf[2]
HHbf = Hbf[1]*Hbf[2]
HHabf = (1/2)*(Haf[1]*Hbf[1]+Haf[2]*Hbf[2])


'---------------------- Super amplitude ----------------------'


solutionf = {}


solutionf[1] = 0

solutionf[1] += deltaf12*(etaPs)*(1+(-2*GZ/FZ)*HHabf+(GZ**2/FZ**2)*HHaf*HHbf)
for (A, B) in product(*([range(1, 3)]*2)):
    solutionf[1] += deltaf12*etaPf[A]*(1/FZ)*epsilonL(A, B)*Haf[B]
for A in range(1, 3):
    solutionf[1] += deltaf12*etaPf[A]*(GZ/(FZ**2))*HHaf*Hbf[A]
solutionf[1] += deltaf12*(1/(FZ**2))*HHaf

solutionf[1] *= FZ**2


solutionf[2] = 0

solutionf[2] += deltaf12*(etaPs)*HHbf
for (A, B) in product(*([range(1, 3)]*2)):
    solutionf[2] += deltaf12*etaPf[A]*FZ*epsilonL(A, B)*Haf[B]*HHbf
for A in range(1, 3):
    solutionf[2] += deltaf12*etaPf[A]*(GZ)*Hbf[A]
solutionf[2] += deltaf12*((GZ**2)+(2*FZ*GZ)*HHabf+(FZ**2)*HHaf*HHbf)


solutionf[2] /= x12


'---------------------- Fourier transform ----------------------'


def uSUSYF(*a):
    return SUSY(1, [''.join(str(i) for i in a)])


def POWERft(a, n):
    if (not type(n) == int) or n < 0:
        print('error: the power should be an positive integer.')
        return None
    elif n == 0:
        return 1
    elif n == 1:
        return a
    else:
        return POWERft(a, n-1)*a


def fac(n):
    if (not type(n) == int) or n < 0:
        print('error: the power should be an positive integer.')
        return None
    elif n == 0:
        return 1
    elif n == 1:
        return 1
    else:
        return fac(n-1)*n


power_exp = 4
print('The order of exponentiation is', power_exp)


def EXPft(x, power_max=power_exp):
    temp_power = 0
    result = POWERft(x, temp_power)
    temp_power += 1
    while temp_power <= power_max:
        result += POWERft(x, temp_power)/fac(temp_power)
        temp_power += 1
    return result


exponent = 1
all_eta = 0
for (A, B) in product(*([range(1, 3)]*2)):
    all_eta += epsilonL(A, B)*uSUSYF(0, A, 'f')*uSUSYF(0, B)
exponent *= EXPft(all_eta)
for i in (1, 2):
    all_eta = 0
    for (A, B) in product(*([range(1, 3)]*2)):
        for (I, J) in product(*([range(2)]*2)):
            all_eta += epsilon_U[I, J] * \
                epsilonL(A, B)*uSUSYF(i, I, A, 'f')*uSUSYF(i, J, B)
    exponent *= EXPft(all_eta)
del (all_eta)


def FT(amplitude):

    etaList = [[1, ['01', '02', '101', '111',
                    '102', '112', '201', '211', '202', '212']]]
    result = diff_etalist(exponent*amplitude, etaList)

    return -result*x12


exponent1to2 = 1
all_eta = 0
for (A, B) in product(*([range(1, 3)]*2)):
    all_eta += epsilonL(A, B)*uSUSYF(0, A, 'f')*uSUSYF(0, B)
exponent1to2 *= EXPft(all_eta)
all_eta = 0
for i in (1, 2):
    all_eta = 0
    for (I, J) in product(*([range(2)]*2)):
        all_eta += epsilon_U[I, J]*uSUSYF(i, I, 2, 'f')*uSUSYF(i, J, 1)
    exponent1to2 *= EXPft(all_eta)
del (all_eta)


def FT1to2(amplitude):

    etaList = [[1, ['01', '02', '101', '111', '201', '211']]]
    result = diff_etalist(exponent1to2*amplitude, etaList)

    return result*x12


exponent2to1 = 1
all_eta = 0
for (A, B) in product(*([range(1, 3)]*2)):
    all_eta += epsilonL(A, B)*uSUSYF(0, A, 'f')*uSUSYF(0, B)
exponent2to1 *= EXPft(all_eta)
all_eta = 0
for i in (1, 2):
    all_eta = 0
    for (I, J) in product(*([range(2)]*2)):
        all_eta += epsilon_U[I, J]*uSUSYF(i, I, 1, 'f')*uSUSYF(i, J, 2)
    exponent2to1 *= EXPft(all_eta)
del (all_eta)


def FT2to1(amplitude):

    etaList = [[1, ['01', '02', '102', '112', '202', '212']]]
    result = diff_etalist(exponent2to1*amplitude, etaList)

    return result*x12


if __name__ == '__main__':
    print(len((FT(solution[1])-solutionf[2]).estimate(1e-5).arr) == 0)
    print(len((FT(solution[2])-solutionf[1]).estimate(1e-5).arr) == 0)

del solution
