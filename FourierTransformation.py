from Supersymmetry_nonBPS_3pt import*
from Supersymmetry_nonBPS_3pt_fourier import*


def etaP_Mul(amplitude, multiplier=sqrt(2)):
    result = 0
    a = multiplier
    for susy in amplitude.arr:
        temp = susy.copy()
        if '01' in susy.etalist:
            temp *= a
        if '02' in susy.etalist:
            temp *= a
        result += temp
    return result


def etaPf_Mul(amplitude, multiplier=sqrt(2)):
    result = 0
    a = multiplier
    for susy in amplitude.arr:
        temp = susy.copy()
        if '01f' in susy.etalist:
            temp *= a
        if '02f' in susy.etalist:
            temp *= a
        result += temp
    return result


factorized = {}
factorized[1] = {}
factorized[2] = {}

factorized[1][1] = FZ*etaP[1]*delta[1]-delta[1]*Ha[1] - \
    (FZ*GZ)*delta[1]*Hb[2]*etaP[1]*etaP[2]+GZ*delta[1]*Ha[1]*Hb[2]*etaP[2]
factorized[1][2] = FZ*etaP[2]*delta[2]-delta[2]*Ha[2] + \
    (FZ*GZ)*delta[2]*Hb[1]*etaP[1]*etaP[2]+(GZ)*delta[2]*Ha[2]*Hb[1]*etaP[1]

factorized[2][1] = FZ*delta[1]*Ha[1]*Hb[1]+delta[1]*etaP[1]*Hb[1]+GZ*delta[1]
factorized[2][2] = FZ*delta[2]*Ha[2]*Hb[2]+delta[2]*etaP[2]*Hb[2]+GZ*delta[2]


factorizedf = {}
factorizedf[1] = {}
factorizedf[2] = {}

factorizedf[1][1] = FZ*etaPf[1]*deltaf[1]-deltaf[1]*Haf[1] - \
    (FZ*GZ)*deltaf[1]*Hbf[2]*etaPf[1] * \
    etaPf[2]+GZ*deltaf[1]*Haf[1]*Hbf[2]*etaPf[2]
factorizedf[1][2] = FZ*etaPf[2]*deltaf[2]-deltaf[2]*Haf[2] + \
    (FZ*GZ)*deltaf[2]*Hbf[1]*etaPf[1]*etaPf[2] + \
    (GZ)*deltaf[2]*Haf[2]*Hbf[1]*etaPf[1]

factorizedf[2][1] = FZ*deltaf[1]*Haf[1] * \
    Hbf[1]+deltaf[1]*etaPf[1]*Hbf[1]+GZ*deltaf[1]
factorizedf[2][2] = FZ*deltaf[2]*Haf[2] * \
    Hbf[2]+deltaf[2]*etaPf[2]*Hbf[2]+GZ*deltaf[2]


if __name__ == '__main__':
    print('FZ = {}\nGZ = {}'.format(FZ, GZ))
    print((factorized[1][1]*factorized[1][2]-solution[1]).vanish(1e-8))
    print(((-x12)*factorized[2][1]*factorized[2][2]-solution[2]).vanish(1e-8))
    print((etaPf_Mul(factorizedf[1][2]) -
          FT1to2(etaP_Mul(factorized[1][1]))).vanish(1e-8))
    print((etaPf_Mul(factorizedf[1][1], -sqrt(2)) -
          FT2to1(etaP_Mul(factorized[1][2]))).vanish(1e-8))
