import Supersymmetry_BPS_3pt as BPS
import Supersymmetry_nonBPS_3pt as nBPS
import numpy as np
from numpy import sqrt
print('FZ=%f, GZ=%f' % (nBPS.FZ, nBPS.GZ))


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


if __name__ == '__main__':
    print((BPS.solution[1]*BPS.solution[2] -
          etaP_Mul(nBPS.solution[1])).vanish())
