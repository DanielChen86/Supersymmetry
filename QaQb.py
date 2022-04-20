from Supersymmetry_nonBPS_3pt import*

# write amplitude in terms of QaQb (instead of delta, Ha, Hb)

print_TF=True
list_TF=[]


Qa12={0:{1:None,2:None},1:{1:None,2:None}}
Qb12={0:{1:None,2:None},1:{1:None,2:None}}
QaP={0:{1:None,2:None},1:{1:None,2:None}}
QbP={0:{1:None,2:None},1:{1:None,2:None}}



def AALR(left,leftR,right,rightR):
    result=0
    for (alpha,beta) in product(*([range(2)]*2)):
        result+=(-epsilon_U[alpha,beta])*left[alpha][leftR]*right[beta][rightR]
    return result
def AAL(left,leftR,right):
    result=0
    for (alpha,beta) in product(*([range(2)]*2)):
        result+=(-epsilon_U[alpha,beta])*left[alpha][leftR]*right[beta,0]
    return result
def AAR(left,right,rightR):
    result=0
    for (alpha,beta) in product(*([range(2)]*2)):
        result+=(-epsilon_U[alpha,beta])*left[alpha,0]*right[beta][rightR]
    return result
def AA(left,right):
    result=0
    for (alpha,beta) in product(*([range(2)]*2)):
        result+=(-epsilon_U[alpha,beta])*left[alpha,0]*right[beta,0]
    return result

def BBLR(left,leftR,right,rightR):
    result=0
    for (alpha,beta) in product(*([range(2)]*2)):
        result+=(epsilon_U[alpha,beta])*left[alpha][leftR]*right[beta][rightR]
    return result
def BBL(left,leftR,right):
    result=0
    for (alpha,beta) in product(*([range(2)]*2)):
        result+=(epsilon_U[alpha,beta])*left[alpha][leftR]*right[0,beta]
    return result
def BBR(left,right,rightR):
    result=0
    for (alpha,beta) in product(*([range(2)]*2)):
        result+=(epsilon_U[alpha,beta])*left[0,alpha]*right[beta][rightR]
    return result
def BB(left,right):
    result=0
    for (alpha,beta) in product(*([range(2)]*2)):
        result+=(epsilon_U[alpha,beta])*left[0,alpha]*right[0,beta]
    return result





'---------------------- building blocks ----------------------'


for alpha in range(2):
    for A in range(1,3):
        amplitude=0
        for I in range(2):
            amplitude+=SUSY(RA[1][alpha,I],['1'+str(I)+str(A)])
            amplitude+=SUSY(RA[2][alpha,I],['2'+str(I)+str(A)])
        Qa12[alpha][A]=amplitude




for alpha in range(2):
    for A in range(1,3):
        amplitude=0
        for I in range(2):
            amplitude+=SUSY(LB[1][I,alpha],['1'+str(I)+str(A)])
            amplitude-=SUSY(LB[2][I,alpha],['2'+str(I)+str(A)])
        Qb12[alpha][A]=amplitude

for alpha in range(2):
    for A in range(1,3):
        QaP[alpha][A]=SUSY(RA[0][alpha,0],['0'+str(A)])

for alpha in range(2):
    for A in range(1,3):
        QbP[alpha][A]=SUSY(LB[0][0,alpha],['0'+str(A)])


ETAs={}
for A in range(1,3):
    for B in range(1,3):
        amplitude=0
        for I in range(2):
            for J in range(2):
                amplitude-=SUSY((1/2)*epsilon_U[I,J],['1'+str(I)+str(A),'1'+str(J)+str(B)])
                amplitude+=SUSY((1/2)*epsilon_U[I,J],['2'+str(I)+str(A),'2'+str(J)+str(B)])
        ETAs[A,B]=amplitude





HaEta=0
for (A,B) in product(*([range(1,3)]*2)):
    HaEta+=Ha[A]*epsilonL(A,B)*etaP[B]


EtaHb=0
for A in range(1,3):
    EtaHb+=etaP[A]*Hb[A]




'---------------------- identities ----------------------'


precision=1e-8

if print_TF:

    amplitude=0
    for (A,B) in product(*([range(1,3)]*2)):
        amplitude+=(-1/2)*epsilonL(A,B)*AAL(Qa12,A,RA[0])*AAL(Qa12,B,RA[0])

    list_TF.append( (delta12-amplitude).estimate(precision).size==0 )



    amplitude=0
    for (A,B,C,D) in product(*([range(1,3)]*4)):
        amplitude+=(1/4)*epsilonL(A,C)*epsilonL(B,D)*AALR(Qa12,A,QaP,B)*AALR(Qa12,C,QaP,D)

    list_TF.append( (delta12*etaPs-amplitude).estimate(precision).size==0 )



    amplitude=0
    for (A,B,C,D) in product(*([range(1,3)]*4)):
        amplitude+=(1/12)*epsilonL(A,C)*epsilonL(B,D)*AALR(Qa12,A,Qa12,B)*AALR(Qa12,C,Qa12,D)

    list_TF.append( (delta12*HHa-amplitude).estimate(precision).size==0 )


    amplitude=0
    for (A,B,C,D) in product(*([range(1,3)]*4)):
        amplitude+=(-1/12/x12)*epsilonL(A,C)*epsilonL(B,D)*AALR(Qa12,A,Qa12,B)*BBLR(Qb12,C,Qb12,D)

    list_TF.append( (delta12*HHab-amplitude).estimate(precision).size==0 )




    amplitude=0
    for (A,B,C,D) in product(*([range(1,3)]*4)):
        amplitude+=(-1/12/x12**2)*epsilonL(A,C)*epsilonL(B,D)*BBLR(Qb12,A,Qb12,B)*BBLR(Qb12,C,Qb12,D)

    list_TF.append( (delta12*HHb-amplitude).estimate(precision).size==0 )




    amplitude=0
    for (A,B,C,D) in product(*([range(1,3)]*4)):
        amplitude+=(1/3)*epsilonL(A,C)*epsilonL(B,D)*AALR(Qa12,A,Qa12,B)*AALR(Qa12,C,QaP,D)

    list_TF.append( (delta12*HaEta-amplitude).estimate(precision).size==0 )


    amplitude=0
    for (A,B,C,D) in product(*([range(1,3)]*4)):
        amplitude+=(-1/3/x12)*epsilonL(A,C)*epsilonL(B,D)*BBLR(Qb12,A,Qb12,B)*AALR(Qa12,C,QaP,D)

    list_TF.append( (delta12*EtaHb-amplitude).estimate(precision).size==0 )





    amplitude=0
    for (A,B,C,D,E,F) in product(*([range(1,3)]*6)):
        amplitude+=(-2/9/x12)*epsilonL(A,D)*epsilonL(B,E)*epsilonL(C,F)*AALR(Qa12,A,Qa12,B)*BBLR(Qb12,D,QbP,C)*ETAs[E,F]

    list_TF.append( (delta12*HHa*EtaHb-amplitude).estimate(precision).size==0 )


    amplitude=0
    for (A,B,C,D,E,F) in product(*([range(1,3)]*6)):
        amplitude+=(2/9/x12**2)*epsilonL(A,D)*epsilonL(B,E)*epsilonL(C,F)*BBLR(Qb12,A,Qb12,B)*BBLR(Qb12,D,QbP,C)*ETAs[E,F]

    list_TF.append( (delta12*HHb*HaEta-amplitude).estimate(precision).size==0 )





    amplitude=0
    for (A,B,C,D,E,F) in product(*([range(1,3)]*6)):
        amplitude+=(-1/6/x12**2)*epsilonL(A,D)*epsilonL(B,E)*epsilonL(C,F)*BBL(Qb12,A,LB[0])*BBL(Qb12,D,LB[0])*ETAs[B,C]*ETAs[E,F]

    list_TF.append( (delta12*HHa*HHb-amplitude).estimate(precision).size==0 )



    amplitude=0
    for (A,B,C,D,E,F,G,H) in product(*([range(1,3)]*8)):
        amplitude+=(1/12/x12**2)*epsilonL(A,E)*epsilonL(B,F)*epsilonL(C,G)*epsilonL(D,H)*BBLR(Qb12,A,QbP,B)*BBLR(Qb12,E,QbP,F)*ETAs[C,D]*ETAs[G,H]

    list_TF.append( (delta12*etaPs*HHa*HHb-amplitude).estimate(precision).size==0 )






'---------------------- Super amplitude ----------------------'

solutionQ={}

solutionQ[1]=0

for (A,B,C,D) in product(*([range(1,3)]*4)):
    solutionQ[1]+=(FZ**2/4)*epsilonL(A,C)*epsilonL(B,D)*AALR(Qa12,A,QaP,B)*AALR(Qa12,C,QaP,D)

for (A,B,C,D) in product(*([range(1,3)]*4)):
    solutionQ[1]+=(FZ*GZ/6/x12)*epsilonL(A,C)*epsilonL(B,D)*AALR(Qa12,A,Qa12,B)*BBLR(Qb12,C,Qb12,D)*etaPs

for (A,B,C,D,E,F,G,H) in product(*([range(1,3)]*8)):
    solutionQ[1]+=(GZ**2/12/x12**2)*epsilonL(A,E)*epsilonL(B,F)*epsilonL(C,G)*epsilonL(D,H)*BBLR(Qb12,A,QbP,B)*BBLR(Qb12,E,QbP,F)*ETAs[C,D]*ETAs[G,H]

for (A,B,C,D) in product(*([range(1,3)]*4)):
    solutionQ[1]+=(FZ/3)*epsilonL(A,C)*epsilonL(B,D)*AALR(Qa12,A,Qa12,B)*AALR(Qa12,C,QaP,D)

for (A,B,C,D,E,F) in product(*([range(1,3)]*6)):
    solutionQ[1]+=(-2*GZ/9/x12)*epsilonL(A,D)*epsilonL(B,E)*epsilonL(C,F)*AALR(Qa12,A,Qa12,B)*BBLR(Qb12,D,QbP,C)*ETAs[E,F]

for (A,B,C,D) in product(*([range(1,3)]*4)):
    solutionQ[1]+=(1/12)*epsilonL(A,C)*epsilonL(B,D)*AALR(Qa12,A,Qa12,B)*AALR(Qa12,C,Qa12,D)





solutionQ[2]=0

for (A,B,C,D) in product(*([range(1,3)]*4)):
    solutionQ[2]+=(1/12/x12)*epsilonL(A,C)*epsilonL(B,D)*BBLR(Qb12,A,Qb12,B)*BBLR(Qb12,C,Qb12,D)*etaPs

for (A,B,C,D,E,F) in product(*([range(1,3)]*6)):
    solutionQ[2]+=(-2*FZ/9/x12)*epsilonL(A,D)*epsilonL(B,E)*epsilonL(C,F)*BBLR(Qb12,A,Qb12,B)*BBLR(Qb12,D,QbP,C)*ETAs[E,F]

for (A,B,C,D) in product(*([range(1,3)]*4)):
    solutionQ[2]+=(GZ/3)*epsilonL(A,C)*epsilonL(B,D)*BBLR(Qb12,A,Qb12,B)*AALR(Qa12,C,QaP,D)

for (A,B) in product(*([range(1,3)]*2)):
    solutionQ[2]+=(GZ**2*x12/2)*epsilonL(A,B)*AAL(Qa12,A,RA[0])*AAL(Qa12,B,RA[0])

for (A,B,C,D) in product(*([range(1,3)]*4)):
    solutionQ[2]+=(FZ*GZ/6)*epsilonL(A,C)*epsilonL(B,D)*AALR(Qa12,A,Qa12,B)*BBLR(Qb12,C,Qb12,D)

for (A,B,C,D,E,F) in product(*([range(1,3)]*6)):
    solutionQ[2]+=(FZ**2/6/x12)*epsilonL(A,D)*epsilonL(B,E)*epsilonL(C,F)*BBL(Qb12,A,LB[0])*BBL(Qb12,D,LB[0])*ETAs[B,C]*ETAs[E,F]




'---------------------- verify ----------------------'

if __name__=='__main__':

    precision=1e-7

    print (all(list_TF))
    
    print((
    len(QA[1].operate(solutionQ[1]).estimate(precision).arr)\
    +len(QA[2].operate(solutionQ[1]).estimate(precision).arr)\
    +len(QB[1].operate(solutionQ[1]).estimate(precision).arr)\
    +len(QB[2].operate(solutionQ[1]).estimate(precision).arr)\
    +len(QA[1].operate(solutionQ[2]).estimate(precision).arr)\
    +len(QA[2].operate(solutionQ[2]).estimate(precision).arr)\
    +len(QB[1].operate(solutionQ[2]).estimate(precision).arr)\
    +len(QB[2].operate(solutionQ[2]).estimate(precision).arr))\
    ==0)




