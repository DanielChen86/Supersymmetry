
'---------------------- Examples ----------------------'

'''

Eta([]) , Eta(['111']) , Eta(['101','202','01'])
Eta.etalist , Eta.sign

SUSY(1.,[]) , SUSY(3.+2.j,['111','112'])
SUSY.number , SUSY.etalist , SUSY.order
SUSY.copy() , SUSY.vanish(1e-6)

addSUSY([SUSY,SUSY,SUSY,...])
addSUSY.arr , addSUSY.alleta , addSUSY.size
addSUSY.get_order(3) , addSUSY.get(['101','01']) , addSUSY.contain(['01','02'])
addSUSY.estimate(1e-6) , addSUSY.vanish(1e-6) , addSUSY.copy()

diff(3.,'111') , diff(SUSY,'111')

Diff(3.,'111')
Diff.number , Diff.differentials
Diff.diff(SUSY/addSUSY) , Diff.copy()

addDiff([Diff,Diff,Diff,...])
addDiff.arr , addDiff.alldiff
addDiff.diff(SUSY/addSUSY) , addDiff.copy() , addDiff.estimate(1e-6)

SUSY_zeros() , Diff_zeros()
uSUSY('01') , uDiff('01')

Operator(multi,diff)    # multi=SUSY/addSUSY , diff=Diff/addDiff
Operator.operate(SUSY/addSUSY) , Operator.estimate(1e-6)
Operator = SUSY/addSUSY + Diff/addDiff
Operator_zeros()
Commute(Operator1,Operator2)
AntiCommute(Operator1,Operator2)

'''

'---------------------- Grassmannian numbers ----------------------'


validity=1e-10

# itos={0:'01',1:'02',
#       2:'101',3:'111',4:'102',5:'112',
#       6:'201',7:'211',8:'202',9:'212'}


itos={0:'01',1:'02',2:'03',3:'04',
      4:'101',5:'111',6:'102',7:'112',8:'103',9:'113',10:'104',11:'114',
      12:'201',13:'211',14:'202',15:'212',16:'203',17:'213',18:'204',19:'214',
      20:'301',21:'311',22:'302',23:'312',24:'303',25:'313',26:'304',27:'314',
      28:'401',29:'411',30:'402',31:'412',32:'403',33:'413',34:'404',35:'414'}



stoi={}
for key,value in itos.items():
    stoi[value]=key


def Fstoi(arr):
    result=[]
    for i in arr:
        result.append(stoi[i])
    return result

def Fitos(arr):
    result=[]
    for i in arr:
        result.append(itos[i])
    return result



def sort(arr_input):
    arr=arr_input.copy()
    power=0
    for i in range(len(arr)-1):
        minimum=arr[i]
        minimum_index=i
        for j in range(i,len(arr)):
            if arr[j]<minimum:
                minimum=arr[j]
                minimum_index=j
        arr=arr[:i]+[arr[minimum_index]]+arr[i:minimum_index]+arr[minimum_index+1:]
        power+=minimum_index-i
    return arr,power


class Eta():
    def __init__(self,etalist0):
        if (all(eta in itos.values() for eta in etalist0))==False:
            print ("error")
        self.etas0=etalist0
        self.__etanum=sort(Fstoi(etalist0))[0]
        self.sign=(-1)**(sort(Fstoi(etalist0))[1])
        self.etalist=Fitos(self.__etanum)
    
    def __repr__(self):
        return "Eta("+str(self.etalist)+",sign="+str(self.sign)+")"
    def __str__(self):
        return self.__repr__()

    def __mul__(self,other):
        if isinstance(other,Eta):
            return Eta(self.etas0+other.etas0)
        else:
            print ("error")


class SUSY(Eta):
    def __init__(self,number,etalist):
        super().__init__(etalist)
        self.number=number*self.sign

        if len(self.etalist)>len(set(self.etalist)):
            self.number=0.
        
        if self.number==0.:
            self.etalist=[]
        
        self.order=len(self.etalist)

        del self.sign
        del self.etas0

    def __repr__(self):
        return "SUSY("+str(complex(self.number))+","+str(self.etalist)+")"
    def __str__(self):
        return self.__repr__()
    
    def vanish(self,precision=validity):
        return bool(abs(complex(self.number))<=precision)


    def __neg__(self):
        return (-1)*self
    
    def __eq__(self,other):
        if isinstance(other,(int,float,complex,SUSY,addSUSY)):
            if len((self-other).arr)==0:
                return True
            else:
                return False
        else:
            return False

    def __add__(self,other):
        if isinstance(other,(int,float,complex)):
            return self+SUSY(other,[])
        elif isinstance(other,SUSY):
            return addSUSY([self,other])
        elif isinstance(other,addSUSY):
            return addSUSY([self]+other.arr)
        elif isinstance(other,(Diff,addDiff)):
            return Operator(multi=self,diff=other)
        elif isinstance(other,Operator):
            return other+self
    
    def __radd__(self,other):
        if isinstance(other,(int,float,complex)):
            return SUSY(other,[])+self

    def __sub__(self,other):
        if isinstance(other,(int,float,complex,SUSY,addSUSY,Diff,addDiff,Operator)):
            return self+(-1)*other
    
    def __rsub__(self,other):
        if isinstance(other,(int,float,complex)):
            return other+(-1)*self

    def __mul__(self,other):
        if isinstance(other,SUSY):
            return SUSY(self.number*other.number,self.etalist+other.etalist)
        elif isinstance(other,(int,float,complex)):
            return SUSY(self.number*other,self.etalist)
        elif isinstance(other,addSUSY):
            new_arr=[]
            for i in other.arr:
                new_arr.append(self*i)
            return (addSUSY(new_arr))

    def __rmul__(self,other):
        if isinstance(other,(int,float,complex)):
            return self*other

    def __truediv__(self,other):
        if isinstance(other,(int,float,complex)):
            return self*(1/other)
    
    def __iadd__(self,other):
        if isinstance(other,(int,float,complex,SUSY,addSUSY,Diff,addDiff,Operator)):
            return self+other


    def copy(self):
        return SUSY(self.number,self.etalist)



class addSUSY():
    def __init__(self,arr):

        self.arr=[]

        self.alleta=[]
        for i in range(len(arr)):
            if (arr[i].etalist in self.alleta)==False:
                self.alleta.append(arr[i].etalist)


        for i in self.alleta:
            num_temp=0
            for j in arr:
                if j.etalist==i:
                    num_temp+=j.number
            if num_temp!=0:
                self.arr.append(SUSY(num_temp,i))
        
        self.size=len(self.arr)

    def __repr__(self):
        result="addSUSY([\n"
        for i in range(len(self.arr)):
            if i<len(self.arr)-1:
                result+=(str(self.arr[i])+",\n")
            else:
                result+=(str(self.arr[i])+"\n")
        return result+"])"
    def __str__(self):
        return self.__repr__()

    def get(self,arr):
        thereis=False
        Eta_temp=Eta(arr)
        for i in self.arr:
            if i.etalist==Eta_temp.etalist:
                thereis=True
                return i.number
        if thereis==False:
            return 0.+0.j

    def get_order(self,num):
        result=SUSY_zeros()
        for i in self.arr:
            if i.order==num:
                result+=i
        return result

    

    def contain(self,etalist=[]):
        result=[]
        for susy in self.arr:
            if all (eta in susy.etalist for eta in etalist):
                result.append(susy)
        return addSUSY(result)
    

    def get_count(self,counter):
        result=[]
        for susy in self.arr:
            boolean=True
            for terms in counter:
                list_temp=susy.etalist+terms[0]
                if (len(list_temp)-len(set(list_temp)))!=terms[1]:
                    boolean=False
                    break
            if boolean:
                result.append(susy)
        return addSUSY(result)


    def estimate(self,precision=validity):
        result=[]
        for element in self.arr:
            if abs(complex(element.number))>=precision:
                result.append(SUSY(element.number,element.etalist))
        return addSUSY(result)

    def vanish(self,precision=validity):
        result=bool(True)
        for i in self.arr:
            result=result and i.vanish(precision)
        return (result)
    
    def copy(self):
        return addSUSY(self.arr)


    def __neg__(self):
        return (-1)*self

    def __eq__(self,other):
        if isinstance(other,(int,float,complex,SUSY,addSUSY)):
            if len((self-other).arr)==0:
                return True
            else:
                return False
        else:
            return False
        

    def __add__(self,other):
        if isinstance(other,SUSY):
            return addSUSY(self.arr+[other])
        elif isinstance(other,addSUSY):
            return addSUSY(self.arr+other.arr)
        elif isinstance(other,(int,float,complex)):
            return self+SUSY(other,[])
        elif isinstance(other,(Diff,addDiff)):
            return Operator(multi=self,diff=other)
        elif isinstance(other,Operator):
            return other+self

    def __radd__(self,other):
        if isinstance(other,(int,float,complex)):
            return self+other

    def __sub__(self,other):
        if isinstance(other,(int,float,complex,SUSY,addSUSY,Diff,addDiff,Operator)):
            return self+(-1)*other
    
    def __rsub__(self,other):
        if isinstance(other,(int,float,complex)):
            return other+(-1)*self

    def __mul__(self,other):
        if isinstance(other,(int,float,complex)):
            new_arr=[]
            for i in self.arr:
                new_arr.append(other*i)
            return (addSUSY(new_arr))
        
        elif isinstance(other,SUSY):
            new_arr=[]
            for i in self.arr:
                new_arr.append(i*other)
            return (addSUSY(new_arr))
        
        elif isinstance(other,addSUSY):
            new_arr=[]
            for i in self.arr:
                for j in other.arr:
                    new_arr.append(i*j)
            return (addSUSY(new_arr))

    def __rmul__(self,other):
        if isinstance(other,(int,float,complex)):
            return (self*other)
    
    def __truediv__(self,other):
        if isinstance(other,(int,float,complex)):
            return self*(1/other)

    def __iadd__(self,other):
        if isinstance(other,(int,float,complex,SUSY,addSUSY,Diff,addDiff,Operator)):
            return self+other




'---------------------- Grassmannian operators ----------------------'


def diff(susy,eta):
    if isinstance(susy,(int,float,complex)):
        return SUSY(0.+0.j,[])
    elif isinstance(susy,SUSY):
        etastr_temp=susy.etalist.copy()
        if eta in etastr_temp:
            for i in range(len(etastr_temp)):
                if etastr_temp[i]==eta:
                    index=i
                    break
            del etastr_temp[index]
            return SUSY(susy.number*((-1)**index),etastr_temp)
        else:
            return SUSY(0.+0.j,[])
    elif isinstance(susy,addSUSY):
        result=[]
        for i in susy.arr:
            susy_temp=diff(i,eta)
            if i.number!=0:
                result.append(susy_temp)
        return addSUSY(result)
    else:
        print ("error")


class Diff():
    def __init__(self,number,eta):
        if (isinstance(number,(int,float,complex)))==False:
            print ("error")
        elif eta not in itos.values():
            print ("error")
        self.number=number
        self.differentials=eta
    
    def __repr__(self):
        return "Diff("+str(self.number)+",'"+self.differentials+"')"
    def __str__(self):
        return self.__repr__()

    def diff(self,other):
        return self.number*diff(other,self.differentials)

    def copy(self):
        return Diff(self.number,self.differentials)

    def __neg__(self):
        return self*(-1)

    def __eq__(self,other):
        if (isinstance(other,(Diff,addDiff))):
            SUSY_temp=SUSY(1.,itos.values())
            return self.diff(SUSY_temp)==other.diff(SUSY_temp)

    def __add__(self,other):
        if isinstance(other,Diff):
            return addDiff([self,other])
        elif isinstance(other,addDiff):
            return other+self
        elif (isinstance(other,(int,float,complex,SUSY,addSUSY))):
            return Operator(multi=other,diff=self)
        elif isinstance(other,Operator):
            return other+self
    
    def __radd__(self,other):
        if isinstance(other,(int,float,complex)):
            return self+other
    
    def __sub__(self,other):
        if isinstance(other,(int,float,complex,SUSY,addSUSY,Diff,addDiff,Operator)):
            return self+(-1)*other
    
    def __rsub__(self,other):
        if isinstance(other,(int,float,complex)):
            return other+(-1)*self

    def __mul__(self,other):
        if isinstance(other,(int,float,complex)):
            return Diff(self.number*other,self.differentials)

    def __rmul__(self,other):
        if (isinstance(other,(int,float,complex))):
            return Diff(self.number*other,self.differentials)
    
    def __truediv__(self,other):
        if isinstance(other,(int,float,complex)):
            return self*(1/other)

    def __iadd__(self,other):
        if isinstance(other,(int,float,complex,SUSY,addSUSY,Diff,addDiff,Operator)):
            return self+other


class addDiff():
    def __init__(self,arr):

        self.alldiff=[]
        for i in arr:
            if (i.differentials not in self.alldiff):
                self.alldiff.append(i.differentials)

        self.arr=[]
        for differentials in self.alldiff:
            temp_num=0.
            for i in arr:
                if i.differentials==differentials:
                    temp_num+=i.number
            if temp_num!=0.:
                self.arr.append(Diff(temp_num,differentials))

        
    def __repr__(self):
        result="addDiff([\n"
        for i in range(len(self.arr)):
            result+=(str(self.arr[i])+",\n")
        return result+"])"
    def __str__(self):
        return self.__repr__()
    

    def diff(self,other):
        result=SUSY_zeros()
        for Diffs in self.arr:
            result+=Diffs.diff(other)
        return result
    
    def estimate(self,precision=validity):
        result=[]
        for element in self.arr:
            if abs(complex(element.number))>=precision:
                result.append(Diff(element.number,element.differentials))
        return addDiff(result)
    
    def copy(self):
        return addDiff(self.arr)
    
    def __neg__(self):
        return self*(-1)

    def __eq__(self,other):
        if (isinstance(other,(Diff,addDiff))):
            SUSY_temp=SUSY(1.,itos.values())
            return self.diff(SUSY_temp)==other.diff(SUSY_temp)

    def __add__(self,other):
        if isinstance(other,Diff):
            return addDiff(self.arr+[other])
        elif isinstance(other,addDiff):
            return addDiff(self.arr+other.arr)
        elif (isinstance(other,(int,float,complex,SUSY,addSUSY))):
            return Operator(multi=other,diff=self)
        elif isinstance(other,Operator):
            return other+self
    
    def __radd__(self,other):
        if isinstance(other,(int,float,complex)):
            return self+other
    
    def __sub__(self,other):
        if isinstance(other,(int,float,complex,SUSY,addSUSY,Diff,addDiff,Operator)):
            return self+(-1)*other
    
    def __rsub__(self,other):
        if isinstance(other,(int,float,complex)):
            return other+(-1)*self

    def __mul__(self,other):
        if isinstance(other,(int,float,complex)):
            result=[]
            for i in self.arr:
                result.append(i*other)
            return addDiff(result)
    
    def __rmul__(self,other):
        if isinstance(other,(int,float,complex)):
            return self*other

    def __truediv__(self,other):
        if isinstance(other,(int,float,complex)):
            return self*(1/other)
    
    def __iadd__(self,other):
        if isinstance(other,(int,float,complex,SUSY,addSUSY,Diff,addDiff,Operator)):
            return self+other



def SUSY_zeros():
    return SUSY(0.+0.j,[])

def Diff_zeros():
    return Diff(0.,list(itos.values())[0])


def uSUSY(eta):
    return SUSY(1.,[eta])

def uDiff(eta):
    return Diff(1.,eta)


'---------------------- Operator ----------------------'

class Operator():
    def __init__(self,multi=SUSY_zeros(),diff=Diff_zeros()):
        if isinstance(multi,addSUSY):
            self.multi=multi
        elif isinstance(multi,SUSY):
            self.multi=addSUSY([multi])
        elif (isinstance(multi,(int,float,complex))):
            self.multi=addSUSY([SUSY(multi,[])])
        
        if isinstance(diff,addDiff):
            self.diff=diff
        elif isinstance(diff,Diff):
            self.diff=addDiff([diff])
        
    
    def __repr__(self):
        return "Operator(\nmulti="+str(self.multi)+",\ndiff="+str(self.diff)+")"
    def __str__(self):
        return self.__repr__()
    

    def operate(self,other):
        if isinstance(other,(int,float,complex,SUSY,addSUSY)):
            return self.multi*other+self.diff.diff(other)

    def estimate(self,precision=validity):
        return Operator(multi=self.multi.estimate(precision),diff=self.diff.estimate(precision))

    def __neg__(self):
        return self*(-1)
    
    def __eq__(self,other):
        if isinstance(other,Operator):
            return ((self.multi==other.multi) and (self.diff==other.diff))
        elif isinstance(other,(int,float,complex,SUSY,addSUSY)):
            return ((self.multi==other) and (self.diff==Diff_zeros()))
        elif isinstance(other,(Diff,addDiff)):
            return ((self.multi==SUSY_zeros()) and (self.diff==other))
    
    def __add__(self,other):
        if isinstance(other,(int,float,complex,SUSY,addSUSY)):
            return Operator(multi=self.multi+other,diff=self.diff)
        elif isinstance(other,(Diff,addDiff)):
            return Operator(multi=self.multi,diff=self.diff+other)
        elif isinstance(other,Operator):
            return Operator(multi=self.multi+other.multi,diff=self.diff+other.diff)
    
    def __radd__(self,other):
        if isinstance(other,(int,float,complex)):
            return self+other
    
    def __sub__(self,other):
        if isinstance(other,(int,float,complex,SUSY,addSUSY,Diff,addDiff,Operator)):
            return self+(-1)*other
    
    def __rsub__(self,other):
        if isinstance(other,(int,float,complex)):
            return other+(-1)*self
    
    def __mul__(self,other):
        if isinstance(other,(int,float,complex)):
            return Operator(multi=self.multi*other,diff=self.diff*other)
    
    def __rmul__(self,other):
        if isinstance(other,(int,float,complex)):
            return self*other
    
    def __truediv__(self,other):
        if isinstance(other,(int,float,complex)):
            return self*(1/other)
    
    def __iadd__(self,other):
        if isinstance(other,(int,float,complex,SUSY,addSUSY,Diff,addDiff,Operator)):
            return self+other


def Operator_zeros():
    return Operator(SUSY_zeros(),Diff_zeros())




