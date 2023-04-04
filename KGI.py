import sys
import random
import hashlib
from pypbc import *


class KGI(object):
    def __init__(self,stored_param,t,n) -> None:
        self.t=t
        self.n=n
        self.msk,self.mpk=[],[]
        self.stored_param=stored_param
        self.param=Parameters(param_string=self.stored_param)
        self.pairing = Pairing(self.param)
        

    # msk: (lsit) [r,fai]
    def msk_gen(self):
        for i in range(2):
            self.msk.append(Element.random(self.pairing,Zr))

    # mpk: (list) [g,g',g1,...,gt,h,h_1,....,hn,T,Y] -> [0,1,2,...,t+1,t+2,t+3,...,t+n+2,t+n+3,t+n+4]
    def mpk_gen(self):
        # g,g',g1,...,gt
        for i in range(self.t+2):
            self.mpk.append(Element.random(self.pairing,G1))
        # h
        self.mpk.append(Element.random(self.pairing,G2))
        # e(g,h)
        self.mpk.append(self.pairing.apply(self.mpk[0],self.mpk[len(self.mpk)-1]))
        # h1,h2,....,hn
        for i in range(1,self.n+1):
            self.mpk.append(Element(self.pairing,G2,self.mpk[self.t+2]**(self.msk[0]**i)))
        # T 
        T=Element(self.pairing,G1,self.mpk[0]**self.msk[0])
        self.mpk.append(Element(self.pairing,G1,self.mpk[0]**self.msk[0]))
        # Y
        Y=Element(self.pairing,G1,T**self.msk[1])
        self.mpk.append(Y)
        
    # {0,1}^* -> Zp
    @classmethod
    def H(self,pairing,message):
        return Element.from_hash(pairing, Zr, message)

    # {0,1}^* -> Zp
    @classmethod
    def H_0(self,pairing,message):
        return Element.from_hash(pairing,Zr,message)

    # GT -> {0,1}^2l l=128bits
    # return 2进制
    @classmethod
    def H_1(self,GT):
        obj = hashlib.md5()
        obj.update(str(GT).encode('utf-8'))
        result = obj.hexdigest()
        result_0b=bin(int(result,16))[2:].zfill(len(result)*4)
        return result_0b

    def ZKPok1_ver(self,proof_i):
        temp1=Element(self.pairing,G1,proof_i[2]**proof_i[3])
        temp2=Element(self.pairing,G1,self.mpk[1]**proof_i[4])
        new_R_i=Element(self.pairing,G1,temp1*Element.__invert__(temp2))
        information=proof_i[0]+str(proof_i[2])+str(new_R_i)
        for i in range(len(proof_i[1])):
            information += str(proof_i[1][i])
        new_c_i=self.H(self.pairing,information)
        if new_c_i==proof_i[3]:
            return True
        else:
            return False

    def Attri_compute(self,proof_i):
        temp=Element.one(self.pairing,G1)
        Reg=[]
        id_i=proof_i[0]
        Att_i=proof_i[1]
        E_i=proof_i[2]
        for i in range(len(Att_i)):
            temp= Element(self.pairing,G1,temp* Element(self.pairing,G1,self.mpk[i+2]**Att_i[i]))
        temp1=Element(self.pairing,G1,E_i*temp)
        A_i=Element(self.pairing,G1,temp1**(Element.__invert__((self.msk[0]+self.H_0(self.pairing,id_i)))))
        Reg.append(id_i)
        Reg.append(Element(self.pairing,G1,self.mpk[2]**self.H_0(self.pairing,id_i)))
        Reg.append(A_i)
        return Reg
        
    def D_v_gen(self,id_v):
        D_v=[]
        for i in range(len(id_v)):
            temp=Element(self.pairing,G1,self.mpk[0]**Element.__invert__((self.msk[0]+Element(self.pairing,Zr,id_v[i]))))
            D_v.append(temp)
        return D_v
    

    def Trace(self,C,Reg):
        temp=Element(self.pairing,G1,C[4]**self.msk[1])
        flag=Element(self.pairing,G1,C[5]*Element.__invert__(temp))
        if flag==Reg[1]:
            return Reg[0]
        else:
            return False

# if __name__=='__main__':
#     stored_params='''type f
# q 205523667896953300194896352429254920972540065223
# r 205523667896953300194895899082072403858390252929
# b 40218105156867728698573668525883168222119515413
# beta 115334401956802802075595682801335644058796914268
# alpha0 191079354656274778837764015557338301375963168470
# alpha1 71445317903696340296199556072836940741717506375
# '''
#     test=KGI(stored_params,5,10)
#     test.msk_gen()
#     print("The maste private key list\n")
#     print(test.msk)
#     test.mpk_gen()
#     print("The master public key list\n")
#     print(test.mpk)

#     # Hash test
#     message="test"
#     print("The hash function H\n")
#     print(test.H(message))
#     print("The hash function H_0\n")
#     print(test.H_0(message))
#     print("The hash function H_1\n")
#     data=Element(test.pairing,GT)
#     print(test.H_1(data))
#     print(sys.getsizeof(test.H_1(data)))
