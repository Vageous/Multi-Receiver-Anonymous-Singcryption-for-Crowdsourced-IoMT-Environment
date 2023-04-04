import time
import random
import logging
import itertools
import numpy as np
from pypbc import *
from KGI import KGI



class Worker(object):
    def __init__(self,pairing,id_i,Att_i,Att_i_leak,S,mpk) -> None:
        self.pairing=pairing
        self.mpk=mpk
        self.S=S
        self.l=64
        self.id_i=id_i
        # Att_i (list): Zp
        self.Att_i=Att_i
        self.s_i=None
        self.E_i=None
        self.k_i=None
        self.Att_i_leak=Att_i_leak
        self.t=len(self.Att_i)
        self.theta=len(self.Att_i_leak)
    
    def ZKPoK1_gen(self):
        proof_i=[]
        proof_i.append(self.id_i)
        proof_i.append(self.Att_i)
        self.s_i=Element.random(self.pairing,Zr)
        r_i=Element.random(self.pairing,Zr)
        self.E_i=Element(self.pairing,G1,self.mpk[1]**self.s_i)
        proof_i.append(self.E_i)
        R_i=Element(self.pairing,G1,self.mpk[1]**r_i)
        information=self.id_i+str(self.E_i)+str(R_i)
        for i in range(len(self.Att_i)):
            information += str(self.Att_i[i])
        c_i=KGI.H(self.pairing,information)
        u_i=Element(self.pairing,Zr,self.s_i*c_i-r_i)
        proof_i.append(c_i)
        proof_i.append(u_i)
        return proof_i
   
    def Attri_ver(self,A_i):
        cred_i=[]
        cred_i.append(A_i)
        temp1=Element(self.pairing,G2, self.mpk[self.t+4]*(Element(self.pairing,G2, self.mpk[self.t+2]**KGI.H_0(self.pairing,self.id_i))))
        temp=Element.one(self.pairing,G1)
        for i in range(len(self.Att_i)):
            temp= Element(self.pairing,G1,temp* Element(self.pairing,G1,self.mpk[i+2]**self.Att_i[i]))
        temp2=Element(self.pairing,G1,self.E_i*temp)
        e1=self.pairing.apply(A_i,temp1)
        e2=self.pairing.apply(temp2,self.mpk[self.t+2])
        if e1==e2:
            B_i=Element(self.pairing,G1,(A_i**(-KGI.H_0(self.pairing,self.id_i))*temp2))
            cred_i.append(B_i)
            cred_i.append(self.s_i)
            return cred_i
        else:
            return False

    def int_to_bits(self,m):
        m_0x=format(m,"b")
        length=len(m_0x)
        if len(m_0x)<self.l:
            for i in range(self.l-len(m_0x)):
                m_0x +='0'
            return m_0x,length
        else:
            print("The length of m is out of the limation")
    

    def str_XOR(self,str1,str2):
        return "".join([str(ord(a)^ord(b)) for a,b in zip(str1,str2)])
    

    def bits_to_int(self,m_0x,length):
        return int(m_0x[0:length])

    def list_mul(self,id_s,num):
        result=[]
        combination=list(itertools.combinations(id_s,num))
        for i in combination:
            result.append(np.prod(list(i)))
        return sum(result)

    def list_com_mul(self,id_s):
        result=[]
        for i in range(len(id_s)+1):
            result.append(Element(self.pairing,Zr,int(self.list_mul(id_s,i))))
        return result

    def signcrypt(self,m,cred_i):
        C=[]
        # 2进制字符串
        m_0x,length=self.int_to_bits(m)
        self.k_i=Element.random(self.pairing,Zr)
        # 2进制字符串
        r_i=self.int_to_bits(random.randint(100,300))[0]
        C_i_0=self.str_XOR(m_0x+r_i,KGI.H_1(Element(self.pairing,GT,self.mpk[self.t+3]**self.k_i)))
        # print(KGI.H_1(Element(self.pairing,GT,self.mpk[self.t+3]**self.k_i)))
        # print("\n")
        C.append(C_i_0)
        ran_Att_i_str="".join(list(map(str,self.Att_i_leak)))
        TS=str(time.time())
        cm_i=KGI.H(self.pairing,r_i+ran_Att_i_str+TS)
        S_combinations=self.list_com_mul(self.S)
        # print(S_combinations)
        # t+2,t+3,...,t+s+2
        mpk_h=self.mpk[self.t+2:]
        mpk_h.pop(1)
        C_i_1=Element.one(self.pairing,G2)
        for i in range(len(S_combinations)):
            C_i_1=Element(self.pairing,G2,C_i_1*Element(self.pairing,G2,mpk_h[i]**(S_combinations[len(S_combinations)-1-i])))
        C.append(Element(self.pairing,G2,C_i_1**self.k_i))
        C_i_2=Element(self.pairing,G1,cred_i[0]**self.k_i)
        C.append(C_i_2)
        C_i_3=Element(self.pairing,G1,cred_i[1]**self.k_i)
        C.append(C_i_3)
        C_i_4=Element(self.pairing,G1,self.mpk[len(self.mpk)-2]**(-self.k_i))
        C.append(C_i_4)
        C_i_5=Element(self.pairing,G1,(self.mpk[2]**KGI.H_0(self.pairing,self.id_i))*(self.mpk[len(self.mpk)-1]**(-self.k_i)))
        C.append(C_i_5)
        return C,cm_i,TS

    def ZKPoK2_gen(self,C,cm_i,TS):
        # r_id,r_S,r_k,r_j1,r_j2,...,r_jt-theta sum=t-theta+3
        random_number,proof_i=[],[]
        for i in range(self.t-self.theta+3):
            random_number.append(Element.random(self.pairing,Zr))
        temp=Element.one(self.pairing,G1)
        for i in range(len(random_number[3:])):
            temp = Element(self.pairing,G1,temp*Element(self.pairing,G1,self.mpk[2+i]**random_number[i+3]))
        mpk_g_i=self.mpk[2:self.t+1]
        Att_i_leak=random.sample(self.Att_i,self.theta)
        Att_i_secret=[]
        for i in self.Att_i:
            if i not in Att_i_leak:
                Att_i_secret.append(i)
        V=Element.one(self.pairing,G1)
        for i in range(self.theta):
            V=Element(self.pairing,G1,V*Element(self.pairing,G1,mpk_g_i[i]**Att_i_leak[i]))
        temp1=Element(self.pairing,G1,C[2]**(-random_number[0]))
        temp2=Element(self.pairing,G1,self.mpk[1]**random_number[1])
        temp3=Element(self.pairing,G1,V**random_number[2])
        R_i_0=Element(self.pairing,G1,temp1*temp2*temp3*temp)
        R_i_1=Element(self.pairing,G1,self.mpk[len(self.mpk)-2]**(-random_number[2]))
        R_i_2=Element(self.pairing,G1,(self.mpk[2]**random_number[0])*(self.mpk[len(self.mpk)-1]**(-random_number[2])))
        information="".join(list(map(str,C)))+str(R_i_0)+str(R_i_1)+str(R_i_2)+str(cm_i)+str(TS)+"".join(list(map(str,Att_i_leak)))
        c_i=KGI.H(self.pairing,information)
        proof_i.append(c_i)
        u_id=Element(self.pairing,Zr,KGI.H_0(self.pairing,self.id_i)*c_i-random_number[0])
        proof_i.append(u_id)
        u_s=Element(self.pairing,Zr,self.k_i*self.s_i*c_i-random_number[1])
        proof_i.append(u_s)
        u_k=Element(self.pairing,Zr,self.k_i*c_i-random_number[2])
        proof_i.append(u_k)
        u_j=[]
        for i in range(len(random_number[3:])):
            u_j.append(Element(self.pairing,Zr,self.k_i*Att_i_secret[i]*c_i-random_number[i+3]))
        proof_i.append(u_j)
        return proof_i
        

# test
# if __name__=='__main__':
#     stored_params='''type f
# q 205523667896953300194896352429254920972540065223
# r 205523667896953300194895899082072403858390252929
# b 40218105156867728698573668525883168222119515413
# beta 115334401956802802075595682801335644058796914268
# alpha0 191079354656274778837764015557338301375963168470
# alpha1 71445317903696340296199556072836940741717506375
# '''
#     t=5
#     n=10
#     test=KGI(stored_params,t,n)
#     test.msk_gen()
#     test.mpk_gen()
#     id_i='user1'
#     Att_i=[]
#     for i in range(t):
#         Att_i.append(Element.random(test.pairing,Zr))
#     test1=Worker(test.pairing,id_i,Att_i,test.mpk)
#     proof_i=test1.ZKPoK1_gen()
#     print("Zero Knowledge Generation\n")
#     print(proof_i)
#     print("Zero Knowledge Verification\n")
#     ver_result=test.ZKPok1_ver(proof_i)
#     print(ver_result)
#     print("Reg Generation\n")
#     Reg=test.Attri_compute(ver_result,proof_i)
#     print(Reg)
    
#     print("Reg Verification\n")
#     cred_i=test1.Attri_ver(Reg[2])
#     print(cred_i)

