import numpy as np
from pypbc import *
from KGI import KGI
from copy import deepcopy

class Professional(object):
    def __init__(self,pairing,mpk,t,S,Att_i_leak) -> None:
        self.pairing=pairing
        self.mpk=mpk
        self.t=t
        self.S=S
        self.Att_i_leak=Att_i_leak

    def D_v_i_ver(self,D_v_i,id_v_i):
        temp1=Element(self.pairing,G2, self.mpk[self.t+4]*(Element(self.pairing,G2, self.mpk[self.t+2]**id_v_i)))
        e1=self.pairing.apply(D_v_i,temp1)
        e2=self.mpk[self.t+3]
        if e1==e2:
            return True
        else:
            return False
    
    def D_v_ver(self,D_v):
        ver_result=[]
        for i in range(len(self.S)):
            ver_result.append(self.D_v_i_ver(D_v[i],self.S[i]))
        
        if sum(ver_result)==len(self.S):
            return True
        else:
            return False

    def str_XOR(self,str1,str2):
        return "".join([str(ord(a)^ord(b)) for a,b in zip(str1,str2)])

    def decrypt(self,C_i,p_v_s,D_v,v,cm_i,TS):
        E_v=self.pairing.apply(C_i[4],p_v_s)
        F_v=self.pairing.apply(D_v,C_i[1])
        id_s=deepcopy(self.S)
        id_s.pop(v)
        item=Element(self.pairing,GT,E_v*F_v)
        T_v=Element(self.pairing,GT,item**(Element.__invert__(Element(self.pairing,Zr,int(np.prod(id_s))))))
        new_m_r=self.str_XOR(C_i[0],KGI.H_1(T_v))
        new_m=new_m_r[0:63]
        new_r=new_m_r[64:]
        temp=KGI.H(self.pairing,new_r+"".join(list(map(str,self.Att_i_leak)))+TS)
        if cm_i==temp:
            return True
        else:
            return False

        
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
#     id_v=['professional 1','professional 2','professional 3']
#     test1=Professional(test.pairing,id_v,test.mpk,t)

#     D_v=test.D_v_gen(id_v)
#     test1.D_v_ver(D_v)
    

    
            
            
