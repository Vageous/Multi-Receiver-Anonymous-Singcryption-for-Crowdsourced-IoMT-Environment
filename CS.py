import itertools
import numpy as np
from pypbc import *
from KGI import KGI
from copy import deepcopy
class CS(object):
    def __init__(self,pairing,proof_i,C,cm_i,TS,Att_i_leak,mpk,t,S) -> None:
        self.pairing=pairing
        # [c_i,u_id,u_s,u_k,u_j]
        # [0,1,2,3,4]
        self.proof_i=proof_i
        # [C_i_0,C_i_1,C_i_2,C_i_3,C_i_4,C_i_5]
        # [G2,G1,G1,G1,G1,G1]
        self.C=C
        self.cm_i=cm_i
        self.TS=TS
        self.mpk=mpk
        self.t=t
        self.S=S
        self.Att_i_leak=Att_i_leak
        self.theta=len(Att_i_leak)


    def ZKPoK2_ver(self):
        temp=Element(self.pairing,G1,self.C[3]**self.proof_i[0])
        mpk_g_i=self.mpk[2:self.t+1]
        V=Element.one(self.pairing,G1)
        for i in range(self.theta):
            V=Element(self.pairing,G1,V*Element(self.pairing,G1,mpk_g_i[i]**self.Att_i_leak[i]))
        temp1=Element(self.pairing,G1,V**self.proof_i[3])
        temp2=Element.one(self.pairing,G1)
        for i in range(self.t-self.theta):
            temp2 = Element(self.pairing,G1,temp2*Element(self.pairing,G1,self.mpk[2+i]**self.proof_i[4][i]))
        temp3=Element(self.pairing,G1,self.mpk[1]**self.proof_i[2])
        temp4=Element(self.pairing,G1,self.C[2]**(-self.proof_i[1]))
        new_temp=Element(self.pairing,G1,temp1*temp2*temp3*temp4)
        new_R_i_0=Element(self.pairing,G1,temp*Element.__invert__(new_temp))
        temp5=Element(self.pairing,G1,self.C[4]**self.proof_i[0])
        temp6=Element(self.pairing,G1,self.mpk[len(self.mpk)-2]**(-self.proof_i[3]))
        new_R_i_1=Element(self.pairing,G1,temp5*Element.__invert__(temp6))
        temp7=Element(self.pairing,G1,self.C[5]**self.proof_i[0])
        temp8=Element(self.pairing,G1,(self.mpk[2]**self.proof_i[1])*(self.mpk[len(self.mpk)-1]**(-self.proof_i[3])))
        new_R_i_2=Element(self.pairing,G1,temp7*Element.__invert__(temp8))
        information="".join(list(map(str,self.C)))+str(new_R_i_0)+str(new_R_i_1)+str(new_R_i_2)+str(self.cm_i)+str(self.TS)+"".join(list(map(str,self.Att_i_leak)))
        new_c_i=KGI.H(self.pairing,information)
        if new_c_i==self.proof_i[0]:
            return True
        else:
            return False
    

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
    # {id_1,...,id_s}
    # p_1_s: {id_2,id_3,...,id_s}的排列组合
    # p_1_s=r^(s-2)+r^(s-3)*(C_(s-1)^1)+r^(s-4)*(C_(s-1)^2)+...+(C_(s-1)^(s-2)))
    def p_v_s_compute(self,v):
       id_s=deepcopy(self.S)
       id_s.pop(v)
       mpk_h=self.mpk[self.t+2:]
       mpk_h.pop(1)
       id_s_com=self.list_com_mul(id_s)
       id_s_com.pop(len(id_s_com)-1)
       p_v_s=Element.one(self.pairing,G2)
       for i in range(len(id_s_com)):
            p_v_s=Element(self.pairing,G2,p_v_s*Element(self.pairing,G2,mpk_h[i]**id_s_com[len(id_s_com)-1-i]))
       return p_v_s
    
    # v denote v-th professional
    def Unsigncrypt(self,v):
        e1=self.pairing.apply(self.C[3],self.mpk[self.t+2])
        e2=self.pairing.apply(self.C[2],self.mpk[self.t+4])
        if e1==e2:
            if self.ZKPoK2_ver():
                p_v_s=self.p_v_s_compute(v)
                return p_v_s
            else:
                print("The second zero knowledge verification cannot pass!")
        else:
            print("The billinear maps verification cannot pass!")
        
    
