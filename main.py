import time
import random
from CS import CS
from pypbc import *
from KGI import KGI
from Worker import Worker
from option import parser_args
from compute_time import compute_time
from Professional import Professional

def main(stored_params,t,n,num,theta,index):

    # KGI classs
    time13=time.time()
    kgi=KGI(stored_params,t,n)
    # the master keys gen
    kgi.msk_gen()
    kgi.mpk_gen()
    time14=time.time()-time13
    # worker's id
    id_i='worker'
    # worker's attribute list
    Att_i=[]
    for i in range(t):
        Att_i.append(Element.random(kgi.pairing,Zr))
    # worker's attribute leakage
    # the number of attribute leakage
    Att_i_leak=random.sample(Att_i,theta)
    # professional number
    # the receiver list
    S=[]
    for i in range(1,num+1):
        S.append(i)

    # worker class
    time11=time.time()
    worker=Worker(kgi.pairing,id_i,Att_i,Att_i_leak,S,kgi.mpk)
    # zero knowledge generation
    proof_1=worker.ZKPoK1_gen()
    # zero knowledge verification
    if kgi.ZKPok1_ver(proof_1):
        print("The first zero knowledge verification pass!\n")
        # Reg generation
        Reg=kgi.Attri_compute(proof_1)
        cred_i=worker.Attri_ver(Reg[2])
    else:
        print("The first zero knowledge verification cannot pass!\n")
    time12=time.time()-time11
    # professional class
    professional=Professional(kgi.pairing,kgi.mpk,t,S,Att_i_leak)
    # all professional in S private key gen
    time9=time.time()
    D_v_list=kgi.D_v_gen(S)
    # all professionals' private key verification
    if professional.D_v_ver(D_v_list):
        print("All professional's private key generated by KGI are valid!\n")
    else:
        print("All professional's private key generated by KGI are invalid!\n")
    time10=time.time()-time9
    # worker's message m
    m=1000
    # worker signcryption
    # C=[C_i_0,C_i_1,C_i_2,C_i_3,C_i_4,C_i_5],cm_i,TS
    time1=time.time()
    C,cm_i,TS=worker.signcrypt(m,cred_i)
    time2=time.time()-time1
    # zero knowledge generation
    proof_2=worker.ZKPoK2_gen(C,cm_i,TS)

    # CS class
    cs=CS(kgi.pairing,proof_2,C,cm_i,TS,Att_i_leak,kgi.mpk,t,S)
    # some information verfication and p_v_s generation
    p_v_s=[]
    time3=time.time()
    for i in range(len(S)):
        p_v_s.append(cs.Unsigncrypt(i))
    time4=time.time()-time3
    # professional S decryption 
    decrypt_result=[]
    time5=time.time()
    for i in range(len(S)):
        decrypt_result.append(professional.decrypt(C,p_v_s[i],D_v_list[i],i,cm_i,TS))
    if sum(decrypt_result)==len(S):
        print("All professionals' decryption are correct!\n")
    else:
        print("One or Some professionals' decryption is/are incorrect!\n")
    time6=time.time()-time5
    # Trace
    time7=time.time()
    id=kgi.Trace(C,Reg)
    print("The worker's identity is {}\n".format(id))
    time8=time.time()-time7

    compute={"signcrypt":time2,"unsigncrypt":time4,"decrypt":time6,"trace":time8,"dekeygen":time10,"ckeygen":time12,"setup":time14}
    compute_time(compute,index)
 

if __name__=='__main__':
    # stored_params: curve parameter
    # Attribute number: t
    # h_i number: n
    # max receiver number: num
    # Attribute leakage number: theta
    # main(t,n,num,theta)
    args=parser_args()
    stored_params=args.cure_param
    t=args.t
    n=args.n
    num=args.num
    theta=args.theta
    for i in range(10):
        main(stored_params,t,n,num,theta,i)


# python main.py --t=2 --n=2 --num=2 --theta=1



