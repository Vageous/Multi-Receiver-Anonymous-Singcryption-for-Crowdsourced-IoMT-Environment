import argparse

def parser_args():
    parser=argparse.ArgumentParser()

    parser.add_argument("--cure_param",type=str,default='''type f
q 205523667896953300194896352429254920972540065223
r 205523667896953300194895899082072403858390252929
b 40218105156867728698573668525883168222119515413
beta 115334401956802802075595682801335644058796914268
alpha0 191079354656274778837764015557338301375963168470
alpha1 71445317903696340296199556072836940741717506375
''',help="The parameters for the type f curve")
    parser.add_argument("--t",type=int,default=5,help="Attribute number")
    parser.add_argument("--n",type=int,default=10,help="max receiver number")
    parser.add_argument("--num",type=int,default=5,help="receiver number")
    parser.add_argument("--theta",type=int,default=3,help="Attribute leakage number")


    args=parser.parse_args()
    return args