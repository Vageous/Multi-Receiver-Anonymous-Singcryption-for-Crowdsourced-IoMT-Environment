from option import parser_args
def compute_time(time,i):
    args=parser_args()
    for key,value in time.items():
        with open("./compute_time/{}.txt".format(key),'a+') as f:
            if i==0:
                f.write("t:{},n:{},s:{},theta:{}\n".format(args.t,args.n,args.num,args.theta))
                f.write("{}\n".format(value))
            else:
                f.write("{}\n".format(value))