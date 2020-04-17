import os

directory = r''

filecounter = 0
with open("parsed.dat","w+") as out:
    for filename in os.listdir(directory):
        filecounter+=1
        if filename.endswith(".dat"):
            print(filename)
            with open(directory+'/'+filename,'r',buffering=1000000) as inp:
                i = 0
                switch = 0
                conf = []
                for line in inp:               
                    if ('t' in line):
                        if int(line[4:-1])%2000000==0:
                            switch = 1
                            line = line[0:4]+'%s'%filecounter+'%s'%filecounter \
                            +'%s'%filecounter+line[4:-1]+line[-1]
                        elif switch == 0 :
                            pass
                        else :
                            switch = 0
                            out.writelines("%s" % lin for lin in conf)
                            print(conf[0])
                            conf = []
                            i += 1
                    
                    if switch == 1:
                        conf.append(line)
                    
                    if i == 100:
                        break
