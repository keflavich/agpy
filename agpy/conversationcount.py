"""
Just for fun, histogram up number of conversations & number of messages
exchanged with each person on your buddy list.  Meant for adium logs.
"""
import glob
import subprocess
import os
import gzip
import pwd

path = '/Users/%s/Library/Application Support/Adium 2.0/Users/Default/Logs/' % (pwd.getpwuid( os.getuid() )[ 0 ])

if __name__ == "__main__":
    username = ""
    countdict = {}
    nmsgdict = {}
    for root,dirs,files in os.walk(path):
        if username =="" or "chatlog" not in root:
            if username != "":
                print "User %s had %i conversations and %i lines" % (username,countdict[username],nmsgdict[username])
            username = os.path.split(root)[1]
            countdict[username] = 0
            nmsgdict[username] = 0
        else:
            countdict[username] += 1
            for filename in files:
                p = subprocess.Popen("grep -c message '%s'" % (root+"/"+filename),shell=True,stdout=subprocess.PIPE)
                nlines = int(p.stdout.readline())
                #nlines = os.system("grep -c message '%s'" % (root+"/"+filename))
                nmsgdict[username] += nlines


    from pylab import *
    figure(1)
    N = array(nmsgdict.values())
    hist(log10(N[N>0]))
    xlabel("log Number of messages")
    ylabel("Number of users")

    NN = array(countdict.keys())
    print "Top 10 most messages: \n"+"\n".join("%s: %i" % (a,b) for a,b in zip(NN[argsort(N)[-10:]],sort(N)[-10:]))

    figure(2)
    C = array(countdict.values())
    hist(log10(C[C>0]))
    xlabel("log Number of conversations")
    ylabel("Number of users")

    CN = array(countdict.keys())
    print "Top 10 most conversations \n"+"\n".join("%s: %i" % (a,b) for a,b in zip(CN[argsort(C)[-10:]],sort(C)[-10:]))

    figure(5) # idea courtesy Jordan Mirocha
    loglog(C,N,'ko')
    xlabel("Number of conversations")
    ylabel("Number of messages exchanged")
    for ii in unique(concatenate([argsort(C)[-11:],argsort(N)[-11:]])):
        text(C[ii],N[ii],CN[ii].split('@')[0])
    title("Adium conversation history")

    try:
        import plfit
        pn = plfit.plfit(N[N>0])
        figure(3)
        pn.plotcdf()

        pc = plfit.plfit(C[C>0])
        figure(4)
        pc.plotcdf()

    except ImportError:
        # if you haven't installed the plfit code from agpy
        pass
