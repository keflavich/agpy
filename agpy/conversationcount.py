import glob
import subprocess
import os
import gzip
import pwd

path = '/Users/%s/Library/Application Support/Adium 2.0/Users/Default/Logs/' % (pwd.getpwuid( os.getuid() )[ 0 ])

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

if __name__ == "__main__":
    from pylab import *
    figure(1)
    N = array(nmsgdict.values())
    hist(log10(N[N>0]))
    xlabel("log Number of messages")
    ylabel("Number of users")

    figure(2)
    C = array(countdict.values())
    hist(log10(C[C>0]))
    xlabel("log Number of conversations")
    ylabel("Number of users")
    
