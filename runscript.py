# This python script should run inside the folder "~~/flagellarobot"
# automatically increase omega value from 30~60. needs to start from 30 for now.
# Returns the omega value back to 30rpm after the script finishes running

import os
for i in range(30,70,10):
    os.system('./simDER option.txt')
    fin = open("omega.txt","rt")
    data = fin.read()
    data = data.replace(str(i),str(i+10))
    fin.close()
    fin = open("omega.txt","wt")
    fin.write(data)
    fin.close()
fin = open("omega.txt", "rt")
data = fin.read()
data = data.replace(str(i+10), str(30))
fin.close()
fin = open("omega.txt", "wt")
fin.write(data)
fin.close()
