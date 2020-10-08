import os
from shutil import copyfile

# Create Omega file ranging from 20~50 rpm (20,30,40,50)
l = 1
for i in range(20, 60, 10):
    copyfile('omega.txt', ('omega'+str(l)+'.txt'))
    l=l+1
    fin = open("omega.txt", "rt")
    data = fin.read()
    data = data.replace(' '+str(i),' '+str(i+10))
    fin.close()
    fin = open("omega.txt", "wt")
    fin.write(data)
    fin.close()

# resets omega.txt ( original omega file to initial 20 rpm value)
fin = open("omega.txt", "rt")
data = fin.read()
data = data.replace(' '+str(i+10),' '+str(20))
fin.close()
fin = open("omega.txt", "wt")
fin.write(data)
fin.close()

# copy option file with changes in pitch/length and omegainput file.
helixpitch = [0.019, 0.032, 0.045, 0.057, 0.070, 0.019]
axisLengthInput = [0.032, 0.064, 0.095, 0.127, 0.159, 0.032]
omegainput = [1,2,3,4,1]
count = 1 

for j in range(0, 5, 1):
    for k in range (0,4,1):
        fin = open('option.txt','rt')
        data = fin.read()
        data = data.replace(
            ('helixpitch '+str(helixpitch[j])), ('helixpitch '+str(helixpitch[j+1])))
        data = data.replace(
            ('axisLengthInput '+str(axisLengthInput[j])), ('axisLengthInput '+str(axisLengthInput[j+1])))
        data = data.replace(('input-file omega'+str(omegainput[k])+'.txt'),
                            ('input-file omega'+str(omegainput[k+1])+'.txt'))
        fin.close()
        fin = open('option.txt', "wt")
        fin.write(data)
        fin.close()
        copyfile('option.txt', ('option'+str(count)+'.txt'))
        count = count+1

