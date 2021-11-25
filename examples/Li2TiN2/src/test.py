import sys
import os
import shutil
import yaml
from module_input import *
import numpy as np
import subprocess

# input
a = []
a.append([1,2,3])
a.append([2,3,4])
a.append([3,4,5])

a = np.array(a)

try:
    c = subprocess.check_output(['grep',"'free  '",'main.py'])
except:
    c = "1"

c = os.listdir('./')
print (c)
  
num_sample = 1
dat_sorted=[1,2,3,4]

for j in range(10,1,-1):
  print (j)

dat = []
print (len(dat))


a = 0
for i in range(1000):
    if a == 0:
        print (i)
        if i == 10:
            break


dat = [[1,2], [62,3], [78,34],[234,32]]

dat = sorted(dat, key = lambda x : x[1])[0:2]

print (dat)

for i in range(10,1,-1):
    print(i)

with open("dat","r") as f:
    d = []
    for fstr in f.readlines():
        d.append(fstr.split()[0])

print (d)


for i in range(len(d)-1,-1,-1):
   if i ==  3 or i == 5 or i == 7:
       del d[i]

print (d)
