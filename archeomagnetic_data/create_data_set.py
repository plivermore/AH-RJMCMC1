import numpy as np
min_age = 1000.0
max_age = 2000.0
#
f = open('Lubeck.txt', 'r')
a=f.readline()
output = open('Lubeck_Paris700.txt', "w")
output.write('# Lubeck+Paris, age range ' + str(min_age)+ ' - ' + str(max_age) + '   Lat  Long  Age   dt    Int-Paris  sd    rep\n')

while 1:
    a=f.readline()
    if a=='': break
    a=a.split()
    a.append('\n')
    if a[0] != "#":
        if float(a[2]) <= max_age  and float(a[2]) >= min_age:
            output.writelines( " ".join(a) )
f.close()

f = open('Paris700.txt', 'r')
a=f.readline()
while 1:
    a=f.readline()
    if a=='': break
    a=a.split()
    a.append('\n')
    if a[0] != "#":
        if float(a[2]) <= max_age  and float(a[2]) >= min_age:
            output.writelines( " ".join(a) )
f.close()

output.close()

