import numpy as np
#

f = open('original_data/Grouplevel.txt', 'r')
output = open('Grouplevel2.txt', "w")
for i in range(0,12):
    a=f.readline()
    output.write('# ' + a)

while 1:
    a=f.readline()
    if a=='': break
    a=a.split()
    a.append('\n')
    if a[0] != "#":
        a[3] = str(-float(a[3]))
    if a[7] == 'G':
        a[4] = str(0.5 * float(a[4]) )
        a[7] = 'N'
    output.writelines( " ".join(a) )

f.close()
output.close()



f = open('original_data/Mixed_Strat.txt', 'r')
output = open('Mixed_Strat2.txt', "w")
for i in range(0,12):
    a=f.readline()
    output.write('# ' + a)

while 1:
    a=f.readline()
    if a=='': break
    a=a.split()
    a.append('\n')
    # convert BC age to date
    a[3] = str(-float(a[3]))
   # a[6] = str(max(2,5*float(a[6])))
    if a[7] == 'G':
        a[4] = str(0.5 * float(a[4]) )
        a[7] = 'N'

    if a[0] != "#":
        output.writelines( " ".join(a) )

f.close()
output.close()

