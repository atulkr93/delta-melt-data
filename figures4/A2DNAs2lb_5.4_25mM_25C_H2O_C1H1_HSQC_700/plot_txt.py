import numpy as np
import matplotlib.pyplot as plt

x = []
y = []
z = []
x_master = []
y_master = []
z_master = []
yval = 1

f = open('see.txt', 'r')
for line in f.readlines():
    line = line.strip('\n').split(' ')
    line = [ele for ele in line if ele != '']
    x.append(int(line[0]))
    y.append(int(line[1]))
    z.append(float(line[2]))
    if len(x) == 821:
        x_master.append(x)
        y_master.append(y)
        z_master.append(z)
        x = []
        y = []
        z = []
f.close()
print "done"
print len(x_master)

x_master = np.array(x_master)
y_master = np.array(y_master)
z_master = np.array(z_master)
print x_master.shape, y_master.shape, z_master.shape

plt.figure(1)
plt.contour(y_master, x_master, z_master)
plt.show()
