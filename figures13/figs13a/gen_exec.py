import os

listdir = os.listdir('.')
f = open('get_ft2.sh', 'w')
for name in listdir:
    if name[:2] == 'at':
        f.write('cd ' + name + '/2/' + '\n')
        f.write('fid.com\n')
        f.write('nmrproc_1d.com\n')
        f.write('cd ../../' + '\n\n')
f.close() 
