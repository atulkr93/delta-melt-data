import os

listdir = os.listdir('.')
f = open('get_ft2.sh', 'w')
for name in listdir:
    if name[:4] == 'atul':
        f.write('cd ' + name + '/4/' + '\n')
        f.write('fid.com\n')
        f.write('nmrproc_1d.com\n')
        f.write('cd ../../' + '\n\n')
f.close() 
