'''
Created on Apr 15, 2014

@author: mrakitin
'''

import os, sys

#-------------------------------------------------------------------------------
# Arguments processing:

input_xyz = ''
if len(sys.argv) < 2:
    print 'Tinker input XYZ file is not specified! Exit.'
    exit(1)
else:
    input_xyz = sys.argv[1]

if not os.path.exists(input_xyz):
    print 'File ' + input_xyz + ' does not exist! Exit.'
    exit(2)
#-------------------------------------------------------------------------------


# Read input xyz file from Tinker:
f = open(input_xyz, 'rb')
content = f.readlines()
f.close()

row       = content[0]
atoms_num = row.split()[0]
name      = ' '.join(row.split()[1:])

atoms_dict = {
              'num' : [],
              'atom': [],
              'x'   : [],
              'y'   : [],
              'z'   : [],
              }

for i in range(1, len(content)):
    row = content[i]
    
    row_list = row.split()
    atoms_dict['num'].append(int(row_list[0]))
    atoms_dict['atom'].append(row_list[1][0])
    atoms_dict['x'].append(float(row_list[2]))
    atoms_dict['y'].append(float(row_list[3]))
    atoms_dict['z'].append(float(row_list[4]))


# Find max abs value to make POSCAR coordinates <= 1:
tmp_crds = atoms_dict['x'] + atoms_dict['y'] + atoms_dict['z']
tmp_crds.sort()
min_value = tmp_crds[0]
max_value = tmp_crds[-1]
alat = max_value - min_value


# Find unique atoms:
uniq_atoms = []
for atom in atoms_dict['atom']:
    if atom not in uniq_atoms:
        uniq_atoms.append(atom)


# Find number of unique atoms and their coordinates:
uniq_nums = []
crds      = {
             'x': [],
             'y': [],
             'z': [],
             }

for uniq_atom in uniq_atoms:
    number = 0
    count  = 0
    for atom in atoms_dict['atom']:
        if uniq_atom == atom: 
            number += 1
            crds['x'].append(atoms_dict['x'][count])
            crds['y'].append(atoms_dict['y'][count])
            crds['z'].append(atoms_dict['z'][count])
        count += 1
    uniq_nums.append(str(number))


# Prepare POSCAR contents:
poscar = 'EA ' + name + '''
''' + '%1.6f' % (1) + '''
    ''' + '%5.8f %5.8f %5.8f' % (1*alat, 0,      0     ) + '''
    ''' + '%5.8f %5.8f %5.8f' % (0,      1*alat, 0     ) + '''
    ''' + '%5.8f %5.8f %5.8f' % (0,      0,      1*alat) + '''
    ''' + '  '.join(uniq_atoms) + '''
    ''' + ' '.join(uniq_nums)  + '''
Direct
'''

crds_fmt = ''
for i in range(len(atoms_dict['num'])):
    crds_fmt += '%12.8f %12.8f %12.8f\n' % ((crds['x'][i] - min_value)/alat + abs(min_value)*0.05/alat,\
                                            (crds['y'][i] - min_value)/alat + abs(min_value)*0.05/alat,\
                                            (crds['z'][i] - min_value)/alat + abs(min_value)*0.05/alat)
poscar += crds_fmt

#print poscar


# Write POSCAR:
#outfile = name.replace(' ', '_') + '_POSCAR.vasp'
outfile = 'POSCAR'
f = open(outfile, 'wb')
f.writelines(poscar)
f.close()
