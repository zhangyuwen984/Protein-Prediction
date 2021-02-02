from __future__ import print_function
import sys
import sqlite3
import string
import numpy as np
from math import sin,cos,radians,sqrt
import scipy.io as sio

def partitions(set_):
    if not set_:
        yield []
        return
    for i in range(2**len(set_)//2):
        parts = [set(), set()]
        for item in set_:
            parts[i&1].add(item)
            i >>= 1
        for b in partitions(parts[1]):
            yield [parts[0]]+b

def get_apropriate_nod_partitions(atom_numbers_by_type, nod_numbers_by_type):
    if nod_numbers_by_type.shape[0] < 13:
        nod_types = range(nod_numbers_by_type.shape[0])
        for nod_partition in partitions(nod_types):
            if len(nod_partition) <= atom_numbers_by_type.shape[0]:
                nod_numbers_by_type_for_nod_partition = np.array(list(np.array(nod_numbers_by_type[list(nod_type_group)]).sum() for nod_type_group in nod_partition)) 
                if np.all(np.sort(nod_numbers_by_type_for_nod_partition)[::-1] <= np.sort(atom_numbers_by_type)[::-1][:nod_numbers_by_type_for_nod_partition.shape[0]]):
                    yield nod_partition
    else:
        yield


def get_atom_indices_by_type(atom_type_permutation,nod_numbers_by_type,nod_partition):
    for i in np.argsort(atom_type_permutation):
        if i < len(nod_partition):
            indices = []
            for nod_type in list(nod_partition[i]):
                start_index = nod_numbers_by_type[:nod_type].sum()
                end_index = start_index + nod_numbers_by_type[nod_type]
                indices.extend(range(start_index,end_index))
            yield indices
        else:
            yield []

class Multicell:
    def __init__(self,multiplicity):
        self.multiplicity = multiplicity
        self.perm = np.random.permutation(np.arange(3))
        m = self.multiplicity - 1
        self.mults = np.dot(np.array([m&1,m&2,m&4]),np.identity(3)[self.perm]) + np.ones(3)

    def get_lattice_vectors(self,lattice_volume,cell_info):
        a,b,c,alpha,beta,gamma = cell_info
        alpha,beta,gamma = radians(alpha),radians(beta),radians(gamma)
        vnf = self.mults.reshape([3,1])*np.array([[a,0.,0.],[b*cos(gamma),b*sin(gamma),0.0],[c*cos(beta),c*(cos(alpha)-cos(beta)*cos(gamma))/(sin(gamma)),c*(sqrt(sin(beta)**2*sin(gamma)**2)-(cos(alpha)-cos(beta)*cos(gamma))**2)/(sin(gamma))]],np.float64)
        factor = ((lattice_volume)/np.dot(np.cross(vnf[0],vnf[1]),vnf[2]))**(1.0/3.0)
        return factor*vnf

    def get_coordinates(self,nod_coordinates,bonds,atom_indices_by_type,atom_2_coordinated_numbers_by_type):
        bond_center_coordinates = ((nod_coordinates[bonds[0]] + nod_coordinates[bonds[1]] + bonds[2:].T)/2)[np.random.permutation(np.arange(bonds.shape[1]))]
        shift = np.identity(3)[self.perm]
        coordinates = []
        offset_2_coordinated = 0
        for atom_type in range(len(atom_indices_by_type)):
            for i in range(self.multiplicity):
                coordinates.extend(nod_coordinates[atom_indices_by_type[atom_type]] + (i&1)*shift[0] + (i&2)*shift[1] + (i&4)*shift[2])
                a2cn = atom_2_coordinated_numbers_by_type[atom_type]
                coordinates.extend(bond_center_coordinates[offset_2_coordinated:offset_2_coordinated+a2cn] + (i&1)*shift[0] + (i&2)*shift[1] + (i&4)*shift[2])
                offset_2_coordinated += a2cn
        return np.modf(np.array(coordinates/self.mults) + np.array([3.0001,3.0001,3.0001]))[0]

def match_coordinating_numbers(atom_coordinating_numbers,nod_numbers_by_type,nod_coordinating_numbers,atom_indices_by_type,atom_2_coordinated_numbers_by_type):
    nod_coordinating_numbers_all = np.repeat(nod_coordinating_numbers,nod_numbers_by_type)
    for atom_type in range(len(atom_coordinating_numbers)):
        match_vector = np.array(list((nod_coordinating_numbers_all[atom_index] in atom_coordinating_numbers[atom_type]) for atom_index in atom_indices_by_type[atom_type]))
        if (len(atom_coordinating_numbers[atom_type]) != 0) and (not (np.all(match_vector) and ((2 in atom_coordinating_numbers[atom_type]) or (atom_2_coordinated_numbers_by_type[atom_type] == 0)))):
            return False
    return True


individual_number = int(sys.argv[1])
lattice_volume = float(sys.argv[2])
symmetry_simbol = sys.argv[3]
atom_numbers_by_type = np.array(sys.argv[4:],np.int32)
total_atom_number = atom_numbers_by_type.sum()

atom_coordinating_numbers = [[],[],[]]

contents = sio.loadmat('../../coordinatingNumbers.mat')
coordinatingNumbers = contents['coordinatingNumbers']
if coordinatingNumbers != 'NONE':
    for i in range(3):
        if i < coordinatingNumbers.shape[1]:
            atom_coordinating_numbers[i][0:0] = range(coordinatingNumbers[0][i],coordinatingNumbers[1][i])

fixed_atom_number = 0
for i in range(3):
    if (len(atom_coordinating_numbers[i]) > 0) and (min(atom_coordinating_numbers[i]) > 2):
        fixed_atom_number += atom_numbers_by_type[i]
not_fixed_atom_number = total_atom_number - fixed_atom_number

print('<CALLRESULT>')

conn = sqlite3.connect('idealnets.db')
c = conn.cursor()
structs = list(c.execute("SELECT * FROM idealnets"))# WHERE size<=?",(total_atom_number,)))
conn.close()

count = 0
while count < 1000:
    count += 1
    structure = structs[np.random.randint(len(structs))]
    size = structure[1]
    bonds_number = structure[4]
    multiplicity = fixed_atom_number//size
    reminder = fixed_atom_number%size
    if (multiplicity in [1,2,4,8]) and (reminder <= not_fixed_atom_number) and (reminder + bonds_number >= not_fixed_atom_number):
        nod_numbers_by_type = np.array(structure[2].translate({ord('['):None,ord(']'):None}).split(),np.int32)
        apropriate_nod_partitions = list(get_apropriate_nod_partitions(atom_numbers_by_type,multiplicity*nod_numbers_by_type))
        if (len(apropriate_nod_partitions) > 0)and(apropriate_nod_partitions[0]!=None):
            multicell = Multicell(multiplicity)
            lattice_vectors = multicell.get_lattice_vectors(lattice_volume,np.array(structure[6].split(),np.float64))
            nod_partition = apropriate_nod_partitions[np.random.randint(len(apropriate_nod_partitions))]
            nod_numbers_by_type_for_nod_partition = np.array(list(np.array(nod_numbers_by_type[list(nod_type_group)]).sum() for nod_type_group in nod_partition)) 
            nod_coordinates = np.array(list(line.split() for line in structure[7].translate({ord('['):u' ',ord(']'):u' '}).split("\n")),np.float64)
            bonds = np.array(list((line.split() for line in structure[8].translate({ord('['):u' ',ord(']'):u' '}).split("\n"))),np.int32).T
            nod_coordinating_numbers = np.array(structure[3].translate({ord('['):None,ord(']'):None}).split(),np.int32)
            atom_type_permutation = np.random.permutation(np.arange(atom_numbers_by_type.shape[0]))
            if np.all(multiplicity*nod_numbers_by_type_for_nod_partition <= atom_numbers_by_type[atom_type_permutation][:nod_numbers_by_type_for_nod_partition.shape[0]]):
                atom_indices_by_type = list(get_atom_indices_by_type(atom_type_permutation,nod_numbers_by_type,nod_partition))
                atom_2_coordinated_numbers_by_type = (atom_numbers_by_type - multiplicity*np.array(list(len(indices) for indices in atom_indices_by_type)))
                if match_coordinating_numbers(atom_coordinating_numbers,nod_numbers_by_type,nod_coordinating_numbers,atom_indices_by_type,atom_2_coordinated_numbers_by_type):
                    coordinates = multicell.get_coordinates(nod_coordinates,bonds,atom_indices_by_type,atom_2_coordinated_numbers_by_type)
                    if total_atom_number == len(coordinates):
                        coordinates_str = np.array_str(lattice_vectors) + u' ' + np.array_str(coordinates)
                        struct_name = structure[0][5:structure[0].find('-')]
                        print(struct_name + ' ' + coordinates_str.translate({ord('['):None,ord(']'):None,ord('\n'):None}))
                        sys.exit(0)
print(list(0 for i in range(1+3*(total_atom_number+3))))
