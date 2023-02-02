# from Bio import PDB
# import numpy as np
# import math

# parser = PDB.PDBParser()  # для вызова методов
# io = PDB.PDBIO()  # создает пдб обьект
# chain_ids = []  # создает array id
# chain_com = []  # array координат лиг
# lig_vectors =[] # все векторы лиганды
# # all_lig_vectors =[] # все векторы лиганды

# cdx_coor =[]  # array координат cdx
# lig_coor =[]
# good_vectors =[] 
# check_structure =[] # все векторы cdx-lig
# rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0]))  # matrix rotate moving vector onto fixed
  
# structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name

# for model in structure:
#     for chain in model:
#         for residue in chain:  # res =атомы
#             chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
#             if residue.get_resname() == 'DOP':  # возврат имени = dop
#                 for atom in residue:
#                     chain_com.append(atom.coord)  # arr координат[x,y,z]
# # print(chain_ids)
# # print(len(chain_com))  # 15шт

# for i in range(len(chain_com)):  # 23-33 добавляет в массив все возможные векторы
#         for j in range(len(chain_com)):
#             if i!=j:
#                 vector = chain_com[i]-chain_com[j]
#                 b = float(1) * (1)
#                 a = np.linalg.norm(vector)
#                 xc = (vector[0]/a)*b # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
#                 yc = (vector[1]/a)*b
#                 zc = (vector[2]/a)*b
#                 new_vector = np.array([xc, yc, zc])
                
#                 lig_vectors.append(new_vector)
# # print(lig_vectors) # 210
# print(lig_vectors[0])
# new_vector=lig_vectors[0]

# structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name
# for model in structure:
#     for chain in model:
#         for residue in chain:
#             if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
#                 for atom in residue:
#                     atom.transform(rotation_matrix, new_vector)
# io.set_structure(structure)  # create file
# io.save('Cyclodextrine_ligand_names' + '_' + str(b) + '.pdb')  # save file



from Bio import PDB
import numpy as np

from math import cos, sin

import numpy as np

theta = - np.deg2rad(30)

rotation_matrix = np.array([[1, 0, 0], [0, cos(theta), -sin(theta)], [0, sin(theta), cos(theta)]])


# Xrotation_matrix = np.array([[1, 0, 0], [0, cos(theta), -sin(theta)], [0, sin(theta), cos(theta)]])
# Yrotation_matrix = np.array([[cos(theta), 0, sin(theta)], [0, 1, 0], [-sin(theta), 0, cos(theta)]])
# Zrotation_matrix = np.array([[cos(theta), -sin(theta), 0], [sin(theta), cos(theta), 0], [0, 0, 1]])

parser = PDB.PDBParser()  # для вызова методов
io = PDB.PDBIO()  # создает пдб обьект
chain_ids = []  # создает array id
chain_com = []  # array координат лиг
lig_vectors =[] # все векторы лиганды
# all_lig_vectors =[] # все векторы лиганды

cdx_coor =[]  # array координат cdx
lig_coor =[]
good_vectors =[] 
check_structure =[] # все векторы cdx-lig
# rotation_matrix = PDB.rotmat(PDB.Vector([0, 1, 0]), PDB.Vector([0, 0, 1]))  # matrix rotate moving vector onto fixed
# ab,bc,cd = 0, 0, 0 
# ba,cb,dc = 0,0,0
# rotation_matrix = PDB.rotmat(PDB.Vector([ab, bc, cd]), PDB.Vector([ba, cb , dc]))
# rotation_matrix = PDB.rotaxis2m(360, PDB.Vector([1, 1, 0]))  # matrix rotate moving vector onto fixed
# p=PDB.Vector([0, 0, 0])
# q=PDB.Vector([0, 0, 0])
# rotation_matrix = PDB.rotmat(p, q)  # matrix rotate moving vector onto fixed


# rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0]))  # matrix rotate moving vector onto fixed    
structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name

for model in structure:
    for chain in model:
        for residue in chain:  # res =атомы
            chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
            if residue.get_resname() == 'DOP':  # возврат имени = dop
                for atom in residue:
                    chain_com.append(atom.coord)  # arr координат[x,y,z]
# print(chain_ids)
# print(len(chain_com))  # 15шт
# for ab in range(int(3)):
#     for bc in range(int(3)):
#         for cd in range(int(3)):
#             for ba in range(int(3)):
#                 for cb in range(int(3)):
#                     for dc in range(int(3)):

for i in range(len(chain_com)):  # 23-33 добавляет в массив все возможные векторы
        for j in range(len(chain_com)):
            if i!=j:
                vector = chain_com[i]-chain_com[j]
                b = float(1) * (1)
                a = np.linalg.norm(vector)
                xc = (vector[0]/a)*b # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
                yc = (vector[1]/a)*b
                zc = (vector[2]/a)*b
                new_vector = np.array([xc, yc, zc])
                
                lig_vectors.append(new_vector)
# print(lig_vectors) # 210
print(lig_vectors[0])
new_vector=lig_vectors[5]

# ab,bc,cd = 0, 0, -1 
# ba,cb,dc = 0,0,0
# rotation_matrix = PDB.rotmat(PDB.Vector([ab, bc, cd]), PDB.Vector([ba, cb , dc]))
# print(rotation_matrix)

# for ab in range(3):
#     w=ab
#     for bc in range(3):
#         e=bc
#         for cd in range(3):
#             r=cd
#             for ba in range(3):
#                 t=ba
#                 for cb in range(3):
#                     y=cb
#                     for dc in range(3):
#                         u=dc

structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
                for atom in residue:
                    atom.transform(rotation_matrix, new_vector)
io.set_structure(structure)  # create file
io.save('Cyclodextrine_ligand_names' + '_' + str(80) + '.pdb')  # save file