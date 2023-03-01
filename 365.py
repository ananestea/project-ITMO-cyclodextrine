from Bio import PDB
import numpy as np
import math
from math import cos, sin

import numpy as np

def magnitude(vector):  # функция длины вектора
    return math.sqrt(sum(pow(element, 2) for element in vector))

parser = PDB.PDBParser()  # для вызова методов
io = PDB.PDBIO()  # создает пдб обьект
chain_ids = []  # создает array id
chain_com = []  # array координат лиг
lig_vectors =[] # все векторы лиганды
cdx_coor =[]  # array координат cdx
lig_coor =[]
good_vectors =[] 
check_structure =[] # все векторы cdx-lig
# rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0]))  # matrix rotate moving vector onto fixed    
structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name

degree=0  # градус
theta = np.deg2rad(degree) # градус смещени
# Xrotation_matrix = np.array([[1, 0, 0], [0, cos(theta), -sin(theta)], [0, sin(theta), cos(theta)]])
# Yrotation_matrix = np.array([[cos(theta), 0, sin(theta)], [0, 1, 0], [-sin(theta), 0, cos(theta)]])
# Zrotation_matrix = np.array([[cos(theta), -sin(theta), 0], [sin(theta), cos(theta), 0], [0, 0, 1]])
rotation_matrix = np.array([[1, 0, 0], [0, cos(theta), -sin(theta)], [0, sin(theta), cos(theta)]])

for model in structure:
    for chain in model:
        for residue in chain:  # res =атомы
            chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
            if residue.get_resname() == 'DOP':  # возврат имени = dop
                for atom in residue:
                    chain_com.append(atom.coord)  # arr координат[x,y,z]

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
# new_vector=lig_vectors[0]
good_rotation = []
new_vector = np.array([0,0,0])
for degreex in range(36): # будет смещение на 3 градуса
    structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name
    theta = np.deg2rad(degreex*10) # градус смещени 
    rotation_matrix = np.array([[1, 0, 0], [0, cos(theta), -sin(theta)], [0, sin(theta), cos(theta)]])

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
                    for atom in residue:
                        atom.transform(rotation_matrix, new_vector)
                        
    for degreey in range(36): # будет смещение на 3 градуса
        theta = np.deg2rad(degreey*10) # градус смещени 
        rotation_matrix = np.array([[cos(theta), 0, sin(theta)], [0, 1, 0], [-sin(theta), 0, cos(theta)]])

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
                        for atom in residue:
                            atom.transform(rotation_matrix, new_vector)
        for degreez in range(36): # будет смещение на 3 градуса
            theta = np.deg2rad(degreez*10) # градус смещени 
            rotation_matrix = np.array([[cos(theta), -sin(theta), 0], [sin(theta), cos(theta), 0], [0, 0, 1]])

            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
                            for atom in residue:
                                atom.transform(rotation_matrix, new_vector)
                                
            for model in structure:
                for chain in model:
                    for residue in chain:  # res =атомы
                        chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
                        if residue.get_resname() == 'DOP':  # возврат имени = dop
                            for atom in residue:
                                lig_coor.append(atom.coord)

                        else:
                            for atom in residue:
                                cdx_coor.append(atom.coord)  # arr координат[x,y,z]

            for i in range(len(lig_coor)):  #  все возможные векторы между cdx и lig
                for j in range(len(cdx_coor)):
                    xvector = lig_coor[i]-cdx_coor[j]
                    if magnitude(xvector) > 2:
                        check_structure.append(1)
                    else:
                        check_structure.append(0)
                        
            lig_coor.clear()
            cdx_coor.clear()
            
            if 0 not in check_structure:
                good_rotation.append([degreex*10,degreey*10,degreez*10])
            check_structure.clear()
            
            
            # io.set_structure(structure)  # create file
            # io.save('Cyclodextrine_ligand_names' + '_' + str(degreex) + str(degreey) + str(degreez) + '.pdb')  # save file
print(len(good_rotation))
print(good_rotation)


# для вывода очень хороших векторов и очень хороших ротаций - помни про имена массивов - это ОЧЕНЬ хорошие массивы, а не просто хорошие
# if 0 not in check_structure:
#     good_vectors.append(vector)
#     good_rotation.append([degreex*10,degreey*10,degreez*10])
# check_structure.clear()