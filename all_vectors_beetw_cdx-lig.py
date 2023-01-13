from Bio import PDB
import numpy as np
import math

def magnitude(vector):  # функция длины вектора
    return math.sqrt(sum(pow(element, 2) for element in vector))

parser = PDB.PDBParser()  # для вызова методов
io = PDB.PDBIO()  # создает пдб обьект
chain_ids = []  # создает array id
chain_com = []  # array координат лиг
cdx_coor =[]  # array координат cdx
lig_coor =[] 
ligcdx_vectors =[] # все векторы cdx-lig
rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0]))  # matrix rotate moving vector onto fixed
structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name

for model in structure:
    for chain in model:
        for residue in chain:  # res =атомы
            chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
            if residue.get_resname() == 'DOP':  # возврат имени = dop
                for atom in residue:
                    chain_com.append(atom.coord)  # arr координат[x,y,z]
                    lig_coor.append(atom.coord)

            else:
                for atom in residue:
                    cdx_coor.append(atom.coord)  # arr координат[x,y,z]
print(chain_ids)
print(len(chain_com))  # 15шт
print(len(lig_coor))  # 15шт
print(len(cdx_coor))  # 86шт

for i in range(len(lig_coor)):  #  все возможные векторы между cdx и lig
        for j in range(len(cdx_coor)):
            # if i!=j:
                vector = lig_coor[i]-cdx_coor[j]
                b = float(1) * (i + 1)
                a = np.linalg.norm(vector)
                xc = -(chain_com[0][0] + (chain_com[1][0] - chain_com[0][0]) * float(b)) / a  # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
                yc = -(chain_com[0][1] + (chain_com[1][1] - chain_com[0][1]) * float(b)) / a
                zc = -(chain_com[0][2] + (chain_com[1][2] - chain_com[0][2]) * float(b)) / a
                new_vector = np.array([xc, yc, zc])
                if magnitude(new_vector) > 1.5:
                    ligcdx_vectors.append(new_vector)
print(len(ligcdx_vectors))
