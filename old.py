from Bio import PDB
import numpy as np


parser = PDB.PDBParser()  # для вызова методов
io = PDB.PDBIO()  # создает пдб обьект
chain_ids = []  # создает array id
chain_com = []  # array
rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0]))  # matrix rotate moving vector onto fixed
structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name
for model in structure:
    for chain in model:
        for residue in chain:  # res =атомы
            chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
            if residue.get_resname() == 'DOP':  # возврат имени = dop
                for atom in residue:
                    if atom.name == 'C6' or atom.name == 'H1': # вектор от точки C6 до точки H1
                        chain_com.append(atom.coord)  # arr [x,y,z]
print(chain_ids)
print(chain_com)

for i in range(int(20)):
    # structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")
    b = float(1) * (i + 1)
    vector = chain_com[1] - chain_com[0]
    a = np.linalg.norm(vector)
    xc = -(chain_com[0][0] + (chain_com[1][0] - chain_com[0][0]) * float(b)) / a  # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
    yc = -(chain_com[0][1] + (chain_com[1][1] - chain_com[0][1]) * float(b)) / a
    zc = -(chain_com[0][2] + (chain_com[1][2] - chain_com[0][2]) * float(b)) / a
    new_vector = np.array([xc, yc, zc])
    # print(new_vector)


    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
                    for atom in residue:
                        atom.transform(rotation_matrix, new_vector)  # right mrx, tran mx

    io.set_structure(structure)  # create file
    io.save('Cyclodextrine_ligand_names' + '_' + str(b) + '.pdb')  # save file