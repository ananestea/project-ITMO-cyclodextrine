# from Bio import PDB
# import numpy as np


# parser = PDB.PDBParser()  # для вызова методов
# io = PDB.PDBIO()  # создает пдб обьект
# chain_ids = []  # создает array id
# chain_com = []  # array
# rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0]))  # matrix rotate moving vector onto fixed
# structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name
# for model in structure:
#     for chain in model:
#         for residue in chain:  # res =атомы
#             chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
#             # chain_com.append(residue.center_of_mass())
#             if residue.get_resname() == 'DOP':  # возврат имени = dop
#                 for atom in residue:
#                     if atom.name == 'C6' or atom.name == 'H1': # вектор от точки C6 до точки H1
#                         chain_com.append(atom.coord)  # arr [x,y,z]
# print(chain_ids)
# print(chain_com)
# for i in range(int(20)):
#     structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")
#     b = float(1) * (i + 1)
#     vector = chain_com[1] - chain_com[0]
#     a = np.linalg.norm(vector)
#     xc = -(chain_com[0][0] + (chain_com[1][0] - chain_com[0][0]) * float(b)) / a  # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
#     yc = -(chain_com[0][1] + (chain_com[1][1] - chain_com[0][1]) * float(b)) / a
#     zc = -(chain_com[0][2] + (chain_com[1][2] - chain_com[0][2]) * float(b)) / a
#     new_vector = np.array([xc, yc, zc])

#     for model in structure:
#         for chain in model:
#             for residue in chain:
#                 if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
#                     for atom in residue:
#                         atom.transform(rotation_matrix, new_vector)  # right mrx, tran mx

#     io.set_structure(structure)  # create file
#     io.save('Cyclodextrine_ligand_names' + '_' + str(b) + '.pdb')  # save file


# 1 version


# from Bio import PDB
# import numpy as np
# import math

# def magnitude(vector):
#     return math.sqrt(sum(pow(element, 2) for element in vector))

# parser = PDB.PDBParser()  # для вызова методов
# io = PDB.PDBIO()  # создает пдб обьект
# chain_ids = []  # создает array id
# chain_com = []  # array
# lig_coor = []  
# cdx_coor = []
# approved_vectors = []

# rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0]))  # matrix rotate moving vector onto fixed
# structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name
# for model in structure:
#     for chain in model:
#         for residue in chain:  # res =атомы
#             chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
#             if residue.get_resname() == 'DOP':  # возврат имени = dop
#                 for atom in residue:
#                     chain_com.append(atom.coord)  # arr [x,y,z]
# # print(chain_ids)
# # print(chain_com)
# for i in range(len(chain_com)):
#     print(chain_com[i])
#     for j in range(len(chain_com)):
#         print(chain_com[j])
#         if i!=j:
#             b = float(1) * (i + 1)
#             vector = chain_com[i]-chain_com[j]
#             a = np.linalg.norm(vector)
#             xc = -(chain_com[0][0] + (chain_com[1][0] - chain_com[0][0]) * float(b)) / a  # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
#             yc = -(chain_com[0][1] + (chain_com[1][1] - chain_com[0][1]) * float(b)) / a
#             zc = -(chain_com[0][2] + (chain_com[1][2] - chain_com[0][2]) * float(b)) / a
#             new_vector = np.array([xc, yc, zc])
#             # print(new_vector)

#             for model in structure:
#                 for chain in model:
#                     for residue in chain:
#                         if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
#                             for atom in residue:
#                                 atom.transform(rotation_matrix, new_vector)  # right mrx, tran mx
#                                 for model in structure:
#                                     for chain in model:
#                                         for residue in chain:  # res =атомы
#                                             chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
#                                             if residue.get_resname() == 'DOP':  # возврат имени = dop
#                                                 for atom in residue:
#                                                     lig_coor.append(atom.coord)  # arr [x,y,z]
#                                             else:
#                                                 for atom in residue:
#                                                     cdx_coor.append(atom.coord)  # arr [x,y,z]

#                                 for i in range(len(lig_coor)):
#                                     for j in range(len(lig_coor)):
#                                         # if i!=j:
#                                             b = float(1) * (i + 1)
#                                             vector = lig_coor[i]-lig_coor[j]
#                                             a = np.linalg.norm(vector)
#                                             xc = -(lig_coor[0][0] + (lig_coor[1][0] - lig_coor[0][0]) * float(b)) / a  # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
#                                             yc = -(lig_coor[0][1] + (lig_coor[1][1] - lig_coor[0][1]) * float(b)) / a
#                                             zc = -(lig_coor[0][2] + (lig_coor[1][2] - lig_coor[0][2]) * float(b)) / a
#                                             new_vector = np.array([xc, yc, zc])
#                                             if magnitude(new_vector) > 1.5:
#                                                 approved_vectors.append(new_vector)

#                                             for model in structure:
#                                                 for chain in model:
#                                                     for residue in chain:
#                                                         if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
#                                                             for atom in residue:
#                                                                 atom.transform(rotation_matrix, new_vector)  # right mrx, tran mx
#                                 for i in range(len(cdx_coor)):
#                                     for j in range(len(cdx_coor)):
#                                         # if i!=j:
#                                             b = float(1) * (i + 1)
#                                             vector = cdx_coor[i]-cdx_coor[j]
#                                             a = np.linalg.norm(vector)
#                                             xc = -(cdx_coor[0][0] + (cdx_coor[1][0] - cdx_coor[0][0]) * float(b)) / a  # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
#                                             yc = -(cdx_coor[0][1] + (cdx_coor[1][1] - cdx_coor[0][1]) * float(b)) / a
#                                             zc = -(cdx_coor[0][2] + (cdx_coor[1][2] - cdx_coor[0][2]) * float(b)) / a
#                                             new_vector = np.array([xc, yc, zc])
#                                             if magnitude(new_vector) > 1.5:
#                                                 approved_vectors.append(new_vector)

#                                             for model in structure:
#                                                 for chain in model:
#                                                     for residue in chain:
#                                                         if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
#                                                             for atom in residue:
#                                                                 atom.transform(rotation_matrix, new_vector)  # right mrx, tran mx
#                                 print(approved_vectors)

# for i in range(int(20)):
#     structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")
#     b = float(1) * (i + 1)
#     vector = approved_vectors[1] - approved_vectors[0]
#     a = np.linalg.norm(vector)
#     xc = -(approved_vectors[0][0] + (approved_vectors[1][0] - approved_vectors[0][0]) * float(b)) / a  # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
#     yc = -(approved_vectors[0][1] + (approved_vectors[1][1] - approved_vectors[0][1]) * float(b)) / a
#     zc = -(approved_vectors[0][2] + (approved_vectors[1][2] - approved_vectors[0][2]) * float(b)) / a
#     new_vector = np.array([xc, yc, zc])

#     for model in structure:
#         for chain in model:
#             for residue in chain:
#                 if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
#                     for atom in residue:
#                         atom.transform(rotation_matrix, new_vector)  # right mrx, tran mx
                        

#     io.set_structure(structure)  # create file
#     io.save('Cyclodextrine_ligand_names' + '_' + str(b) + '.pdb')  # save file






# new version


from Bio import PDB
import numpy as np
import math

def magnitude(vector):  # функция длины вектора
    return math.sqrt(sum(pow(element, 2) for element in vector))

parser = PDB.PDBParser()  # для вызова методов
io = PDB.PDBIO()  # создает пдб обьект
chain_ids = []  # создает array id
chain_com = []  # array lig
lig_coor = []  # array lig
cdx_coor = []  # array cdx
approved_vectors = [] # array 
not_approved_vectors = [] # array 

rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0]))  # matrix rotate moving vector onto fixed
structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name
for model in structure:
    for chain in model:
        for residue in chain:  # res =атомы
            chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
            if residue.get_resname() == 'DOP':  # возврат имени = dop
                for atom in residue:
                    chain_com.append(atom.coord)  # arr [x,y,z]


for model in structure:
    for chain in model:
        for residue in chain:  # res =атомы
            chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
            if residue.get_resname() == 'DOP':  # возврат имени = dop
                for atom in residue:
                    lig_coor.append(atom.coord)  # arr [x,y,z]
            else:
                for atom in residue:
                    cdx_coor.append(atom.coord)  # arr [x,y,z]
# print(len(cdx_coor))
# print(len(lig_coor))
# print(len(chain_com))

for x in range(int(1)):
    for i in range(len(chain_com)):
        for j in range(len(chain_com)):
            if i!=j:
                vector = chain_com[i]-chain_com[j]
                b = float(1) * (i + 1)
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
                

    for i in range(len(cdx_coor)):
        for j in range(len(lig_coor)):
                b = float(1) * (i + 1)
                vector = cdx_coor[i]-lig_coor[j]
                a = np.linalg.norm(vector)
                xc = -(chain_com[0][0] + (chain_com[1][0] - chain_com[0][0]) * float(b)) / a  # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
                yc = -(chain_com[0][1] + (chain_com[1][1] - chain_com[0][1]) * float(b)) / a
                zc = -(chain_com[0][2] + (chain_com[1][2] - chain_com[0][2]) * float(b)) / a
                new_vector = np.array([xc, yc, zc])

                for model in structure:
                    for chain in model:
                        for residue in chain:
                            if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
                                for atom in residue:
                                    atom.transform(rotation_matrix, new_vector)  # right mrx, tran mx
                                    
                if magnitude(new_vector) > 1.5:
                    approved_vectors.append(new_vector)
                else:
                    not_approved_vectors.append(new_vector)
print(len(approved_vectors))
print(len(not_approved_vectors))

# for x in range(int(20)):
for i in range(len(approved_vectors)):
    for j in range(len(approved_vectors)):
        if i!=j:
            vector = approved_vectors[i]-approved_vectors[j]
            # structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")
            b = float(1) * (i + 1)
            a = np.linalg.norm(vector)
            xc = -(chain_com[0][0] + (chain_com[1][0] - chain_com[0][0]) * float(b)) / a  # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
            yc = -(chain_com[0][1] + (chain_com[1][1] - chain_com[0][1]) * float(b)) / a
            zc = -(chain_com[0][2] + (chain_com[1][2] - chain_com[0][2]) * float(b)) / a
            new_vector = np.array([xc, yc, zc])

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id() !=chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
                        for atom in residue:
                            atom.transform(rotation_matrix, new_vector)  # right mrx, tran mx

        io.set_structure(structure)  # create file
        io.save('Cyclodextrine_ligand_names' + '_' + str(b) + '.pdb')  # save file


                                # for i in range(len(cdx_coor)):
                                #     for j in range(len(cdx_coor)):
                                #         # if i!=j:
                                #             b = float(1) * (i + 1)
                                #             vector = cdx_coor[i]-cdx_coor[j]
                                #             a = np.linalg.norm(vector)
                                #             xc = -(cdx_coor[0][0] + (cdx_coor[1][0] - cdx_coor[0][0]) * float(b)) / a  # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
                                #             yc = -(cdx_coor[0][1] + (cdx_coor[1][1] - cdx_coor[0][1]) * float(b)) / a
                                #             zc = -(cdx_coor[0][2] + (cdx_coor[1][2] - cdx_coor[0][2]) * float(b)) / a
                                #             new_vector = np.array([xc, yc, zc])
                                #             if magnitude(new_vector) > 1.5:
                                #                 approved_vectors.append(new_vector)

                                #             for model in structure:
                                #                 for chain in model:
                                #                     for residue in chain:
                                #                         if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
                                #                             for atom in residue:
                                #                                 atom.transform(rotation_matrix, new_vector)  # right mrx, tran mx
                                # print(approved_vectors)

# for i in range(int(20)):
#     structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")
#     b = float(1) * (i + 1)
#     vector = approved_vectors[1] - approved_vectors[0]
#     a = np.linalg.norm(vector)
#     xc = -(chain_com[0][0] + (chain_com[1][0] - chain_com[0][0]) * float(b)) / a  # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
#     yc = -(chain_com[0][1] + (chain_com[1][1] - chain_com[0][1]) * float(b)) / a
#     zc = -(chain_com[0][2] + (chain_com[1][2] - chain_com[0][2]) * float(b)) / a
#     new_vector = np.array([xc, yc, zc])

#     for model in structure:
#         for chain in model:
#             for residue in chain:
#                 if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
#                     for atom in residue:
#                         atom.transform(rotation_matrix, new_vector)  # right mrx, tran mx
                        

    # io.set_structure(structure)  # create file
    # io.save('Cyclodextrine_ligand_names' + '_' + str(b) + '.pdb')  # save file




                         # here we go again

