# from Bio import PDB
# import numpy as np
# import math
# from math import cos, sin

# def magnitude(vector):  # функция длины вектора
#     return math.sqrt(sum(pow(element, 2) for element in vector))

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
# # rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0]))  # matrix rotate moving vector onto fixed
# degree=0
# theta = np.deg2rad(degree) # градус смещения
# rotation_matrix =  np.array([[1, 0, 0], [0, cos(theta), -sin(theta)], [0, sin(theta), cos(theta)]])
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
# print(lig_vectors) # 210

# for x in range (len(lig_vectors)):  # 37-51  строит все возможные векторы из массива (210шт)
#     vector = lig_vectors[x]
#     # print(vector)
#     structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name
#     for i in range(5):
#         structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")
#         length = float(1) * (1+i)
#         a = np.linalg.norm(vector)
#         xc = (vector[0]/a)*length # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
#         yc = (vector[1]/a)*length
#         zc = (vector[2]/a)*length
#         new_vector = np.array([xc, yc, zc])
#         for degree in range(3):
#             structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name
#             theta = np.deg2rad(degree) # градус смещения
#             rotation_matrix =  np.array([[1, 0, 0], [0, cos(theta), -sin(theta)], [0, sin(theta), cos(theta)]])

#             for model in structure:
#                 for chain in model:
#                     for residue in chain:
#                         if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
#                             for atom in residue:
#                                 atom.transform(rotation_matrix, new_vector)

#             for model in structure:
#                 for chain in model:
#                     for residue in chain:  # res =атомы
#                         chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
#                         if residue.get_resname() == 'DOP':  # возврат имени = dop
#                             for atom in residue:
#                                 lig_coor.append(atom.coord)

#                         else:
#                             for atom in residue:
#                                 cdx_coor.append(atom.coord)  # arr координат[x,y,z]

#             for i in range(len(lig_coor)):  #  все возможные векторы между cdx и lig
#                 for j in range(len(cdx_coor)):
#                     xvector = lig_coor[i]-cdx_coor[j]
#                     if magnitude(xvector) > 1.5:
#                         check_structure.append(1)
#                     else:
#                         check_structure.append(0)
                        
#             lig_coor.clear()
#             cdx_coor.clear()
 
#     for i in range(5):
#         structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")
#         length = float(1) * (1+i)
#         a = np.linalg.norm(vector)
#         xc = (vector[0]/a)*length # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
#         yc = (vector[1]/a)*length
#         zc = (vector[2]/a)*length
#         new_vector = -np.array([xc, yc, zc])
#         for degree in range(3):
#             structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name
#             theta = np.deg2rad(degree) # градус смещения
#             rotation_matrix =  np.array([[1, 0, 0], [0, cos(theta), -sin(theta)], [0, sin(theta), cos(theta)]])

#             for model in structure:
#                 for chain in model:
#                     for residue in chain:
#                         if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
#                             for atom in residue:
#                                 atom.transform(rotation_matrix, new_vector)

#             for model in structure:
#                 for chain in model:
#                     for residue in chain:  # res =атомы
#                         chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
#                         if residue.get_resname() == 'DOP':  # возврат имени = dop
#                             for atom in residue:
#                                 lig_coor.append(atom.coord)

#                         else:
#                             for atom in residue:
#                                 cdx_coor.append(atom.coord)  # arr координат[x,y,z]

#             for i in range(len(lig_coor)):  #  все возможные векторы между cdx и lig
#                 for j in range(len(cdx_coor)):
#                     xvector = lig_coor[i]-cdx_coor[j]
#                     if magnitude(xvector) > 1.5:
#                         check_structure.append(1)
#                     else:
#                         check_structure.append(0)
                        
#             lig_coor.clear()
#             cdx_coor.clear()
# # #    print(check_structure)    
#     if 0 not in check_structure:
#         good_vectors.append(vector)
#     check_structure.clear()

# print(len(good_vectors))
# print(check_structure[1])
# print(len(lig_vectors)) # 210 шт
# print(lig_vectors) 
# print(len(good_vectors)) #1 но он соответствует lig_vetors[0]
# print(good_vectors)

# for x in range (len(good_vectors)):  # 37-51  строит все возможные векторы из массива (210шт)
#     trans_vector =good_vectors[x]
#     structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name
#     b = float(1) * (x + 1)
#     for model in structure:
#         for chain in model:
#             for residue in chain:
#                 if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
#                     for atom in residue:
#                         atom.transform(rotation_matrix, trans_vector)  # right mrx, tran mx
                            
#     io.set_structure(structure)  # create file
#     io.save('Cyclodextrine_ligand_names' + '_' + str(b) + '.pdb')  # save file
    
# for i in range(int(20)):
#     structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")
#     b = float(1) * (i + 1)
#     vector = good_vectors[0]
#     a = np.linalg.norm(vector)
#     xc = (vector[0]/a)*b # изменяет    # print(new_vector)
#     yc = (vector[1]/a)*b
#     zc = (vector[2]/a)*b
#     new_vector = np.array([xc, yc, zc])
    
#     for model in structure:
#         for chain in model:
#             for residue in chain:
#                 if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
#                     for atom in residue:
#                         atom.transform(rotation_matrix, new_vector)  # right mrx, tran mx

#     io.set_structure(structure)  # create file
#     io.save('Cyclodextrine_ligand_names' + '_' + str(b) + '.pdb')  # save file
    
# for i in range(int(20)):
#     structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")
#     b = float(1) * (i + 1)
#     vector = good_vectors[0]
#     a = np.linalg.norm(vector)
#     xc = (vector[0]/a)*b # изменяет    # print(new_vector)
#     yc = (vector[1]/a)*b
#     zc = (vector[2]/a)*b
#     new_vector = -np.array([xc, yc, zc])
    
#     for model in structure:
#         for chain in model:
#             for residue in chain:
#                 if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
#                     for atom in residue:
#                         atom.transform(rotation_matrix, new_vector)  # right mrx, tran mx

#     io.set_structure(structure)  # create file
#     io.save('Cyclodextrine_ligand_names' + '_' + str(-b) + '.pdb')  # save file
    
# print(len(good_vectors))


# гуд

# from Bio import PDB
# import numpy as np
# import math
# from math import cos, sin

# def magnitude(vector):  # функция длины вектора
#     return math.sqrt(sum(pow(element, 2) for element in vector))

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
# # rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0]))  # matrix rotate moving vector onto fixed
# degree=0
# theta = np.deg2rad(degree) # градус смещения
# rotation_matrix =  np.array([[1, 0, 0], [0, cos(theta), -sin(theta)], [0, sin(theta), cos(theta)]])
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
# print(lig_vectors) # 210

# for x in range (len(lig_vectors)):  # 37-51  строит все возможные векторы из массива (210шт)
#     vector = lig_vectors[x]
#     structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name

#     for degree in range(20):
#         structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name
#         theta = np.deg2rad(degree) # градус смещения
#         rotation_matrix = np.array([[cos(theta), -sin(theta), 0], [sin(theta), cos(theta), 0], [0, 0, 1]])
    
#         for i in range(5):
#             structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")
#             length = float(1) * (1+i)
#             a = np.linalg.norm(vector)
#             xc = (vector[0]/a)*length # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
#             yc = (vector[1]/a)*length
#             zc = (vector[2]/a)*length
#             new_vector = np.array([xc, yc, zc])
#             for model in structure:
#                 for chain in model:
#                     for residue in chain:
#                         if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
#                             for atom in residue:
#                                 atom.transform(rotation_matrix, new_vector)

#             for model in structure:
#                 for chain in model:
#                     for residue in chain:  # res =атомы
#                         chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
#                         if residue.get_resname() == 'DOP':  # возврат имени = dop
#                             for atom in residue:
#                                 lig_coor.append(atom.coord)

#                         else:
#                             for atom in residue:
#                                 cdx_coor.append(atom.coord)  # arr координат[x,y,z]

#             for i in range(len(lig_coor)):  #  все возможные векторы между cdx и lig
#                 for j in range(len(cdx_coor)):
#                     xvector = lig_coor[i]-cdx_coor[j]
#                     if magnitude(xvector) > 1.5:
#                         check_structure.append(1)
#                     else:
#                         check_structure.append(0)
                        
#             lig_coor.clear()
#             cdx_coor.clear()
    
#         # for i in range(5):
#         #     structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")
#         #     length = float(1) * (1+i)
#         #     a = np.linalg.norm(vector)
#         #     xc = (vector[0]/a)*length # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
#         #     yc = (vector[1]/a)*length
#         #     zc = (vector[2]/a)*length
#         #     new_vector = -np.array([xc, yc, zc])

#         #     for model in structure:
#         #         for chain in model:
#         #             for residue in chain:
#         #                 if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
#         #                     for atom in residue:
#         #                         atom.transform(rotation_matrix, new_vector)

#         #     for model in structure:
#         #         for chain in model:
#         #             for residue in chain:  # res =атомы
#         #                 chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
#         #                 if residue.get_resname() == 'DOP':  # возврат имени = dop
#         #                     for atom in residue:
#         #                         lig_coor.append(atom.coord)

#         #                 else:
#         #                     for atom in residue:
#         #                         cdx_coor.append(atom.coord)  # arr координат[x,y,z]

#         #     for i in range(len(lig_coor)):  #  все возможные векторы между cdx и lig
#         #         for j in range(len(cdx_coor)):
#         #             xvector = lig_coor[i]-cdx_coor[j]
#         #             if magnitude(xvector) > 1.5:
#         #                 check_structure.append(1)
#         #             else:
#         #                 check_structure.append(0)
                        
#         #     lig_coor.clear()
#         #     cdx_coor.clear()
# # #    print(check_structure)    
#     if 0 not in check_structure:
#         good_vectors.append(vector)
#     check_structure.clear()

# print(len(good_vectors))

# гуд конец

# from Bio import PDB
# import numpy as np
# import math
# from math import cos, sin

# def magnitude(vector):  # функция длины вектора
#     return math.sqrt(sum(pow(element, 2) for element in vector))

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
# # rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0]))  # matrix rotate moving vector onto fixed
# degree=0
# theta = np.deg2rad(degree) # градус смещения
# rotation_matrix =  np.array([[1, 0, 0], [0, cos(theta), -sin(theta)], [0, sin(theta), cos(theta)]])
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
# print(lig_vectors) # 210

# for x in range (len(lig_vectors)):  # 37-51  строит все возможные векторы из массива (210шт)
#     vector = lig_vectors[x]
#     structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name

#     for degree in range(20):
#         structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name
#         theta = np.deg2rad(degree) # градус смещения
#         rotation_matrix = np.array([[1, 0, 0], [0, cos(theta), -sin(theta)], [0, sin(theta), cos(theta)]])
    
#         for i in range(5):
#             structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")
#             length = float(1) * (1+i)
#             a = np.linalg.norm(vector)
#             xc = (vector[0]/a)*length # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
#             yc = (vector[1]/a)*length
#             zc = (vector[2]/a)*length
#             new_vector = np.array([xc, yc, zc])
#             for model in structure:
#                 for chain in model:
#                     for residue in chain:
#                         if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
#                             for atom in residue:
#                                 atom.transform(rotation_matrix, new_vector)

#             for model in structure:
#                 for chain in model:
#                     for residue in chain:  # res =атомы
#                         chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
#                         if residue.get_resname() == 'DOP':  # возврат имени = dop
#                             for atom in residue:
#                                 lig_coor.append(atom.coord)

#                         else:
#                             for atom in residue:
#                                 cdx_coor.append(atom.coord)  # arr координат[x,y,z]

#             for i in range(len(lig_coor)):  #  все возможные векторы между cdx и lig
#                 for j in range(len(cdx_coor)):
#                     xvector = lig_coor[i]-cdx_coor[j]
#                     if magnitude(xvector) > 1.5:
#                         check_structure.append(1)
#                     else:
#                         check_structure.append(0)
                        
#             lig_coor.clear()
#             cdx_coor.clear()
    
#         # for i in range(5):
#         #     structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")
#         #     length = float(1) * (1+i)
#         #     a = np.linalg.norm(vector)
#         #     xc = (vector[0]/a)*length # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
#         #     yc = (vector[1]/a)*length
#         #     zc = (vector[2]/a)*length
#         #     new_vector = -np.array([xc, yc, zc])

#         #     for model in structure:
#         #         for chain in model:
#         #             for residue in chain:
#         #                 if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
#         #                     for atom in residue:
#         #                         atom.transform(rotation_matrix, new_vector)

#         #     for model in structure:
#         #         for chain in model:
#         #             for residue in chain:  # res =атомы
#         #                 chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
#         #                 if residue.get_resname() == 'DOP':  # возврат имени = dop
#         #                     for atom in residue:
#         #                         lig_coor.append(atom.coord)

#         #                 else:
#         #                     for atom in residue:
#         #                         cdx_coor.append(atom.coord)  # arr координат[x,y,z]

#         #     for i in range(len(lig_coor)):  #  все возможные векторы между cdx и lig
#         #         for j in range(len(cdx_coor)):
#         #             xvector = lig_coor[i]-cdx_coor[j]
#         #             if magnitude(xvector) > 1.5:
#         #                 check_structure.append(1)
#         #             else:
#         #                 check_structure.append(0)
                        
#         #     lig_coor.clear()
#         #     cdx_coor.clear()
# # #    print(check_structure)    
#     if 0 not in check_structure:
#         good_vectors.append(vector)
#     check_structure.clear()

# print(len(good_vectors))
# print(good_vectors)

# some new

from Bio import PDB
import numpy as np
import math
from math import cos, sin

def magnitude(vector):  # функция длины вектора
    return math.sqrt(sum(pow(element, 2) for element in vector))

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
# rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0]))  # matrix rotate moving vector onto fixed
degree=0
theta = np.deg2rad(degree) # градус смещения
rotation_matrix =  np.array([[1, 0, 0], [0, cos(theta), -sin(theta)], [0, sin(theta), cos(theta)]])
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
print(lig_vectors) # 210

for x in range (len(lig_vectors)):  # 37-51  строит все возможные векторы из массива (210шт)
    vector = lig_vectors[x]
    structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name

    for degree in range(360):
        structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name
        theta = np.deg2rad(degree) # градус смещения
        rotation_matrix = np.array([[cos(theta)*cos(theta), -sin(theta)*cos(theta), sin(theta)],[sin(theta)*sin(theta)*cos(theta) + sin(theta)*cos(theta), -sin(theta)*sin(theta)*sin(theta)+cos(theta)*cos(theta), -sin(theta)*cos(theta)],[sin(theta)*sin(theta)-sin(theta)*cos(theta)*cos(theta), sin(theta)*cos(theta)+sin(theta)*sin(theta)*cos(theta), cos(theta)*cos(theta)]])

    
        for i in range(5):
            structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")
            length = float(1) * (1+i)
            a = np.linalg.norm(vector)
            xc = (vector[0]/a)*length # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
            yc = (vector[1]/a)*length
            zc = (vector[2]/a)*length
            new_vector = np.array([xc, yc, zc])
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
                    if magnitude(xvector) > 1.5:
                        check_structure.append(1)
                    else:
                        check_structure.append(0)
                        
            lig_coor.clear()
            cdx_coor.clear()
    
        for i in range(5):
            structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")
            length = float(1) * (1+i)
            a = np.linalg.norm(vector)
            xc = (vector[0]/a)*length # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
            yc = (vector[1]/a)*length
            zc = (vector[2]/a)*length
            new_vector = -np.array([xc, yc, zc])

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
                    if magnitude(xvector) > 1.5:
                        check_structure.append(1)
                    else:
                        check_structure.append(0)
                        
            lig_coor.clear()
            cdx_coor.clear()
# #    print(check_structure)    
    if 0 not in check_structure:
        good_vectors.append(vector)
    check_structure.clear()

print(len(good_vectors))
print(good_vectors)