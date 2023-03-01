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
cdx_coor = []  # array координат cdx
lig_coor = []
good_vectors = [] 
super_vectors = []
super_rotations = []
check_structure = [] # все векторы cdx-lig  
structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name

degree=0  # градус
theta = np.deg2rad(degree) # градус смещени
rotation_matrix = np.array([[1, 0, 0], [0, cos(theta), -sin(theta)], [0, sin(theta), cos(theta)]])
zero_rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0]))  # matrix rotate moving vector onto fixed

# ТО ГДЕ ИЩЕМ ХОРОШИЕ РОТАЦИИ ПО НУЛЕВОМУ ВЕКТОРУ

for model in structure:
    for chain in model:
        for residue in chain:  # res =атомы
            chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
            if residue.get_resname() == 'DOP':  # возврат имени = dop
                for atom in residue:
                    chain_com.append(atom.coord)  # arr координат[x,y,z]

good_rotations = []
zero_vector = np.array([0,0,0])

for degreex in range(36): # будет смещение на 3 градуса
    structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name
    theta = np.deg2rad(degreex*10) # градус смещени 
    rotation_matrix = np.array([[1, 0, 0], [0, cos(theta), -sin(theta)], [0, sin(theta), cos(theta)]])

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
                    for atom in residue:
                        atom.transform(rotation_matrix, zero_vector)
                        
    for degreey in range(36): # будет смещение на 3 градуса
        theta = np.deg2rad(degreey*10) # градус смещени 
        rotation_matrix = np.array([[cos(theta), 0, sin(theta)], [0, 1, 0], [-sin(theta), 0, cos(theta)]])

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
                        for atom in residue:
                            atom.transform(rotation_matrix, zero_vector)

        for degreez in range(36): # будет смещение на 3 градуса
            theta = np.deg2rad(degreez*10) # градус смещени 
            rotation_matrix = np.array([[cos(theta), -sin(theta), 0], [sin(theta), cos(theta), 0], [0, 0, 1]])

            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
                            for atom in residue:
                                atom.transform(rotation_matrix, zero_vector)
                                
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
                good_rotations.append([degreex*10,degreey*10,degreez*10])
            check_structure.clear()

print(len(good_rotations))
print(good_rotations)


# ТО ГДЕ ИЩЕМ ХОРОШИЕ ВЕКТОРЫ ПО НУЛЕВОМУ РОТАЦИОНУ

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
    # print(vector)
    structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name
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
                            atom.transform(zero_rotation_matrix, new_vector)

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
        good_vectors.append(vector)
    check_structure.clear()
print(len(good_vectors)) 

# ТО ГДЕ ПОВОРАЧИВАЕМ ХОРОШИЕ ВЕКТОРЫ ПО ХОРОШИМ РОТАЦИЯМ И ПОЛУЧАЕМ СУПЕР ХОРОШИЕ ШТУКИ

for rotation in range(len(good_rotations)):
  rotation_matrix = good_rotations[rotation],[0][1][2]

  for good_vector in range(len(good_vectors)):
      vector = good_vectors[good_vector]
      structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")  # id,file-name

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
                  if magnitude(xvector) > 2:
                      check_structure.append(1)
                  else:
                      check_structure.append(0)
                      
          lig_coor.clear()
          cdx_coor.clear()
  
    # for i in range(5):
    #     structure = parser.get_structure('Cyclodextrine_ligand_names', "Cyclodextrine_ligand_names.pdb")
    #     length = float(1) * (1+i)
    #     a = -np.linalg.norm(vector)
    #     xc = (vector[0]/a)*length # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
    #     yc = (vector[1]/a)*length
    #     zc = (vector[2]/a)*length
    #     new_vector = -np.array([xc, yc, zc])
    #     for model in structure:
    #         for chain in model:
    #             for residue in chain:
    #                 if residue.get_id() != chain_ids[0]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
    #                     for atom in residue:
    #                         atom.transform(rotation_matrix, new_vector)

    #     for model in structure:
    #         for chain in model:
    #             for residue in chain:  # res =атомы
    #                 chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
    #                 if residue.get_resname() == 'DOP':  # возврат имени = dop
    #                     for atom in residue:
    #                         lig_coor.append(atom.coord)

    #                 else:
    #                     for atom in residue:
    #                         cdx_coor.append(atom.coord)  # arr координат[x,y,z]

    #     for i in range(len(lig_coor)):  #  все возможные векторы между cdx и lig
    #         for j in range(len(cdx_coor)):
    #             xvector = lig_coor[i]-cdx_coor[j]
    #             if magnitude(xvector) > 1.5:
    #                 check_structure.append(1)
    #             else:
    #                 check_structure.append(0)
                    
    #     lig_coor.clear()
    #     cdx_coor.clear()
# #    print(check_structure)    
      if 0 not in check_structure:
          super_vectors.append(vector)
          super_rotations.append(rotation_matrix)
      check_structure.clear()

