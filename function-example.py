# на основе файла anather_lig.py

from Bio import PDB
import numpy as np
import math

def magnitude(vector):  # функция длины вектора
    return math.sqrt(sum(pow(element, 2) for element in vector))

rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0]))  # matrix rotate moving vector onto fixed
chain_ids = []  # создает array id
chain_com = []  # array координат лиг
lig_vectors =[]  # все векторы лиганды
parser = PDB.PDBParser()

cdx_coor =[]  # array координат cdx
lig_coor =[]
good_vectors =[] 
check_structure =[] # все векторы cdx-lig

def creating_vector(vector,direction,length):

  """
  Функция, создающая векторы
  переменные: 
  """

  xc = (vector[0]/direction)*length # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
  yc = (vector[1]/direction)*length
  zc = (vector[2]/direction)*length
  new_vector = np.array([xc, yc, zc])
  return new_vector

def check_vectors(lig_vectors, structure_name, steps, coeff, name_lig, threshold):
    """
    Функция, проверяющая сообщённые вектора на вшивость
    переменные: lig_vectors - набор векторов для проверки, structure_name - путь к файлу структуры, steps - расстояние (А)
    на которое будет смещаться структура лиганда при проверке, coeff - -1 или 1 показывает направление по которому будет
    смещаться структура лиганда, name_lig - определитель(название) лиганда в файле структуры (pdb), threshold - расстояние
    в А, которое определяет отсутсвие ковалентной связи между атомами, нужно для проверки пересечения атомов
    """
    check_structure=[]
    lig_coor = []
    cdx_coor= []
    global rotation_matrix
    for x in range(len(lig_vectors)):  # 37-51  строит все возможные векторы из массива (210шт)
        vector = lig_vectors[x]
        print(vector)
        structure = parser.get_structure(str(structure_name),structure_name)  # id,file-name
        for i in range(steps):
            structure = parser.get_structure(str(structure_name),structure_name)
            b = float(1) * (1 + i)
            a = int(coeff)*np.linalg.norm(vector)
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.get_id() != chain_ids[1]:  # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
                            for atom in residue:
                                atom.transform(rotation_matrix, creating_vector((vector,a,b)))

            for model in structure:
                for chain in model:
                    for residue in chain:  # res =атомы
                        chain_ids.append(
                            residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле
                        if residue.get_resname() == str(name_lig):  # возврат имени = dop
                            for atom in residue:
                                lig_coor.append(atom.coord)
                        else:
                            for atom in residue:
                                cdx_coor.append(atom.coord)  # arr координат[x,y,z]

            for i in range(len(lig_coor)):  # все возможные векторы между cdx и lig
                for j in range(len(cdx_coor)):
                    xvector = lig_coor[i] - cdx_coor[j]
                    if magnitude(xvector) > threshold:
                        check_structure.append(1)
                    else:
                        check_structure.append(0)

            lig_coor.clear()
            cdx_coor.clear()

        return check_structure
    


def coord_lig_atoms(structure_name):
  """
  Функция, находящая все координаты атомов для лиганда
  переменные: 
  """

  structure = parser.get_structure(str(structure_name),structure_name)
  for model in structure:
      for chain in model:
          for residue in chain:  # res =атомы
              chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 
              if residue.get_resname() == 'LIG':  # возврат имени = dop
                  for atom in residue:
                      chain_com.append(atom.coord)  # arr координат[x,y,z]


def ligand_vectors(coeff):

  """
  Функция, создающая все векторы для лиганда
  переменные: 
  """
  for i in range(len(chain_com)):  # 23-33 добавляет в массив все возможные векторы
    for j in range(len(chain_com)):
      if i!=j:
        vector = chain_com[i]-chain_com[j]
        b = float(1) * (1)
        a = int(coeff)*np.linalg.norm(vector)
        lig_vectors.append(creating_vector(vector,a,b))
        # return lig_vectors


def creating_model(range, structure_name, coeff):
  """
  Функция, прорисовывающая модель
  переменные: 
  """

  io = PDB.PDBIO() 

  for i in range(int(range)):
    structure = parser.get_structure(str(structure_name),structure_name)
    b = float(1) * (i + 1)
    vector = good_vectors[0]
    a = int(coeff)*np.linalg.norm(vector)
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id() != chain_ids[1]: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
                    for atom in residue:
                        atom.transform(rotation_matrix, creating_vector(vector,a,b))  # right mrx, tran mx

    io.set_structure(structure)  # create file
    io.save('Cyclodextrine_ligand_names' + '_' + str(b) + '.pdb')  # save file 