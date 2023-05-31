# на основе файла anather_lig.py

from Bio import PDB
import numpy as np
# from numba import jit

def magnitude(vector):  # функция длины вектора
    return np.linalg.norm(vector)

def creating_vector(vector,direction,length):
    """
    Функция, создающая векторы
    переменные: vector - вектор, direction - направление ,length - длина шага
    """

    xc = (vector[0] / direction) * length  # изменяет длину вектора, сохраняя расстояние (это вот все снизу и сверху до for model in structure)
    yc = (vector[1] / direction) * length
    zc = (vector[2] / direction) * length
    new_vector = np.array([xc, yc, zc])

    return new_vector


def id_separator(structure):
    """
    Функция, разделяющая атомы лиганда и молекулы
    переменные: structure_name - путь к файлу структуры
    """

    # parser = PDB.PDBParser()
    # structure = parser.get_structure(str(structure_name), structure_name)
    atom_count1 =[]
    atom_count2 =[]
    lig_id =[]
    rec_id =[]
    # lig_coor =[]
    # rec_coor =[]
    chain_ids =[]

    for model in structure:
        for chain in model:
            for residue in chain:  # res =атомы
                chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле 

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id() == chain_ids[1]:  # возврат имени = dop
                    for atom in residue:
                        atom_count1.append(atom.name)  # arr координат[x,y,z]
                else:
                    for atom in residue:
                        atom_count2.append(atom.name)

    if len(atom_count1) < len(atom_count2):
        lig_id = chain_ids[1]
        rec_id = chain_ids[0]
    else:
        lig_id = chain_ids[0]
        rec_id = chain_ids[1]

    # for model in structure:
    #     for chain in model:
    #         for residue in chain:
    #             if residue.get_id() == lig_id: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
    #                 for atom in residue:
    #                     lig_coor.append(atom.get_coord())
    #             else:
    #                 for atom in residue:
    #                     rec_coor.append(atom.get_coord())
    
    return lig_id, rec_id




def check_vectors(lig_vector, structure_name, steps, coeff, threshold, chain_ids):
    """
    Функция, проверяющая сообщённые вектора на вшивость
    переменные: lig_vector - вектор для проверки, structure_name - путь к файлу структуры, steps - расстояние (А)
    на которое будет смещаться структура лиганда при проверке, coeff - -1 или 1 показывает направление по которому будет
    смещаться структура лиганда, threshold - расстояние в А, которое определяет отсутсвие ковалентной связи между атомами, 
    нужно для проверки пересечения атомов возвращает массив с проверкой
    """
    parser = PDB.PDBParser()
    check_structure=[]
    coords=[]
    rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0]))
    for i in range(steps):
        structure = parser.get_structure(str(structure_name), structure_name)
        b = float(1) * (1 + i)
        a = int(coeff) * np.linalg.norm(lig_vector)
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id() == chain_ids[0]:  # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
                        for atom in residue:
                            atom.transform(rotation_matrix, creating_vector(lig_vector, a, b))

        # for model in structure:
        #     for chain in model:
        #         for residue in chain:  # res =атомы
        #             #chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле
        #             if residue.get_id() == chain_ids[0]:  # возврат имени = dop
        #                 for atom in residue:
        #                     lig_coor.append(atom.coord)
        #             else:
        #                 for atom in residue:
        #                     cdx_coor.append(atom.coord)  # arr координат[x,y,z]
        coords = coord_lig_atoms(structure, id_separator(structure)[0])

        for i in range(len(coords[0])):  # все возможные векторы между cdx и lig
            for j in range(len(coords[1])):
                xvector = coords[0][i] - coords[1][j]
                if magnitude(xvector) > threshold:
                    check_structure.append(1)
                else:
                    check_structure.append(0)

        coords[0].clear()
        coords[1].clear()

    return check_structure


# def coord_lig_atoms(structure):
#     """
#     Функция, находящая все координаты атомов для лиганда
#     переменные:
#     """

#     chain_ids, chain_com = [], []

#     for model in structure:
#         for chain in model:
#             for residue in chain:  # res =атомы
#                 chain_ids.append(residue.get_id())  # получаем все id в аrr, мы добавляем в массив данные, которые позволяют нам отличить лиганд от протеина в файле
#                 if residue.get_resname() == 'LIG':  # возврат имени = dop
#                     for atom in residue:
#                         chain_com.append(atom.coord)  # arr координат[x,y,z]

#     return chain_ids, chain_com

def coord_lig_atoms(structure, lig_id):
    """
    Функция, находящая все координаты атомов для лиганда и рецептора
    переменные: structure - структура,  lig_id - id цепочки
    """
    lig_coor =[]
    rec_coor =[]

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id() == lig_id: # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
                    for atom in residue:
                        lig_coor.append(atom.get_coord())
                else:
                    for atom in residue:
                        rec_coor.append(atom.get_coord())

    return lig_coor, rec_coor

def ligand_vectors(coeff, chain_com):
    """
    Функция, создающая все векторы для лиганда
    переменные:  coeff - -1 или 1 показывает направление по которому будет
    смещаться структура лиганда, chain_com - координаты атомов лиганда 
    """
    lig_vectors=[]
    for i in range(len(chain_com)):  # 23-33 добавляет в массив все возможные векторы
        for j in range(len(chain_com)):
            if i != j:
                vector = chain_com[i] - chain_com[j]
                b = float(1) * (1)
                a = int(coeff) * np.linalg.norm(vector)
                lig_vectors.append(creating_vector(vector, a, b))
    return lig_vectors

def creating_model(ranges, structure_name, coeff, good_vectors, chain_ids):
    """
    Функция, прорисовывающая модель
    переменные: ranges - кол-во повторов итераций, structure_name - путь к файлу структуры, coeff - -1 или 1 показывает направление по которому будет
    смещаться структура лиганд, good_vectors - хорошие векторы, прошедшие проверку,  chain_ids - id молекул
    """
    parser = PDB.PDBParser()
    rotation_matrix = PDB.rotmat(PDB.Vector([0, 0, 0]), PDB.Vector([0, 0, 0]))
    io = PDB.PDBIO()
    for i in range(ranges):
        structure = parser.get_structure(str(structure_name), structure_name)
        b = float(1) * (i + 1)
        vector = good_vectors[0]
        a = int(coeff) * np.linalg.norm(vector)

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id() == chain_ids[0]:  # если лиганд (в данном случае не протеин), то применить то, что дальше, а дальше цикл для смещения координат атомов лиганда
                        for atom in residue:
                            atom.transform(rotation_matrix, creating_vector(vector, a, b))  # right mrx, tran mx

        io.set_structure(structure)  # create file
        io.save('Cyclodextrine_ligand_names' + '_' + str(b) + '.pdb')  # save file

def make_good_vectors(vectors, structure_name, chains_ids):

    """
    Функция, находящая все хорошие векторы
    переменные: vectors - полученные векторы, structure_name - путь к файлу структуры, chain_ids - id молеку
    """

    good_vectors=[]

    for vector in vectors:
        print(vector)
        #structure = parser.get_structure(structure_name, structure_name)
        check_result_pos = check_vectors(vector, structure_name, 5, 1, 1.75, chains_ids)
        check_result_neg = check_vectors(vector, structure_name, 5, -1, 1.75, chains_ids)
        if 0 not in check_result_pos and check_result_neg:
            good_vectors.append(vector)

    return good_vectors

def main():
    parser = PDB.PDBParser()
    structure_name = "Cyclodextrine_ligand_names_good.pdb"
    structure = parser.get_structure(structure_name, structure_name)
    chains_com = coord_lig_atoms(structure, id_separator(structure)[0])
    chains_ids = id_separator(structure)
    print(len(chains_com[0]))
    print(len(chains_com[1]))
    print(chains_ids)
    vectors = ligand_vectors(1, chains_com[0])
    print(len(vectors))
    good_vectors=make_good_vectors(vectors,structure_name, chains_ids)
    creating_model(20, structure_name, 1, good_vectors, chains_ids)


if __name__ == "__main__":
    main()
