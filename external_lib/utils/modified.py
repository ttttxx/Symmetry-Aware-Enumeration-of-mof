import numpy as np
import os
import re
import math
from pymatgen.core.lattice import Lattice

#缩放因子阈值
scale_threshold = 3.0

vertices = ('H' , 'He', 'Li', 'Be', 'B' , 'C' , 'N' , 'O' , 'F' , 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar',
	  'K' , 'Ca', 'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
	  'Rb', 'Sr', 'Y' , 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I' , 'Xe',
	  'Cs', 'Ba', 'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr',
	  'Ra', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th',
	  'Pa', 'U' , 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'FG', 'X' )

def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

def nn(string):
	#匹配单个非字母，并将单个非字母替换
	return re.sub('[^a-zA-Z]','', string)

def nl(string):
	#匹配单个非数字，并将单个非数字替换
	return re.sub('[^0-9]','', string)

def isvert(line):

    if len(line) >= 5:
        if nn(line[0]) in vertices and line[1] in vertices and False not in map(isfloat, line[2:5]):
            # 此句if语句判断：第一个和第二个是否为集合中的顶点以及三、四、五是否全为浮点数（坐标）
            # 以上所有条件构成了从cif文件中判断是否为顶点坐标的条件
            return True
        else:
            return False

def isedge(line):

    if len(line) >= 5:
        if nn(line[0]) in vertices and nn(line[1]) in vertices and isfloat(line[2]) and line[-1] in ('S', 'D', 'T', 'A'):
            # 此句if语句判断：第一个和第二个是否为集合中的顶点，第三个是否为浮点数（边距），以及最后一个是否为化学键的类型
            # 以上所有条件构成了从cif文件中判断是否为边的信息的条件
            return True
        else:
            return False

#计算改性基团的原子在结果材料中的起始索引
def calculate_index(cif_path,elements_set):

    elements_dict = {}  # 创建一个空字典
    for value in elements_set:
        elements_dict[value] = 1

    with open(cif_path, 'r') as f:
        lines = f.readlines()

    for line in lines:
        tokens = line.split()
        if isvert(tokens) and tokens[1] in elements_dict:#如果是顶点:
            element = re.sub('[^a-zA-Z]', '',tokens[0])
            index = int(''.join(filter(str.isdigit,tokens[0]))) # 提取字符串中的数字部分)
            if index > elements_dict[tokens[1]]:  # 使用正则表达式查找字符串中的所有数字部分
                elements_dict[tokens[1]] = index

    return elements_dict

#统计改性基团中的元素类型
def get_elements_type(cif_path):

    elements_set=set()
    with open(cif_path, 'r') as f:
        lines = f.readlines()
    for line in lines:
        tokens = line.split()
        if isvert(tokens):#如果是顶点:
            # print(tokens)
            if tokens[1].isalpha():
                elements_set.add(tokens[1])

    return elements_set

#读取改性基团的cif文件,分别建立定位点和其余原子的的坐标字典
def read_functional_group(functional_group_file):
    with open(functional_group_file, 'r') as f:
        lines = f.readlines()

    ref_node_dict = {}
    other_node_dict = {}

    flag =False
    for line in lines:
        s = line.split()
        if "_modified_remove_atom_site_label" in line:
            flag = True
            continue
        if flag:
            ref_node_dict[s[0]]=[]

    for line in lines:
        s = line.split()
        if isvert(s) and s[0] not in ref_node_dict:  # 如果是顶点且不是定位点
            other_node_dict[s[0]] = [float(s[2]),float(s[3]),float(s[4])]
        if isvert(s) and s[0] in ref_node_dict:
            ref_node_dict[s[0]] = float(s[2]),float(s[3]),float(s[4])

    return ref_node_dict,other_node_dict

#读取母体材料的cif文件,分别建立定位点和删除原子的坐标字典
def read_mofcof(real_path):
    with open(real_path, 'r') as f:
        lines = f.readlines()

    ref_node_dict = {}
    remove_node_dict = {}
    remove_node_list = []

    flag = False
    for line in lines:
        s = line.split()
        if "_modified_remove_atom_site_label" in line:
            flag = True
            continue
        if flag:
            ref_node_dict[s[0]] = []
            if len(s) > 1:
                remove_node_dict[s[1]] = []
                remove_node_list.append(s[1])

    for line in lines:
        s = line.split()
        if isvert(s) and s[0] in ref_node_dict:
            ref_node_dict[s[0]] = float(s[2]),float(s[3]),float(s[4])
        if isvert(s) and s[0] in remove_node_dict:
            remove_node_dict[s[0]] = float(s[2]),float(s[3]),float(s[4])

    return ref_node_dict,remove_node_dict,remove_node_list

def get_unit_cell(cif_file):

    with open(cif_file, 'r') as cif:
        cif = cif.read()
        cif = filter(None, cif.split('\n'))

    for line in cif:
        l = line.split()
        if '_cell_length_a' in line:
            a = float(l[1])
        if '_cell_length_b' in line:
            b = float(l[1])
        if '_cell_length_c' in line:
            c = float(l[1])
        if '_cell_angle_alpha' in line:
            alpha = float(l[1])
        if '_cell_angle_beta' in line:
            beta = float(l[1])
        if '_cell_angle_gamma' in line:
            gamma = float(l[1])
            break

    unit_cell = [a, b, c, alpha, beta, gamma]
    return unit_cell

def frac_to_cart_ori(unit_cell, frac_coords):
    # Convert degrees to radians
    a, b, c, alpha, beta, gamma = unit_cell
    alpha = np.radians(alpha)
    beta = np.radians(beta)
    gamma = np.radians(gamma)

    volume = np.sqrt(1 - np.cos(alpha) ** 2 - np.cos(beta) ** 2 - np.cos(gamma) ** 2 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma))

    crystal_vectors = np.array([[a, b * np.cos(gamma), c * np.cos(beta)],
                                [0, b * np.sin(gamma), c * (np.cos(alpha) -np.cos(beta) * np.cos(gamma)) / np.sin(gamma)],
                                [0, 0, c * volume / np.sin(gamma)]])

    cart_coords = np.dot(frac_coords, crystal_vectors)

    return cart_coords

def frac_to_cart(unit_cell, frac_coords):

    # 定义晶胞参数
    a, b, c, alpha, beta, gamma = unit_cell

    # 使用晶胞参数创建晶格对象
    lattice = Lattice.from_parameters(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma)
    # 转换为笛卡尔坐标
    cartesian_coords = lattice.get_cartesian_coords(frac_coords)
    # print(cartesian_coords)
    return cartesian_coords

def cart_to_frac_ori(unit_cell, cart_coords):
    # Convert degrees to radians
    a, b, c, alpha, beta, gamma = unit_cell
    alpha = np.radians(alpha)
    beta = np.radians(beta)
    gamma = np.radians(gamma)

    volume = np.sqrt(1 - np.cos(alpha) ** 2 - np.cos(beta) ** 2 -
                     np.cos(gamma) ** 2 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma))

    # Matrix of crystal vectors
    crystal_vectors = np.array([[a, b * np.cos(gamma), c * np.cos(beta)],
                                [0, b * np.sin(gamma), c * (np.cos(alpha) -
                                                            np.cos(beta) * np.cos(gamma)) / np.sin(gamma)],
                                [0, 0, c * volume / np.sin(gamma)]])

    # Invert the crystal vectors matrix
    inv_crystal_vectors = np.linalg.pinv(crystal_vectors)

    frac_coords = np.dot(cart_coords, inv_crystal_vectors)

    return frac_coords

def cart_to_frac(unit_cell, cart_coords):
    # 定义晶胞参数
    a, b, c, alpha, beta, gamma = unit_cell

    # 使用晶胞参数创建晶格对象
    lattice = Lattice.from_parameters(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma)

    fractional_coords2 = lattice.get_fractional_coords(cart_coords)

    return fractional_coords2

def kabsch(P, Q):
    """
    Kabsch算法，寻找P到Q的最优旋转矩阵
    P和Q是NxM的矩阵，表示N个M维的点
    返回平移向量、缩放系数和旋转矩阵
    """
    # 计算质心
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)

    # 中心化
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q

    # 计算协方差矩阵
    C = np.dot(np.transpose(P_centered), Q_centered)
    # print("协方差矩阵如下：",C)

    # 使用SVD进行分解
    V, S, Wt = np.linalg.svd(C)

    # 确保它是右手坐标系
    if (np.linalg.det(V) * np.linalg.det(Wt)) < 0.0:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # 创建旋转矩阵U
    U = np.dot(V, Wt)
    # print("旋转矩阵如下：",U)

    # 计算缩放因子
    # scale = np.sum(S) / np.sum(np.square(P_centered))
    scale = np.sum(np.linalg.norm(Q_centered, axis=1)) / np.sum(np.linalg.norm(P_centered, axis=1))
    # print("缩放因子如下：",scale)
    # 计算平移向量
    translation = centroid_Q - scale * np.dot(centroid_P, U)
    # print("平移向量如下：",translation)

    return translation, scale, U

def find_new_bonds(cif_path,node,func_ref_node):
    with open(cif_path, 'r') as f:
        lines = f.readlines()
    new_bonds=[]
    for line in lines:
        s = line.split()
        if isedge(s):
            if s[0]!=node and s[1]!=node and s[0] not in func_ref_node and s[1] not in func_ref_node:
                new_bonds.append(line)
            if s[0]==node and s[1] not in func_ref_node:
                connect_node=s[1]
            if s[1]==node and s[0] not in func_ref_node:
                connect_node=s[0]
    return new_bonds,connect_node

def modify_cif_file(cif_path,file,remove_node,new_node,new_bonds,output_folder):
    # 创建改性输出文件夹
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    with open(cif_path, 'r') as f:
        lines = f.readlines()

    new_lines = []  # A new list to store lines you want to keep
    flag = False
    flag_ = False

    for line in lines:
        s = line.split()
        if '_atom_site_charge' in line:
            flag = True
            new_lines.append(line)
            # 将列表中的元素转换为字符串类型
            result = ' '.join(new_node)
            new_lines.append(result)
            continue
        if flag and 'loop_' not in line:
            if s[0] in remove_node:
                continue  # Skip the current line
        if 'loop_' in line:
            flag = False
            flag_ = False
        if '_ccdc_geom_bond_type' in line:
            flag_ = True
            new_lines.append(line)
            # 将列表中的元素转换为字符串类型
            result = ' '.join(new_bonds)
            new_lines.append(result)
            continue
        if flag_ and 'loop_' not in line:
            if s[0] in remove_node or s[1] in remove_node:
                continue
        new_lines.append(line)  # Add line to new_lines if it should be kept

    cif_name = os.path.basename(cif_path)  # 获取文件名（带扩展名）
    cif_name_without_ext = os.path.splitext(cif_name)[0]  # 去除扩展名

    new_file_path = os.path.join(output_folder,cif_name_without_ext+"-"+file)
    # print(new_lines)
    with open(new_file_path, 'w') as f:
        f.writelines(new_lines)

    return new_file_path

def get_key_index(dictionary, target_key):
    keys = list(dictionary.keys())
    for i, key in enumerate(keys):
        if key == target_key:
            return i
    return -1  # 如果未找到目标键，返回 -1 或其他适当的值

def distance(point1, point2):
    # 计算两个点之间的欧氏距离
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(point1, point2)))

def find_point_with_max_distance(points):
    # 遍历每个点，计算其与其他两个点的距离，并找到最大距离对应的点
    max_distance = 0
    max_distance_point = None

    for i, point in enumerate(points):
        # 计算与其他两个点的距离总和
        total_distance = 0
        for j, other_point in enumerate(points):
            if j != i:
                total_distance += distance(point, other_point)

        # 更新最大距离和对应的点
        if total_distance > max_distance:
            max_distance = total_distance
            max_distance_point = point

    return max_distance_point

def find_dimension_with_max_difference(outlier_point, point1, point2):
    dimensions = len(outlier_point)

    max_difference = 0
    max_difference_dimension = None

    for dim in range(dimensions):
        difference = abs(outlier_point[dim] - (point1[dim] + point2[dim]) / 2)
        if difference > max_difference:
            max_difference = difference
            max_difference_dimension = dim

    return max_difference_dimension
    
def renew_mofcoffile(real_path):
    
    with open(real_path, 'r') as f:
        lines = f.readlines()
        
    # print(type(lines))
    
    modified_lines = lines[:-3]
    
    f.close()
    
    with open(real_path, 'w') as f:
        f.writelines(modified_lines)

    f.close()

# 函数参数分别为：存放改性后的文件的文件夹、改性基团文件夹、材料结构文件夹
def modify_MOF_COF(result_path, functions_path, mofcof_path):

    functions_file = os.listdir(functions_path)
    functions_file.sort()
    print('Start Modify Function')
    print('funcions_file: ', functions_file)
    print()

    # mofcof_file = os.listdir(mofcof_path)
    mofcof_file = mofcof_path
    # mofcof_file.sort()
    print('mofcof_file: ', mofcof_file)
    print()

    # 遍历每一个mof/cof结构
    # for mcfile in mofcof_file:

    # print('MOF/COF: ', mcfile)
    # 材料cif文件的真正路径
    # real_path = os.path.join(mofcof_path, mcfile)
    real_path = mofcof_file

    # if os.path.isfile(real_path) == False or mcfile == 'debug.log':
    #     continue

    # 遍历每一个改性基团文件
    for file in functions_file:

        new_node = []
        new_bonds = []

        if file == '.DS_Store':
            continue
        print('Function Group: ', file)

        # 获取改性基团中所有元素种类
        elements_set = get_elements_type(os.path.join(functions_path, file))
        # 获取将改性基团加入母体材料中时的起始下标
        elements_index_dict = calculate_index(real_path, elements_set)
        # print(elements_index_dict)

        # 读取改性基团的定位原子，以及剩余原子
        func_ref_node, func_other_node = read_functional_group(os.path.join(functions_path, file))
        func_ref_node_list = [node for node in func_ref_node]
        # print("func_ref_node",func_ref_node)
        # print("func_other_node",func_other_node)
        corresponding_atoms = {key: "" for key in func_other_node}
        # print("corresponding_atoms",corresponding_atoms)

        mof_ref_node, remove_node, remove_node_list= read_mofcof(real_path)
        mof_ref_node_list = [node for node in mof_ref_node]
        # print("mof_ref_node",mof_ref_node)
        # print("remove_node",remove_node)
        # print("remove_node",len(remove_node))
        # print("remove_node_list",remove_node_list)

        mof_index_nodes_num = len(mof_ref_node)
        func_index_nodes_num = len(func_ref_node)
        if mof_index_nodes_num == 0 :
            renew_mofcoffile(real_path)
            print("母体材料中不存在改性位点，改性失败！")
            break

        if mof_index_nodes_num % func_index_nodes_num != 0:
            print("改性基团和母体材料结构不匹配，改性失败")
            continue

        # 计算得到改性基团定位点原子的笛卡尔坐标
        func_unitcell = get_unit_cell(functions_path + '/' + file)
        # print("func_unit_cell",func_unitcell)
        ref_fcoords = [value for value in func_ref_node.values()]
        # print("ref_fcoords",ref_fcoords)
        ref_fcoords_np = np.array(ref_fcoords)
        ref_ccoords_np = frac_to_cart(func_unitcell, ref_fcoords_np)
        # 计算得到改性基团中剩余原子的笛卡尔坐标
        other_fcoords = [value for value in func_other_node.values()]
        other_fcoords_np = np.array(other_fcoords)
        other_ccoords_np = frac_to_cart(func_unitcell, other_fcoords_np)

        grouped_list = [mof_ref_node_list[i:i + func_index_nodes_num] for i in
                        range(0, mof_index_nodes_num, func_index_nodes_num)]
        # print('grouped_list: ', grouped_list)

        for grouped in grouped_list:
            # print("grouped:",grouped)
            # 计算得到母体材料中每一组定位原子的笛卡尔坐标
            group_unitcell = get_unit_cell(real_path)
            group_fcoord = [mof_ref_node[node] for node in grouped]
            group_np_fcoord = np.array(group_fcoord, dtype=float)
            group_np_ccoord = frac_to_cart(group_unitcell, group_np_fcoord)

            translation, scale, rotation = kabsch(ref_ccoords_np, group_np_ccoord)
            other_result_ccoords = scale * np.dot(other_ccoords_np, rotation) + translation
            other_result_ffoords = cart_to_frac(group_unitcell, other_result_ccoords)

            # print(scale)
            if scale < scale_threshold:
                # 新原子
                i = 0
                for node in func_other_node.keys():
                    element = re.sub('[^a-zA-Z]', '', node)
                    new_index = elements_index_dict[element]
                    new_index += 1
                    elements_index_dict[element] = new_index
                    element_with_index = element + str(new_index)
                    corresponding_atoms[node] = element_with_index
                    result = np.array2string(other_result_ffoords[i], separator=', ')[1:-1]
                    result = result.replace(',', '  ')
                    line_in_cif = element_with_index + "\t" + element + "  " + result + '\t0\n'
                    new_node.append(line_in_cif)
                    i += 1
                # print(corresponding_atoms)
                # 新化学键
                func_new_bonds, connect_node = find_new_bonds(functions_path + '/' + file, func_ref_node_list[1],func_ref_node)
                # 替换列表中存在于字典键中的部分为对应的值
                for bond in func_new_bonds:
                    words = bond.split()
                    new_words = []
                    for word in words:
                        if word in corresponding_atoms:
                            new_words.append(corresponding_atoms[word])
                        else:
                            new_words.append(word)
                    new_bond = ' '.join(new_words)
                    new_bonds.append(new_bond + '\n')
                # C-C键
                mof_c_coord = group_np_ccoord[1]  # 母体材料中定位点原子（连接点）的笛卡尔坐标
                fun_c_index = get_key_index(func_other_node, connect_node)
                fun_c_coord = other_result_ccoords[fun_c_index]  # 改性基团中连接点的笛卡尔坐标
                # print("分数坐标如下：")
                # print(grouped[1], "mof_c_coord", group_np_fcoord[1])
                # print(corresponding_atoms[connect_node], "fun_c_coord", other_result_ffoords[fun_c_index])
                # print("笛卡尔坐标如下：")
                # print(grouped[1], "mof_c_coord", mof_c_coord)
                # print(corresponding_atoms[connect_node], "fun_c_coord", fun_c_coord)
                dist = np.linalg.norm(mof_c_coord - fun_c_coord)   #计算化学键长
                line_in_cif = grouped[1] + "  " + corresponding_atoms[connect_node] + "  "+f"{dist:.3f}"+" " + ' . ' + ' S\n'
                new_bonds.append(line_in_cif)
                # print(line_in_cif)
                # print("-------------------------------------------------")
            else:
                # 三个点的三维坐标
                group_fcoord = [mof_ref_node[node] for node in grouped]
                #取代位点
                substitutional_point = list(group_fcoord[1])
                #定位点
                positioning_point1 = list(group_fcoord[0])
                positioning_point2 = list(group_fcoord[2])

                # 遍历每一个维度  
                for dim in range(3):  
                    # 计算差值  
                    diff = abs(positioning_point1[dim] - substitutional_point[dim])  
                    
                    # 如果差值大于0.7  
                    if diff > 0.7:  
                        pos = substitutional_point[dim]  
                        # 根据条件判断增加还是减少 
                        if abs(1 - pos) < abs(0 - pos):
                            positioning_point1[dim] += 1
                        else:
                            positioning_point1[dim] -= 1
                
                # 遍历每一个维度  
                for dim in range(3):  
                    # 计算差值  
                    diff = abs(positioning_point2[dim] - substitutional_point[dim])  
                    
                    # 如果差值大于0.7  
                    if diff > 0.7:  
                        pos = substitutional_point[dim]  
                        # 根据条件判断增加还是减少，并设置对应的标志位  
                        if abs(1 - pos) < abs(0 - pos):
                            positioning_point2[dim] += 1
                        else:
                            positioning_point2[dim] -= 1
                
                new_group_fcoord = [tuple(positioning_point1),tuple(substitutional_point),tuple(positioning_point2)]

                # 计算重心
                center_of_mass = np.mean(new_group_fcoord, axis=0)
                # print(center_of_mass)

                cell_center = [0.5, 0.5, 0.5]  # 晶胞中心坐标

                # 计算平移向量
                translation_vector = cell_center - center_of_mass

                # 将边缘原子团应用平移向量进行平移
                translated_atoms = new_group_fcoord + translation_vector
                # print("平移后原子的坐标:", translated_atoms)

                # 利用平移后的坐标进行点云配对
                group_np_fcoord = np.array(translated_atoms, dtype=float)
                group_np_ccoord = frac_to_cart(group_unitcell, group_np_fcoord)

                translation, scale, rotation = kabsch(ref_ccoords_np, group_np_ccoord)
                other_result_ccoords = scale * np.dot(other_ccoords_np, rotation) + translation
                other_result_ffoords = cart_to_frac(group_unitcell, other_result_ccoords)
            
                i = 0
                for node in func_other_node.keys():
                    element = re.sub('[^a-zA-Z]', '', node)
                    new_index = elements_index_dict[element]
                    new_index += 1
                    elements_index_dict[element] = new_index
                    element_with_index = element + str(new_index)
                    corresponding_atoms[node] = element_with_index
                    other_result_ffoords[i] -= translation_vector
                    result = np.array2string(other_result_ffoords[i], separator=', ')[1:-1]
                    result = result.replace(',', '  ')
                    line_in_cif = element_with_index + "\t" + element + "  " + result + '\t0\n'
                    new_node.append(line_in_cif)
                    i += 1
                # print(corresponding_atoms)
                # 新化学键
                func_new_bonds, connect_node = find_new_bonds(functions_path + '/' + file, func_ref_node_list[1],
                                                            func_ref_node)
                # 替换列表中存在于字典键中的部分为对应的值
                for bond in func_new_bonds:
                    words = bond.split()
                    new_words = []
                    for word in words:
                        if word in corresponding_atoms:
                            new_words.append(corresponding_atoms[word])
                        else:
                            new_words.append(word)
                    new_bond = ' '.join(new_words)
                    new_bonds.append(new_bond + '\n')
                # C-C键
                mof_c_coord = group_np_ccoord[1]  # 母体材料中定位点原子（连接点）的笛卡尔坐标
                fun_c_index = get_key_index(func_other_node, connect_node)
                fun_c_coord = other_result_ccoords[fun_c_index]  # 改性基团中连接点的笛卡尔坐标
                dist = np.linalg.norm(mof_c_coord - fun_c_coord)  # 计算化学键长
                line_in_cif = grouped[1] + "  " + corresponding_atoms[
                    connect_node] + "  " + f"{dist:.3f}" + " " + ' . ' + ' S\n'
                new_bonds.append(line_in_cif)

        # print("new_node",new_node)
        # print("new_bonds",new_bonds)
        modify_cif_file(real_path, file, remove_node, new_node, new_bonds, os.path.join(result_path, file))
        print(real_path + "利用" + file + "改性完成")

# if __name__ == '__main__':
    
#     # 改性基团文件夹
#     functions_path = 'functional_groups'
#     # 材料结构文件夹
#     mofcof_path = '/mnt/data/outputs/20240123_124153_new'
#     # 用来存储改性后的材料cif文件
#     result_path = os.path.join(mofcof_path,"outputs_modified")

#     modify_MOF_COF(result_path, functions_path, mofcof_path)