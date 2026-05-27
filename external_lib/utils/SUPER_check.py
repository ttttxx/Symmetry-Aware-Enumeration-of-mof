import os
import re
import numpy as np
import logging
import logging.handlers

pi = np.pi
PERIOD_MIN = 0.90
PERIOD_MAX = 1.0 - PERIOD_MIN
BOND_MAX = 4.0
BOND_MIN = 0.8
BOND_LENGTH={'CC':1.54,'CN':1.48,'CH':1.11,'CB':1.56,'CO':1.43,'NB':1.45,'BB':1.62,'NN':1.46,'BO':1.38,'NO':1.42}

vertices = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
      'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
      'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
      'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr',
      'Ra', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th',
      'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'FG', 'X']


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

def process_cif(cif_path):
    
    with open(cif_path, 'r') as f:
        lines = f.readlines()
        # 创建节点坐标字典
    node_dict = {}
    isPass = True
    # f.close()

    print("母体材料的文件长度为：",len(lines))

    for line in lines:
        l = line.split()
        if '_cell_length_a' in line:
            aL = l[1]
        if '_cell_length_b' in line:
            bL = l[1]
        if '_cell_length_c' in line:
            cL = l[1]
        if '_cell_angle_alpha' in line:
            alpha = l[1]
        if '_cell_angle_beta' in line:
            beta = l[1]
        if '_cell_angle_gamma' in line:
            gamma = l[1]
        # flag = is_atom_or_bond(line)
        if isvert(l):
            node_dict[l[0]] = [l[2], l[3], l[4]]

    aL, bL, cL, alpha, beta, gamma = list(map(float, (aL, bL, cL, alpha, beta, gamma)))
    # 根据几何关系，计算出晶胞向量坐标
    ax = aL
    ay = 0.0
    az = 0.0
    bx = bL * np.cos(gamma * pi / 180.0)
    by = bL * np.sin(gamma * pi / 180.0)
    bz = 0.0
    cx = cL * np.cos(beta * pi / 180.0)
    cy = (cL * bL * np.cos(alpha * pi / 180.0) - bx * cx) / by
    cz = (cL ** 2.0 - cx ** 2.0 - cy ** 2.0) ** 0.5
    unit_cell = np.asarray([[ax, ay, az], [bx, by, bz], [cx, cy, cz]]).T

    for line in lines:
        l = line.split()
        if isedge(l):
            # 原子坐标合法化判断
            flag0 = True
            element0 = nn(l[0])
            element1 = nn(l[1])
            s0 = element0 + element1
            s1 = element1 + element0
            for value0, value1 in zip(node_dict[l[0]], node_dict[l[1]]):
                if float(value0) < PERIOD_MIN or float(value0) > PERIOD_MAX or float(value1) < PERIOD_MIN or float(value1) > PERIOD_MAX:
                    flag0 = False
                    break
            if flag0:
                f_nvec0 = np.asarray(list(map(float, node_dict[l[0]])))
                f_nvec1 = np.asarray(list(map(float, node_dict[l[1]])))
                c_nvec0 = np.dot(unit_cell, f_nvec0)
                c_nvec1 = np.dot(unit_cell, f_nvec1)
                dist = np.linalg.norm(c_nvec0 - c_nvec1)
                if s0 in BOND_LENGTH or s1 in BOND_LENGTH:
                    bond_min=BOND_LENGTH[s0]*0.85 if s0 in BOND_LENGTH else BOND_LENGTH[s1]*0.85
                    bond_max=BOND_LENGTH[s0]*1.13 if s0 in BOND_LENGTH else BOND_LENGTH[s1]*1.13
                    # bond_standard=BOND_LENGTH[s0] if s0 in BOND_LENGTH else BOND_LENGTH[s1]
                else:
                    bond_min=BOND_MIN
                    bond_max=BOND_MAX
                # if dist > bond_max or dist < bond_min or float(l[2])>bond_standard*1.05 or float(l[2])<bond_standard*0.9:
                if dist > bond_max or dist < bond_min:
                    isPass = False
                    break
            else:
                bond_len = float(l[2])
                if s0 in BOND_LENGTH or s1 in BOND_LENGTH:
                    bond_min=BOND_LENGTH[s0]*0.85 if s0 in BOND_LENGTH else BOND_LENGTH[s1]*0.85
                    bond_max=BOND_LENGTH[s0]*1.13 if s0 in BOND_LENGTH else BOND_LENGTH[s1]*1.13
                    # bond_standard=BOND_LENGTH[s0] if s0 in BOND_LENGTH else BOND_LENGTH[s1]
                else:
                    bond_min=BOND_MIN
                    bond_max=BOND_MAX
                if bond_len > bond_max or bond_len < bond_min:
                    isPass = False
                    break

    if isPass != True:
        logging.error('atom distance check failed!')
        print("atom distance check failed!")
        os.remove(cif_path)
    
    return isPass
'''
if __name__ == "__main__":

    cifs = os.listdir('../../output_cifs')

    print(len(cifs))

    for cif in cifs:

        process_cif('../../output_cifs/'+cif)

    print(len(cifs))

'''

