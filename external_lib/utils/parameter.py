import re

import numpy as np

PARALLEl = True

# MOF_MODE = False
# COF_MODE = True
# MODIFIED_MODE = False

ONLY_MOF = False
SINGAL_METAL = False
IGNORE_WRONG = False
DEBUG = True
MAX_LENGTH = 1.0e6
CHARGES = True

metals = ['Mg', 'Al', 'Ca', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Sr', 'Y', 'Zr', 'Nb', 'Mo',
          'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl',
          'Pb', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'U']

vertices = ('V', 'Er', 'Ti', 'Ce', 'S',
            'H', 'He', 'Li', 'Be', 'B',
            'C', 'N', 'O', 'F', 'Ne',
            'Na', 'Mg', 'Al', 'Si', 'P',
            'Cl', 'Ar', 'K', 'Ca', 'Sc',
            'Cr', 'Mn', 'Fe', 'Co', 'Ni')

vertex_elements = ('V', 'Er', 'Ti', 'Ce', 'S',
                   'H', 'He', 'Li', 'Be', 'B',
                   'C', 'N', 'O', 'F', 'Ne',
                   'Na', 'Mg', 'Al', 'Si', 'P',
                   'Cl', 'Ar', 'K', 'Ca', 'Sc',
                   'Cr', 'Mn', 'Fe', 'Co', 'Ni')

vname_dict = {'V': 1, 'Er': 2, 'Ti': 3, 'Ce': 4, 'S': 5,
              'H': 6, 'He': 7, 'Li': 8, 'Be': 9, 'B': 10,
              'C': 11, 'N': 12, 'O': 13, 'F': 14, 'Ne': 15,
              'Na': 16, 'Mg': 17, 'Al': 18, 'Si': 19, 'P': 20,
              'Cl': 21, 'Ar': 22, 'K': 23, 'Ca': 24, 'Sc': 24,
              'Cr': 26, 'Mn': 27, 'Fe': 28, 'Co': 29, 'Ni': 30}

bond_type = ('S', 'D', 'T', 'A')

metal_elements = ['Ac', 'Ag', 'Al', 'Am', 'Au', 'Ba', 'Be', 'Bi',
                  'Bk', 'Ca', 'Cd', 'Ce', 'Cf', 'Cm', 'Co', 'Cr',
                  'Cs', 'Cu', 'Dy', 'Er', 'Es', 'Eu', 'Fe', 'Fm',
                  'Ga', 'Gd', 'Hf', 'Hg', 'Ho', 'In', 'Ir',
                  'K', 'La', 'Li', 'Lr', 'Lu', 'Md', 'Mg', 'Mn',
                  'Mo', 'Na', 'Nb', 'Nd', 'Ni', 'No', 'Np', 'Os',
                  'Pa', 'Pb', 'Pd', 'Pm', 'Pr', 'Pt', 'Pu', 'Ra',
                  'Rb', 'Re', 'Rh', 'Ru', 'Sc', 'Sm', 'Sn', 'Sr',
                  'Ta', 'Tb', 'Tc', 'Th', 'Ti', 'Tl', 'Tm', 'U',
                  'V', 'W', 'Y', 'Yb', 'Zn', 'Zr']

PT = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
      'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
      'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
      'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr',
      'Ra', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th',
      'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'FG', 'X']

class ShapeMismatchError(ValueError):
    """配体连接点数量与拓扑锚点不匹配时抛出。"""
    pass


def fractional_to_cartesian_matrix(a, b, c, alpha, beta, gamma):
    """将晶胞参数转换为笛卡尔变换矩阵 (3x3)。

    Returns:
        3x3矩阵 M, 满足 Cartesian = M @ Fractional。
        M 的列向量为晶格基矢 a, b, c 在笛卡尔坐标系中的表示。
    """
    pi = np.pi
    ax = a
    ay = 0.0
    az = 0.0
    bx = b * np.cos(gamma * pi / 180.0)
    by = b * np.sin(gamma * pi / 180.0)
    bz = 0.0
    cx = c * np.cos(beta * pi / 180.0)
    cy = (c * b * np.cos(alpha * pi / 180.0) - bx * cx) / by
    cz = (c ** 2.0 - cx ** 2.0 - cy ** 2.0) ** 0.5
    return np.asarray([[ax, ay, az], [bx, by, bz], [cx, cy, cz]]).T


def isfloat(value):
    """
    function: judge whether a string is a float number
    :param value: string
    :return: True or False
    """
    try:
        float(value)
        return True
    except ValueError:
        return False


def judge_vertex(line):

    if nn(line[0]) in PT and line[1] in PT and False not in map(isfloat, line[2:5]):

        return True
    else:
        return False


def judge_edge(line):

    if nn(line[0]) in PT and nn(line[1]) in PT and isfloat(line[2]) and line[-1] in ('S', 'D', 'T', 'A'):

        return True
    else:
        return False


def is_metalvertex(line):


    if len(line) >= 5:
        line_state = map(isfloat, line[2:5])
        replaced_line0 = re.sub('[^a-zA-Z]', '', line[0])
        if replaced_line0 in vertices and line[1] in vertices and False not in line_state:

            return True
        else:
            return False


def is_metaledge(line):
    #
    if len(line) >= 5:
        replaced_line0 = re.sub('[^a-zA-Z]', '', line[0])
        replaced_line1 = re.sub('[^a-zA-Z]', '', line[1])
        line2_state = isfloat(line[2])
        if replaced_line0 in vertices and replaced_line1 in vertices and line2_state and line[-1] in (
                'S', 'D', 'T', 'A'):

            return True
        else:
            return False

def nn(string):
    return re.sub('[^a-zA-Z]', '', string)