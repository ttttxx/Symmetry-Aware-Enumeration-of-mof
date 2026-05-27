import itertools
import os
import re

import numpy as np

from external_lib.utils.crystal_structure_optimization import adjust_coordinates_periodic_boundary
from external_lib.utils.parameter import judge_edge, judge_vertex, fractional_to_cartesian_matrix


class Edge():

    def __init__(self, edge_path):
        self.edge_path = edge_path
        self.edge_name = edge_path.split('/')[-1][:-len('.cif')]


    def generate_file_combinations(self, combination_length):
        # 这个函数的作用是生成给定长度的文件组合。它首先从指定的`path`（在这个例子中是`edges`）
        # 中列出所有文件，然后使用`itertools.product`生成所有可能的组合，重复次数等于`combination_length`的长度。
        # 最后返回这些文件组合列表。
        '''
        path = 'edges'
        edge_file_combinaitons = list(itertools.product([e for e in os.listdir(path)], repeat=len(combination_length)))

        return edge_file_combinaitons
        '''
        path = 'edges'
        edge_file_combinaitons = list(itertools.product([e for e in os.listdir(path)], repeat=len(combination_length)))

        return edge_file_combinaitons
    
    def calculate_unitcell_fcoords(self):
        """
        function: reads a CIF file from the file path specified in edge_path, extracts the crystal cell parameters and the fractional coordinates of atoms with the "X" type label, and returns the unit cell matrix and the fractional coordinates of the "X" atoms
        :return:  unit_cell, a 3x3 numpy array representing the unit cell matrix, and fcoords, a list of tuples where each tuple contains a label and a 3-element numpy array representing the fractional
        :parameter: self.edge_path: the path to the CIF file containing the edge information.
        """

        path = self.edge_path
        with open(path, 'r') as cif:
            cif = cif.read()
            cif = filter(None, cif.split('\n'))

        fractional_coords = []
        add_fractional_coords = fractional_coords.append
        for line in cif:
            if '_modify_baseline_atom_site_label' in line:
                break
            s = line.split()
            if '_cell_length_a' in line:
                a = s[1]
            if '_cell_length_b' in line:
                b = s[1]
            if '_cell_length_c' in line:
                c = s[1]
            if '_cell_angle_alpha' in line:
                alpha = s[1]
            if '_cell_angle_beta' in line:
                beta = s[1]
            if '_cell_angle_gamma' in line:
                gamma = s[1]
            if judge_vertex(s) and 'X' in s[0]:
                fractionsl_vector = np.array([float(q) for q in s[2:5]])
                add_fractional_coords([s[0], fractionsl_vector])


        a, b, c, alpha, beta, gamma = list(map(float, (a, b, c, alpha, beta, gamma)))
        unit_cell = fractional_to_cartesian_matrix(a, b, c, alpha, beta, gamma)

        return unit_cell,fractional_coords

    def calculate_XX_edgedistance(self):
        """
        function: calculates the Cartesian coordinates of each vertex of type 'X' in the crystal structure and then calculates the Euclidean distance between the Cartesian coordinates of the first two vertices of type 'X' in the crystal structure.
        :return:  returns a single floating-point value, which represents the Euclidean distance between the Cartesian coordinates of the first two vertices of type 'X' in the crystal structure.
        :parameter: edge_path: the path to the cif file of the crystal structure
        """
        
        path = self.edge_path
        with open(path, 'r') as cif:
            cif = cif.read()
            cif = filter(None, cif.split('\n'))

        fractional_coords = []
        add_fractional_coord = fractional_coords.append
        for line in cif:
            if '_modify_baseline_atom_site_label' in line:
                break
            s = line.split()
            if '_cell_length_a' in line:
                a = s[1]
            if '_cell_length_b' in line:
                b = s[1]
            if '_cell_length_c' in line:
                c = s[1]
            if '_cell_angle_alpha' in line:
                alpha = s[1]
            if '_cell_angle_beta' in line:
                beta = s[1]
            if '_cell_angle_gamma' in line:
                gamma = s[1]
            if judge_vertex(s) and 'X' in s[0]:
                fvec = np.array([float(q) for q in s[2:5]])
                add_fractional_coord([s[0], fvec])


        a, b, c, alpha, beta, gamma = list(map(float, (a, b, c, alpha, beta, gamma)))
        unit_cell = fractional_to_cartesian_matrix(a, b, c, alpha, beta, gamma)

        wrapped_fractional_coords = [adjust_coordinates_periodic_boundary(fractional_coords[0][1], vec[1]) for vec in fractional_coords]  # 求第一个X类型的节点与所有X类型节点的PBC3DF
        cartesian_coords = [np.dot(unit_cell, vec) for vec in wrapped_fractional_coords]

        return np.linalg.norm(cartesian_coords[0] - cartesian_coords[1])

    def calculate_cartesian_coords(self, label):

        """
        function: calculates the Cartesian coordinates of each vertex in a crystal structure using the fractional coordinates of the vertices and the lattice vectors of the crystal. It returns a list of the Cartesian coordinates, centered around the center of mass of all the vertices if the label input is true.
        :param : label: boolean, if true, the coordinates are centered around the center of mass of all the vertices
        :return: a list of the Cartesian coordinates of the vertices
        """

        unit_cell, fcoords = self.calculate_unitcell_fcoords()

        wrapped_fractional_coords = [[vec[0], adjust_coordinates_periodic_boundary(fcoords[0][1], vec[1])] for vec in fcoords]

        if label:
            cartesian_coords = [[vec[0], np.dot(unit_cell, vec[1])] for vec in wrapped_fractional_coords]
            center_of_mass = np.average(np.asarray([vec[1] for vec in cartesian_coords]), axis=0)
            shifted_ccoords = [[vec[0], vec[1] - center_of_mass] for vec in cartesian_coords]
        else:
            cartesian_coords = [np.dot(unit_cell, vec[1]) for vec in wrapped_fractional_coords]
            center_of_mass = np.average(cartesian_coords, axis=0)
            shifted_ccoords = [vec - center_of_mass for vec in cartesian_coords]

        return shifted_ccoords

    def crystal_to_array(self):
        """
        这个函数 crystal_to_array 是将一个晶体结构（从CIF文件读取）转换为一个坐标数组。它首先打开并读取CIF文件中的内容，
        然后提取晶体的参数（a，b，c，alpha，beta，gamma）。然后，它遍历CIF文件的每一行，如果这一行包含一个顶点（通过函数 judge_vertex 判断），
        则提取并保存其分数坐标。然后，根据晶体的参数，计算出晶格向量，并将其存储在一个矩阵中。然后，它通过将分数坐标乘以晶格向量，将分数坐标转换为笛卡尔坐标。
        最后，它计算出所有节点类型为X的坐标的平均值（即质心），并将这个值从所有的笛卡尔坐标中减去，以将质心移动到原点。最终返回的是一个包含了节点标签和质心为原点的笛卡尔坐标的列表。
        :return:
        """
        path = self.edge_path

        with open(path, 'r') as cif:
            cif = cif.read()
            cif = filter(None, cif.split('\n'))

        fractional_coords = []
        add_fractional_coord = fractional_coords.append
        modified_point = []

        for line in cif:
            if '_modify_baseline_atom_site_label' in line:
                break
            s = line.split()

            if '_cell_length_a' in line:
                a = s[1]
            if '_cell_length_b' in line:
                b = s[1]
            if '_cell_length_c' in line:
                c = s[1]
            if '_cell_angle_alpha' in line:
                alpha = s[1]
            if '_cell_angle_beta' in line:
                beta = s[1]
            if '_cell_angle_gamma' in line:
                gamma = s[1]
            '''
                        if '_modified_remove_atom_site_label' in line:
                flag = True
                continue
            '''
            if judge_vertex(s):
                fvec = np.array([float(q) for q in s[2:5]])
                add_fractional_coord([s[0], fvec])

        a, b, c, alpha, beta, gamma = list(map(float, (a, b, c, alpha, beta, gamma)))
        lattice_vectors = fractional_to_cartesian_matrix(a, b, c, alpha, beta, gamma)

        reference_vector = fractional_coords[0][1]

        cartesian_coords = [[n[0], np.dot(lattice_vectors, adjust_coordinates_periodic_boundary(reference_vector, n[1]))] for n in fractional_coords]

        center_of_mass = np.average(np.array([n[1] for n in cartesian_coords if re.sub('[0-9]', '', n[0]) == 'X']), axis=0)
        centered_cart_coords = [[n[0], n[1] - center_of_mass] for n in cartesian_coords]

        return centered_cart_coords, modified_point

    def extract_edge_bonds(self):
        """

        """
        path = self.edge_path

        with open(path, 'r') as cif:
            cif = cif.read()
            cif = filter(None, cif.split('\n'))

        bond_info = []
        add_edge_info = bond_info.append
        for line in cif:
            if '_modify_baseline_atom_site_label' in line:
                break
            s = line.split()
            if judge_edge(s):
                add_edge_info(s)

        return bond_info

    def extract_edge_info(self):
        '''

        '''
        path = self.edge_path

        with open(path, 'r') as cif:
            cif = cif.read()
            cif = filter(None, cif.split('\n'))

        vertex_charges = []
        charges_append = vertex_charges.append
        vertex_elements = []
        elements_append = vertex_elements.append
        for line in cif:
            if '_modify_baseline_atom_site_label' in line:
                break
            line_token = line.split()
            if judge_vertex(line_token):
                charges_append(line_token[-1])
                elements_append(line_token[1])

        return vertex_charges, vertex_elements

    def identify_modified_points(self):

        path = self.edge_path

        with open(path, 'r') as cif:
            cif = cif.read()
            cif = filter(None, cif.split('\n'))

        modified_points = []
        remove_points = []
        flag = False

        for line in cif:
            s = line.split()
            if "_modified_remove_atom_site_label"  in line:
                flag = True
                continue
            if flag:
                modified_points.append(s[0])
                if len(s) > 1:
                    remove_points.append(s[1])
                continue

        return modified_points, remove_points
