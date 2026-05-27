import re

import numpy as np

from external_lib.utils.crystal_structure_optimization import adjust_coordinates_periodic_boundary
from external_lib.utils.parameter import judge_vertex, judge_edge, fractional_to_cartesian_matrix


class Node():

    def __init__(self, node_path):
        self.node_path = node_path
        self.node_name = node_path.split('/')[-1][:-len('.cif')]

    def identify_node_connection_points(self):
        '''
        function:  reads a CIF file and identifies the number of connection points (nodes) that have a label containing the letter 'X'.
        :return: the number of nodes in the node file that have a label containing the letter 'X'.
        :parameter: None
        '''
        path = self.node_path
        with open(path, 'r') as cif:
            cif = cif.read()
            cif = filter(None, cif.split('\n'))

        num_x_nodes = 0
        for line in cif:
            line_tokens = line.split()
            if judge_vertex(line_tokens):
                if re.sub('[^a-zA-Z]', '', line_tokens[0]) == 'X':
                    num_x_nodes += 1
        return num_x_nodes

    def calculate_X_cart_vectors(self, label):
        '''
        function: extracts coordinate information for nodes labeled 'X' then calculates the Cartesian coordinates for each node based on the unit cell dimensions provided in the file.
        :return: a list of Cartesian coordinates for each labeled node in the node file, with the center of mass shifted to the origin if label is True, or a list of Cartesian coordinates for all nodes in the node file with the center of mass shifted to the origin if label is False.
        :parameter: label: if True, the coordinates of the nodes are shifted to the center of mass of the unit cell.
        '''
        node_file_path = self.node_path
        with open(node_file_path, 'r') as cif:
            cif = cif.read()
            cif = filter(None, cif.split('\n'))

        labeled_node_coords = []
        labeled_node_coords_append = labeled_node_coords.append

        for line in cif:  # 获取文件中的相关信息
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
            if judge_vertex(s) and 'X' in s[0]:  # 如果存在某条边的结点类型包含X，则提取坐标
                coord_vector = np.array([float(q) for q in s[2:5]])
                labeled_node_coords_append([s[0], coord_vector])

        a, b, c, alpha, beta, gamma = list(map(float, (a, b, c, alpha, beta, gamma)))
        unit_cell = fractional_to_cartesian_matrix(a, b, c, alpha, beta, gamma)

        # 计算最小分数坐标
        mic_labeled_node_coords = [[vec[0], adjust_coordinates_periodic_boundary(labeled_node_coords[0][1], vec[1])] for vec in labeled_node_coords]

        # 计算笛卡尔坐标
        if label:
            cart_coords = [[vec[0], np.dot(unit_cell, vec[1])] for vec in mic_labeled_node_coords]
            center_of_mass = np.average(np.asarray([vec[1] for vec in cart_coords]), axis=0)
            shifted_cart_coords = [[vec[0], vec[1] - center_of_mass] for vec in cart_coords]
        else:
            cart_coords = [np.dot(unit_cell, vec[1]) for vec in mic_labeled_node_coords]
            center_of_mass = np.average(cart_coords, axis=0)
            shifted_cart_coords = [vec - center_of_mass for vec in cart_coords]

        return shifted_cart_coords

    def get_all_elements(self):
        '''
        function: extracts all elements from the cif file
        :return: a list of all elements in the cif file
        :parameter: the path to the cif file
        '''
        path = self.node_path
        with open(path, 'r') as cif:
            cif = cif.read()
            cif = filter(None, cif.split('\n'))

        elems = []
        elems_append = elems.append
        for line in cif:
            s = line.split()

            if judge_vertex(s):
                elems_append(s[1])

        return elems

    def calculate_periodic_boundary_conditions(self):
        '''
        function: extracts information about atomic positions then calculates the periodic boundary conditions and centers the positions around the center of mass of all atoms with element type X
        :return: a list of element symbols and their corresponding positions that have been shifted to be centered around the center of mass of all atoms with element type 'X'
        :parameter: the path to the cif file
        '''
        node_file_path = self.node_path

        with open(node_file_path, 'r') as cif:
            cif = cif.read()
            cif = filter(None, cif.split('\n'))

        atom_positions = []
        atom_positions_append = atom_positions.append
        for line in cif:
            line_contents = line.split()
            if '_cell_length_a' in line:
                a = line_contents[1]
            if '_cell_length_b' in line:
                b = line_contents[1]
            if '_cell_length_c' in line:
                c = line_contents[1]
            if '_cell_angle_alpha' in line:
                alpha = line_contents[1]
            if '_cell_angle_beta' in line:
                beta = line_contents[1]
            if '_cell_angle_gamma' in line:
                gamma = line_contents[1]
            if judge_vertex(line_contents):
                fvec = np.array([float(q) for q in line_contents[2:5]])
                atom_positions_append([line_contents[0], fvec])

        a, b, c, alpha, beta, gamma = list(map(float, (a, b, c, alpha, beta, gamma)))
        unit_cell = fractional_to_cartesian_matrix(a, b, c, alpha, beta, gamma)

        cell_normal_vector = atom_positions[0][1]
        crystal_coords = [[n[0], np.dot(unit_cell, adjust_coordinates_periodic_boundary(cell_normal_vector, n[1]))] for n in atom_positions]
        center_of_mass = np.average(np.array([n[1] for n in crystal_coords if re.sub('[0-9]', '', n[0]) == 'X']), axis=0)
        shifted_atom_positions = [[n[0], n[1] - center_of_mass] for n in crystal_coords]

        return shifted_atom_positions

    def get_node_info(self):
        '''
        function: read a file and extracts the bonds of each atom, it then returns a list of bonds
        :return: a list of edge information extracted from the file
        :parameter: path: the path of the node file
        '''
        node_file_path = self.node_path

        with open(node_file_path, 'r') as cif:
            cif = cif.read()
            cif = filter(None, cif.split('\n'))

        edge_bond_info = []
        bonds_append = edge_bond_info.append
        for line in cif:
            line_contents = line.split()
            if judge_edge(line_contents):
                bonds_append(line_contents)

        return edge_bond_info

    def get_node_charges(self):
        '''
        function: read a file and extracts the charges of each atom and the elements of each atom, it then returns a list of charges and a list of elements
        :return: two lists, one is the charges of each atom, the other is the elements of each atom
        :parameter: path: the path of the node file
        '''
        node_file_path = self.node_path

        with open(node_file_path, 'r') as cif:
            cif = cif.read()
            cif = filter(None, cif.split('\n'))

        node_charges = []
        charges_append = node_charges.append
        node_elements = []
        elements_append = node_elements.append
        for line in cif:
            line_contents = line.split()
            if judge_vertex(line_contents):
                charges_append(line_contents[-1])
                elements_append(line_contents[1])

        return node_charges, node_elements
