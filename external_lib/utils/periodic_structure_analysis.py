import logging
import re
from itertools import combinations
from random import choice
import math
import networkx as nx
import numpy as np
import os
from external_lib.utils.crystal_structure_optimization import (adjust_coordinates_periodic_boundary)
from external_lib.utils.crystal_structure_optimization import minimize_rms_superimpose
from external_lib.utils.parameter import metals, fractional_to_cartesian_matrix
from external_lib.edge import Edge
from external_lib.node import Node

def adjust_placed_edges_to_fit_periodicity(placed_edges, placed_nodes, sc_unit_cell):
    adjusted_placed_edges = []
    adjusted_placed_edges_extend = adjusted_placed_edges.extend

    placed_edges = np.asarray(placed_edges)
    edge_labels = set(map(int, placed_edges[:, -1]))

    edge_dict = dict((k, []) for k in edge_labels)

    node_connection_points = [list(map(float, i[1:4])) for i in placed_nodes if re.sub('[0-9]', '', i[5]) == 'X']

    for edge in placed_edges:
        ty = int(edge[-1])
        edge_dict[ty].append(edge)

    for k in edge_dict:

        edge = np.asarray(edge_dict[k])
        elems = edge[:, 0]
        evecs = [list(map(float, i)) for i in edge[:, 1:4]]
        charges = edge[:, 4]
        cp = edge[:, 5]
        ty = edge[:, 6]

        xvecs = [list(map(float, i)) for (i, j) in zip(evecs, cp) if re.sub('[0-9]', '', j) == 'X']
        relavent_node_xvecs = []
        relavent_node_xvecs_append = relavent_node_xvecs.append

        for ex in xvecs:
            min_dist = (1e6, [], 0)

            if len(node_connection_points) == 0:
                relavent_node_xvecs_append(ex)
                continue

            f_ex = np.dot(np.linalg.inv(sc_unit_cell), ex)
            for i in range(len(node_connection_points)):
                nx = node_connection_points[i]
                f_nx = np.dot(np.linalg.inv(sc_unit_cell), nx)

                fdist_vec, sym = compute_PBC_distance_and_symmetry(f_ex, f_nx)
                cdist = np.linalg.norm(np.dot(sc_unit_cell, fdist_vec))

                if cdist < min_dist[0]:
                    min_dist = (cdist, np.dot(sc_unit_cell, f_nx + sym), i)

            node_connection_points.pop(min_dist[2])
            relavent_node_xvecs_append(min_dist[1])

        ecom = np.average(xvecs, axis=0)
        rnxcom = np.average(relavent_node_xvecs, axis=0)

        evecs = np.asarray(evecs - ecom)
        xvecs = np.asarray(xvecs - ecom)

        relavent_node_xvecs = np.asarray(relavent_node_xvecs - rnxcom)

        trans = rnxcom

        min_dist, rot, tran = minimize_rms_superimpose(xvecs, relavent_node_xvecs)
        adjusted_evecs = np.dot(evecs, rot) + trans

        adjusted_edge = np.column_stack((elems, adjusted_evecs, charges, cp, ty))
        adjusted_placed_edges_extend(adjusted_edge)

    return adjusted_placed_edges

def compute_PBC_distance_and_symmetry(vec1, vec2):

    dX, dY, dZ = vec1 - vec2
    # ndX, ndY, ndZ = 0.0, 0.0, 0.0
    if dX > 0.5:
        s1 = 1.0
        ndX = dX - 1.0
    elif dX < -0.5:
        s1 = -1.0
        ndX = dX + 1.0
    else:
        s1 = 0.0
        ndX = dX

    if dY > 0.5:
        s2 = 1.0
        ndY = dY - 1.0
    elif dY < -0.5:
        s2 = -1.0
        ndY = dY + 1.0
    else:
        s2 = 0.0
        ndY = dY

    if dZ > 0.5:
        s3 = 1.0
        ndZ = dZ - 1.0
    elif dZ < -0.5:
        s3 = -1.0
        ndZ = dZ + 1.0
    else:
        s3 = 0.0
        ndZ = dZ

    sym = np.array([s1, s2, s3])

    return np.array([ndX, ndY, ndZ]), sym

def compute_scaled_node_and_edge_vectors(sc_coords, sc_omega_plus, sc_unit_cell, ea_dict):

    nvecs = []
    evecs = []
    already_placed_edges = []
    nvecs_append = nvecs.append
    evecs_append = evecs.append
    already_placed_edges_append = already_placed_edges.append

    for n in sc_coords:

        vertex, vcif, vfvec, indicent_edges = n
        vcvec = np.dot(sc_unit_cell, vfvec)

        ie = []
        ie_append = ie.append

        for e in indicent_edges:

            ind = e[0]  # index
            positive_direction = e[1]
            ecif = e[2]

            if vertex == positive_direction[0]:
                direction = 1
                on = positive_direction[1]
            else:
                direction = -1
                on = positive_direction[0]

            dxn = ea_dict[vertex][ind][1]
            dxon = ea_dict[on][ind][1]

            ie_append((ind, direction, ecif, dxn, dxon))

        efvec = []
        ecvec = []
        efvec_append = efvec.append
        ecvec_append = ecvec.append

        for e in ie:
            ind, d, ecif, dxn, dxon = e
            cs = np.dot(sc_unit_cell, sc_omega_plus[ind - 1])
            fvec = vfvec + d * sc_omega_plus[ind - 1]
            cvec = vcvec + d * cs

            cs_norm = np.linalg.norm(cs) + 1e-12
            ec1 = vcvec + d * dxn * (cs / cs_norm)
            ec2 = cvec - d * dxon * (cs / cs_norm)

            ecoords = np.average([ec1, ec2], axis=0)

            ecvec_append(cvec)
            efvec_append(fvec)

            if ind not in already_placed_edges:

                evecs_append((ind, ecif, ecoords, np.array([ec1, ec2])))
                already_placed_edges_append(ind)

        nvecs_append((vertex, vcvec, vcif, np.asarray(ecvec)))

    return nvecs, evecs

def adjust_node_charges_to_neutral(placed_all):

    placed_all = np.asarray(placed_all)
    net_charge = np.sum([float(i) for i in placed_all[:, 4]])

    dont_change_charge_atoms = []
    dont_change_charge_atoms_append = dont_change_charge_atoms.append

    for l in placed_all:

        elem = re.sub('[0-9]', '', l[0])

        if elem in metals:
            dont_change_charge_atoms_append(l[0])

        elif 'FG' in elem or 'fg' in elem or 'fG' in elem or 'Fg' in elem:
            dont_change_charge_atoms_append(l[0])

    num_atoms = float(len(placed_all) - len(dont_change_charge_atoms))
    rcb = net_charge / num_atoms

    fc_placed_all = []
    fc_placed_all_append = fc_placed_all.append

    for l in placed_all:

        e, x, y, z, c, oe, i, bbtype = l
        c = float(c)

        if e in dont_change_charge_atoms:

            fc_placed_all_append([e, x, y, z, c, oe, i, bbtype])

        else:

            nc = c - rcb
            fc_placed_all_append([e, x, y, z, nc, oe, i, bbtype])

    for l in fc_placed_all:

        mode = choice(['up', 'down'])
        if mode == 'down':
            l[4] = round_down_to_decimal_places(l[4], 4)
        else:
            l[4] = round_to_decimal_places(l[4], 4)

    nnet_charge = np.sum([l[4] for l in fc_placed_all])

    return fc_placed_all, nnet_charge, net_charge, rcb

def calculate_and_fix_bond_symmetry(bonds_all, placed_all, sc_unit_cell):

    coords_dict = dict((l[0], np.dot(np.linalg.inv(sc_unit_cell), list(map(float, l[1:4])))) for l in placed_all)

    fixed_bonds = []
    fixed_bonds_append = fixed_bonds.append
    for l in bonds_all:
        vec1 = coords_dict[l[0]]
        vec2 = coords_dict[l[1]]

        dist, sym = calculate_periodic_distance_and_symmetry(vec1, vec2)

        fixed_bonds_append([l[0], l[1], l[2], sym, l[4]])

    return fixed_bonds

def build_graph_and_validate_bonds(placed_all, bonds_all, sc_unit_cell, max_length, bond_tol, nconnections_list,
                                   num_possible_XX_bonds):

    G = nx.Graph()

    for n in placed_all:
        isX = (re.sub('[0-9]', '', n[-3]) == 'X')

        G.add_node(n[0], coords=np.array(list(map(float, (n[1:4])))), occ=n[4], isX=isX, bbcode=int(n[6]), nconnect=1,
                   bbtype=n[-1])
    for l in bonds_all:
        G.add_edge(l[0], l[1], length=l[2], sym=l[3], ty=l[4], is_new_XX=False, order=(l[0], l[1]))

    for line in nconnections_list:
        node, nc = line
        G.nodes[node]['nconnect'] -= 1
        G.nodes[node]['nconnect'] += nc

    previous_degrees = dict((node, G.degree(node)) for node in G.nodes())

    ccs = list(nx.connected_components(G))
    count = 0
    bb_tol = max_length + 0.50 * max_length
    logging.debug('distance search tolerance is %s Angstroms', np.round(bb_tol, 3))
    for connect_comp0, connect_comp1 in combinations(ccs, 2):

        xname0 = [n for n in connect_comp0 if G.nodes[n]['isX']]
        xvecs0 = [np.dot(np.linalg.inv(sc_unit_cell), G.nodes[n]['coords']) for n in connect_comp0 if
                  G.nodes[n]['isX']]
        com0 = np.average(xvecs0, axis=0)
        type0 = list(set([G.nodes[n]['bbtype'] for n in connect_comp0]))

        xname1 = [n for n in connect_comp1 if G.nodes[n]['isX']]
        xvecs1 = [np.dot(np.linalg.inv(sc_unit_cell), G.nodes[n]['coords']) for n in connect_comp1 if G.nodes[n]['isX']]
        com1 = np.average(xvecs1, axis=0)
        type1 = list(set([G.nodes[n]['bbtype'] for n in connect_comp1]))

        if len(type0) > 1 or len(type1) > 1:
            logging.error('%s %s', type0, type1)
            logging.error('building block indicated as both node and edge type')
            raise ValueError('building block indicated as both node and edge type')

        if len(xname0) == 0 or len(xname1) == 0:
            logging.error('There are connected components with no connection site atoms')
            raise ValueError('There are connected components with no connection site atoms')
        type0 = type0[0]
        type1 = type1[0]

        if type0 == 'edge' and type1 == 'edge':
            continue

        com_dist = np.linalg.norm(np.dot(sc_unit_cell, com0 - adjust_coordinates_periodic_boundary(com0, com1)))

        if com_dist <= bb_tol:

            min_dist = (100.0, 'no.3', 'the', 'larch')

            for xv0, xn0 in zip(xvecs0, xname0):
                for xv1, xn1 in zip(xvecs1, xname1):

                    DV, sym = calculate_periodic_distance_and_symmetry(xv0, xv1)
                    dist = np.linalg.norm(np.dot(sc_unit_cell, DV))

                    if dist < min_dist[0]:
                        min_dist = (dist, xn0, xn1, sym)

            if min_dist[0] < bond_tol:
                count += 1
                G.add_edge(min_dist[1], min_dist[2], length=min_dist[0], sym=min_dist[3], ty='S', is_new_XX=True,
                           order=(min_dist[1], min_dist[2]))

    connection_nodes = [(n, data) for n, data in G.nodes(data=True) if data['isX']]
    no_bond_connection_sites = []

    for node, data in connection_nodes:

        bbcode0 = data['bbcode']
        nconnections = data['nconnect']
        nbors = list(G.neighbors(node))
        X_nbors = [nbor for nbor in nbors if
                   (G.nodes[nbor]['isX'] and G.nodes[nbor]['bbcode'] != bbcode0)]

        if len(X_nbors) < nconnections:
            no_bond_connection_sites.append(node)

    for n0, n1 in combinations(no_bond_connection_sites, 2):

        vec0 = np.dot(np.linalg.inv(sc_unit_cell), G.nodes[n0]['coords'])
        vec1 = np.dot(np.linalg.inv(sc_unit_cell), G.nodes[n1]['coords'])

        bbcode0 = G.nodes[n0]['bbcode']
        bbcode1 = G.nodes[n1]['bbcode']

        bbtype0 = G.nodes[n0]['bbtype']
        bbtype1 = G.nodes[n1]['bbtype']

        if bbcode0 == bbcode1:
            continue
        if bbtype0 == 'edge' and bbtype1 == 'edge':
            continue

        DV, sym = calculate_periodic_distance_and_symmetry(vec0, vec1)
        dist = np.linalg.norm(np.dot(sc_unit_cell, DV))

        if dist < bond_tol:
            count += 1
            G.add_edge(n0, n1, length=dist, sym=sym, ty='S', is_new_XX=True, order=(n0, n1))

    for node, data in connection_nodes:

        vec0 = np.dot(np.linalg.inv(sc_unit_cell), G.nodes[node]['coords'])
        bbcode0 = data['bbcode']
        bbtype0 = data['bbtype']

        nbors = list(G.neighbors(node))
        X_nbors = [nbor for nbor in nbors if (G.nodes[nbor]['isX'] and G.nodes[nbor]['bbcode'] != bbcode0)]
        nconnections = data['nconnect']

        if len(X_nbors) > nconnections:

            nbor_dists = []
            for nbor in X_nbors:
                vec1 = np.dot(np.linalg.inv(sc_unit_cell), G.nodes[nbor]['coords'])
                dist = np.linalg.norm(np.dot(sc_unit_cell, vec0 - adjust_coordinates_periodic_boundary(vec0, vec1)))
                nbor_dists.append((dist, nbor))
            nbor_dists.sort(key=lambda x: x[0])

            for nbor in nbor_dists[nconnections:]:
                G.remove_edge(node, nbor[1])
                count -= 1

    wrong_connection_nodes = []
    wrong_connection_nodes_append = wrong_connection_nodes.append

    for node, data in connection_nodes:

        nconnect = data['nconnect']
        nbors = list(G.neighbors(node))
        bbcode0 = data['bbcode']

        X_nbors = [nbor for nbor in nbors if (G.nodes[nbor]['isX'] and G.nodes[nbor]['bbcode'] != bbcode0)]
        if len(X_nbors) != nconnect:
            wrong_connection_nodes_append(node)

    XX_bond_count = 0
    for e0, e1, data in G.edges(data=True):
        if data['is_new_XX']:
            XX_bond_count += 1

    bond_check_passed = True

    for node, data in G.nodes(data=True):

        isX = data['isX']
        degree = G.degree(node)
        previous_degree = previous_degrees[node]
        nconnections = data['nconnect']

        if isX and degree != previous_degree + nconnections:
            logging.warning('There are connection sites with too many/not enough bonds formed.')
            bond_check_passed = False
            break
        elif not isX and degree != previous_degree:
            logging.warning('The degree of a non-X atom has changed after XX bond formation.')
            bond_check_passed = False
            break
    if len(wrong_connection_nodes) > 0:
        logging.warning('There are connection sites with too few or too many bonds')
        bond_check_passed = False
    if XX_bond_count < num_possible_XX_bonds:
        logging.warning('The number of XX bonds formed %d is less than the correct number %d', XX_bond_count,
                        num_possible_XX_bonds)
        logging.warning('Try increasing the bond formation distance tolerance.')
        bond_check_passed = False
    if XX_bond_count > num_possible_XX_bonds:
        logging.warning('The number of XX bonds formed %d is more than the correct number %d', XX_bond_count,
                        num_possible_XX_bonds)
        logging.warning('Check your inputs for incorrect X atoms.')
        bond_check_passed = False

    fixed_bonds = []
    fixed_bonds_append = fixed_bonds.append

    for edge in G.edges(data=True):

        edict = edge[2]
        ty = edict['ty']
        leng = edict['length']
        sy = edict['sym']
        order = edict['order']
        fixed_bonds_append([order[0], order[1], leng, sy, ty])

    return fixed_bonds, count, bond_check_passed

def find_and_calculate_bonds(placed_all, bonds_all, sc_unit_cell, tol):

    fixed_bonds = []
    used_bonds = []
    fixed_bonds_append = fixed_bonds.append
    used_bonds_append = used_bonds.append

    for l in bonds_all:
        fixed_bonds_append([l[0], l[1], l[2], '.', l[4]])
        used_bonds_append((l[0], l[1]))

    connection_points = [line for line in placed_all if re.sub('[0-9]', '', line[5]) == 'X']
    nbcount = 0

    for i in range(len(connection_points)):

        ielem = connection_points[i][0]
        ivec = np.dot(np.linalg.inv(sc_unit_cell), np.array([float(q) for q in connection_points[i][1:4]]))
        ibbid = int(connection_points[i][6])

        for j in range(i + 1, len(connection_points)):

            jelem = connection_points[j][0]
            jbbid = int(connection_points[j][6])

            if (ielem, jelem) not in used_bonds and (jelem, ielem) not in used_bonds:

                jvec = np.dot(np.linalg.inv(sc_unit_cell), list(map(float, connection_points[j][1:4])))
                DV, sym = calculate_periodic_distance_and_symmetry(ivec, jvec)
                dist = np.linalg.norm(np.dot(sc_unit_cell, DV))

                if dist < tol and ibbid != jbbid:
                    nbcount += 1
                    fixed_bonds_append([ielem, jelem, dist, sym, 'S'])
                    break
                else:
                    continue

    return fixed_bonds, nbcount

def place_and_calculate_nodes(nvecs):

    placed_nbb_coords = []
    placed_nbb_coords_extend = placed_nbb_coords.extend
    all_bonds = []
    all_bonds_extend = all_bonds.extend
    ind_seg = 0
    bbind = 1

    for n in nvecs:

        bbind = bbind + 1
        name, cvec, cif, nvec = n
        ll = 0

        for v in nvec:
            mag = np.linalg.norm(v - np.average(nvec, axis=0))
            if mag > ll:
                ll = mag

        temp_node_path = os.path.join('nodes', cif)
        temp_node = Node(node_path=temp_node_path)
        bbxvec = np.array(temp_node.calculate_X_cart_vectors(False))

        node_cent = np.mean(bbxvec, axis=0)
        bbxvec_rel = bbxvec - node_cent

        nvec_rel = nvec - cvec

        min_dist, rot, tran = minimize_rms_superimpose(bbxvec_rel, nvec_rel, center=False)

        all_bb = temp_node.calculate_periodic_boundary_conditions()
        all_coords = np.array([v[1] for v in all_bb])
        all_inds = np.array([v[0] for v in all_bb])
        chg, elem = temp_node.get_node_charges()
        all_names = [o + re.sub('[A-Za-z]', '', p) for o, p in
                     zip(elem, all_inds)]

        all_names_indices = np.array([int(re.sub('[A-Za-z]', '', e)) for e in all_names]) + ind_seg

        elem_dict = dict((k, '') for k in all_inds)
        for i, j in zip(all_inds, elem):
            elem_dict[i] = j

        ind_dict = dict((k, '') for k in all_inds)
        for i, j in zip(all_inds, all_names_indices):
            ind_dict[i] = j

        bonds = temp_node.get_node_info()

        anf = [str(elem_dict[n]) + str(ind_dict[n]) for n in all_inds]

        abf = []
        for b in bonds:
            b1 = str(elem_dict[b[0]]) + str(ind_dict[b[0]])
            b2 = str(elem_dict[b[1]]) + str(ind_dict[b[1]])
            abf.append([b1, b2] + b[2:])

        all_coords_rel = all_coords - node_cent
        aff_all = np.dot(all_coords_rel, rot) + cvec

        laff_all = np.c_[anf, aff_all, chg, all_inds, [bbind] * len(anf)]

        placed_nbb_coords_extend(laff_all)
        all_bonds_extend(abf)
        ind_seg = ind_seg + len(all_names)

    return placed_nbb_coords, all_bonds

def place_and_calculate_edges(evecs, nnodes):

    placed_ebb_coords = []
    placed_ebb_coords_extend = placed_ebb_coords.extend
    all_bonds = []
    all_bonds_extend = all_bonds.extend
    ind_seg = nnodes
    bbind = -1

    modified_list = []
    removed_list = []

    for e in evecs:

        bbind = bbind - 1
        index, cif, ecoords, evec = e
        ll = 0

        for v in evec:
            mag = np.linalg.norm(v - np.average(evec, axis=0))
            if mag > ll:
                ll = mag

        temp_edge_path = os.path.join('edges', cif)
        temp_edge = Edge(edge_path=temp_edge_path)

        bbxvec = np.array(temp_edge.calculate_cartesian_coords(False))

        edge_cent = np.mean(bbxvec, axis=0)
        bbxvec_rel = bbxvec - edge_cent

        nbbxvec = []
        for v in bbxvec_rel:
            n_v = np.linalg.norm(v)
            if n_v > 1e-8:
                nbbxvec.append(ll * (v / n_v))
            else:
                nbbxvec.append(np.zeros(3))
        nbbxvec = np.array(nbbxvec)

        min_dist, rot, tran = minimize_rms_superimpose(nbbxvec, evec, center=True)
        all_bb, modified_point = temp_edge.crystal_to_array()
        all_coords = np.array([v[1] for v in all_bb])
        all_inds = np.array([v[0] for v in all_bb])

        mdp, rmp = temp_edge.identify_modified_points()# modified points, removed points

        chg, elem = temp_edge.extract_edge_info()
        all_names = [o + re.sub('[A-Za-z]', '', p) for o, p in zip(elem, all_inds)]

        all_names_indices = np.array([int(re.sub('[A-Za-z]', '', e)) for e in all_names]) + ind_seg

        elem_dict = dict((k, '') for k in all_inds)

        for i, j in zip(all_inds, elem):
            elem_dict[i] = j

        ind_dict = dict((k, '') for k in all_inds)

        for i, j in zip(all_inds, all_names_indices):
            ind_dict[i] = j

        for m in mdp:
            if m in ind_dict:
                if m in elem_dict:
                    modified_list.append([m, elem_dict[m], ind_dict[m]])

        for r in rmp:
            if r in ind_dict:
                if r in elem_dict:
                    removed_list.append([r, elem_dict[r], ind_dict[r]])

        bonds = temp_edge.extract_edge_bonds()
        anf = [str(elem_dict[n]) + str(ind_dict[n]) for n in all_inds]

        abf = []
        for b in bonds:
            b1 = str(elem_dict[b[0]]) + str(ind_dict[b[0]])
            b2 = str(elem_dict[b[1]]) + str(ind_dict[b[1]])
            abf.append([b1, b2] + b[2:])

        all_coords_rel = all_coords - edge_cent

        aff_all = np.dot(all_coords_rel, rot) + tran

        laff_all = np.c_[anf, aff_all, chg, all_inds, [bbind] * len(anf)]

        placed_ebb_coords_extend(laff_all)
        all_bonds_extend(abf)
        ind_seg = ind_seg + len(all_names)

    return placed_ebb_coords, all_bonds, modified_list, removed_list

def get_template_topology_info(alpha, num_edges, num_vertices, sc_covar, Bstar_inv, sc_abc, sc_Roman, start, TG,
                               template_name, g, ea_dict):

    sc_a, sc_b, sc_c = sc_abc
    sc_alpha, sc_beta, sc_gamma = sc_Roman

    sc_Alpha = np.r_[
        alpha[0:num_edges - num_vertices + 1, :], sc_covar]
    sc_omega_plus = np.dot(Bstar_inv, sc_Alpha)

    sc_unit_cell = fractional_to_cartesian_matrix(sc_a, sc_b, sc_c, sc_alpha, sc_beta, sc_gamma)

    scaled_coords = calculate_scaled_coordinates(start, TG, sc_omega_plus,
                                                 (sc_a, sc_b, sc_c, sc_alpha, sc_beta, sc_gamma), num_vertices,
                                                 template_name, g)

    nvecs, evecs = compute_scaled_node_and_edge_vectors(scaled_coords, sc_omega_plus, sc_unit_cell, ea_dict)

    placed_nodes, node_bonds = place_and_calculate_nodes(nvecs)  # CHARGES = True,ORIENTATION_DEPENDENT_NODES = False

    placed_edges, edge_bonds, modified, removed = place_and_calculate_edges(evecs, len(placed_nodes))

    return sc_unit_cell, placed_nodes, node_bonds, placed_edges, edge_bonds, modified, removed

def calculate_X_X_bonds_and_update_charges(placed_all, bonds_all, sc_unit_cell, max_length, nconnections,
                                           num_possible_XX_bonds, bond_tol=2.5):

    BOND_TOL = bond_tol

    logging.debug('start forming X-X bonds...')

    fixed_bonds, nbcount, bond_check_passed = build_graph_and_validate_bonds(placed_all, bonds_all,
                                                                             sc_unit_cell, max_length,
                                                                             BOND_TOL, nconnections,
                                                                             num_possible_XX_bonds)

    logging.debug('there were %d X-X bonds formed', nbcount)

    if bond_check_passed:
        logging.warning('bond check passed,then computing charge parameter...')
        bond_check_code = ''
    else:
        logging.warning('bond check failed, attempting distance search bonding...')
        fixed_bonds, nbcount = find_and_calculate_bonds(placed_all, bonds_all, sc_unit_cell, BOND_TOL)
        bond_check_code = '_BOND_CHECK_FAILED'
        logging.warning('there were %d X-X bonds formed', nbcount)

    fc_placed_all, netcharge, onetcharge, rcb = adjust_node_charges_to_neutral(placed_all)

    fixed_bonds = calculate_and_fix_bond_symmetry(fixed_bonds, placed_all, sc_unit_cell)

    remove_net = choice(range(len(fc_placed_all)))
    fc_placed_all[remove_net][4] -= np.round(netcharge, 4)

    return fixed_bonds, fc_placed_all, bond_check_code

def calculate_periodic_distance_and_symmetry(vec1, vec2):

    dX, dY, dZ = vec1 - vec2

    if dX > 0.5:
        s1 = 1
        ndX = dX - 1.0
    elif dX < -0.5:
        s1 = -1
        ndX = dX + 1.0
    else:
        s1 = 0
        ndX = dX

    if dY > 0.5:
        s2 = 1
        ndY = dY - 1.0
    elif dY < -0.5:
        s2 = -1
        ndY = dY + 1.0
    else:
        s2 = 0
        ndY = dY

    if dZ > 0.5:
        s3 = 1
        ndZ = dZ - 1.0
    elif dZ < -0.5:
        s3 = -1
        ndZ = dZ + 1.0
    else:
        s3 = 0
        ndZ = dZ

    if s1 == 0 and s2 == 0 and s3 == 0:
        sym = '.'
    else:

        sym = '1_' + str(s1 + 5) + str(s2 + 5) + str(s3 + 5)

    return np.array([ndX, ndY, ndZ]), sym

def round_to_decimal_places(n, decimals=0):

    multiplier = 10 ** decimals

    return math.floor(n * multiplier + 0.5) / multiplier

def round_down_to_decimal_places(n, decimals=0):

    multiplier = 10 ** decimals

    return math.ceil(n * multiplier - 0.5) / multiplier

def calculate_scaled_coordinates(start, TG, sc_omega_plus, uc_params, num_vertices, template, g):
    sc_a, sc_b, sc_c, sc_alpha, sc_beta, sc_gamma = uc_params

    shortest_path_dict = dict(nx.shortest_path(TG))

    SN = sorted(TG.nodes(), key=lambda x: int(re.sub('[A-Za-z]', '', str(x))))

    sequential_paths = [(SN[i], SN[i + 1]) for i in range(num_vertices) if i + 1 < num_vertices]

    start = np.asarray(start)
    cnd = nx.get_node_attributes(TG, 'cifname')

    coords = [[sequential_paths[0][0], cnd[sequential_paths[0][0]], start,
               [(e[2]['index'], e[2]['pd'], e[2]['cifname']) for e in TG.edges(data=True) if
                sequential_paths[0][0] in (e[0], e[1])]]]

    coords_append = coords.append
    already_placed = [sequential_paths[0][0]]
    already_placed_append = already_placed.append

    for st in sequential_paths:
        path = shortest_path_dict[st[0]][st[1]]
        lp = len(path)
        traverse = [(path[i], path[i + 1]) for i in range(lp) if i + 1 < lp]
        for e0 in traverse:
            s, e = e0
            edict = TG[s][e]

            key = list(edict.keys())[0]
            ind = edict[key]['index']
            positive_direction = edict[key]['pd']

            if (s, e) == positive_direction:
                direction = 1
            elif (e, s) == positive_direction:
                direction = -1
            else:
                raise ValueError('Error in defining an edge traversal in omega_to_coords.py')
            start = start + direction * sc_omega_plus[ind - 1]

            if e0[1] not in already_placed:
                coords_append([e0[1], cnd[e0[1]], start,
                               [(edge[2]['index'], edge[2]['pd'], edge[2]['cifname']) for edge in TG.edges(data=True) if
                                e0[1] in (edge[0], edge[1])]])
                already_placed_append(e0[1])

    norm_coords = []
    norm_coords_append = norm_coords.append

    for line in coords:
        vec = []
        vec_append = vec.append
        v = line[0]

        for dim in line[2]:
            if dim < 0.0:
                while dim < 0:
                    dim = dim + 1.0
            elif dim >= 1.0:
                while dim >= 1.0:
                    dim = dim - 1.0
            vec_append(dim)
        norm_coords_append([line[0], line[1], vec, line[3]])

    norm_coords = sorted(norm_coords, key=lambda x: int(re.sub('[A-Za-z]', '', str(x[0]))))

    return norm_coords