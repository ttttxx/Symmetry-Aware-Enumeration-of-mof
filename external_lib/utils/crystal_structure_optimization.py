import itertools
import logging
import random

import numpy as np
from numpy.linalg import norm
from scipy.optimize import minimize

from external_lib.utils.Bio import SVDSuperimposer
from external_lib.utils.parameter import ShapeMismatchError

pi = np.pi

def force_shape_match(a1, a2):
    if len(a1) == len(a2):
        return a1

    raise ShapeMismatchError(
        f"配体连接点数 ({len(a1)}) 与拓扑锚点数 ({len(a2)}) 不匹配，无法安全对齐。"
        f"该构象将被丢弃。"
    )

def kabsch_align_core(P, Q, center=True):
    P = np.asarray(P)
    Q = np.asarray(Q)

    if center:
        P_centroid = np.mean(P, axis=0)
        Q_centroid = np.mean(Q, axis=0)
    else:
        P_centroid = np.zeros(3)
        Q_centroid = np.zeros(3)

    P_centered = P - P_centroid
    Q_centered = Q - Q_centroid

    if len(P) == 2:
        v1 = P_centered[1] - P_centered[0]
        v2 = Q_centered[1] - Q_centered[0]
        v1 = v1 / (np.linalg.norm(v1) + 1e-12)
        v2 = v2 / (np.linalg.norm(v2) + 1e-12)

        c = np.dot(v1, v2)
        if c < -0.999999:
            orth = np.array([1.0, 0.0, 0.0]) if abs(v1[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
            axis = np.cross(v1, orth)
            axis = axis / np.linalg.norm(axis)
            vx = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
            R = np.eye(3) + 2 * np.dot(vx, vx)
        else:
            v = np.cross(v1, v2)
            s = np.linalg.norm(v)
            if s < 1e-8:
                R = np.eye(3)
            else:
                vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
                R = np.eye(3) + vx + np.dot(vx, vx) * ((1 - c) / (s ** 2))

        t = Q_centroid - np.dot(P_centroid, R.T)
        rms = np.sqrt(np.mean(np.sum((np.dot(P_centered, R.T) - Q_centered) ** 2, axis=1)))
        return rms, R.T, t

    else:

        H = np.dot(P_centered.T, Q_centered)
        U, S, Vt = np.linalg.svd(H)
        R = np.dot(Vt.T, U.T)

        if np.linalg.det(R) < 0:
            Vt[2, :] *= -1
            R = np.dot(Vt.T, U.T)

        t = Q_centroid - np.dot(P_centroid, R.T)
        rms = np.sqrt(np.mean(np.sum((np.dot(P_centered, R.T) - Q_centered) ** 2, axis=1)))
        return rms, R.T, t

def select_closest_vectors(a1, a2, num):
    dist1 = [(np.linalg.norm(a1[0] - a1[i]), i) for i in range(len(a1))]
    dist2 = [(np.linalg.norm(a2[0] - a2[i]), i) for i in range(len(a2))]

    dist1.sort(key=lambda x: x[0])
    dist2.sort(key=lambda x: x[0])

    vecs1 = np.array([a1[i] for i in [dist1[j][1] for j in range(num)]])
    vecs2 = np.array([a2[i] for i in [dist2[j][1] for j in range(num)]])

    return vecs1, vecs2

def minimize_rms_superimpose(a1, a2, center=True):
    a1 = np.asarray(a1)
    a2 = np.asarray(a2)
    a1 = force_shape_match(a1, a2)

    if len(a1) > 7:
        a1, a2 = select_closest_vectors(a1, a2, 6)

    min_dist = (float('inf'), np.eye(3), np.zeros(3))
    for l in itertools.permutations(a1):
        p = np.asarray(l)
        rms, rot, tran = kabsch_align_core(p, a2, center=center)

        if rms < min_dist[0]:
            min_dist = (rms, rot, tran)

    return min_dist

def calculate_minimum_distance_vectors(a1, a2):
    a1 = np.asarray(a1)
    a2 = np.asarray(a2)

    a1 = force_shape_match(a1, a2)
    mags = [norm(v) for v in a2]

    if len(a1) < 3:
        min_dist = (float('inf'), np.eye(3), np.zeros(3))
        for l in itertools.permutations(a1):

            p = np.array([m * v / (norm(v) + 1e-8) for m, v in zip(mags, l)])
            rms, rot, tran = kabsch_align_core(p, a2)
            if rms < min_dist[0]:
                min_dist = (rms, rot, tran)
        return min_dist

    sup = SVDSuperimposer()
    if len(a1) <= 7:
        min_dist = (float('inf'), np.eye(3), np.zeros(3))
        for l in itertools.permutations(a1):
            p = np.array([m * v / (norm(v) + 1e-8) for m, v in zip(mags, l)])
            sup.set(a2, p)
            sup.run()
            rot, tran = sup.get_rotran()
            rms = sup.get_rms()

            if rms < min_dist[0]:
                min_dist = (rms, rot, tran)
    else:
        a1, a2 = select_closest_vectors(a1, a2, 6)
        mags = [norm(v) for v in a2]
        min_dist = (float('inf'), np.eye(3), np.zeros(3))

        for l in itertools.permutations(a1):
            p = np.array([m * v / (norm(v) + 1e-8) for m, v in zip(mags, l)])
            sup.set(a2, p)
            sup.run()
            rot, tran = sup.get_rotran()
            rms = sup.get_rms()

            if rms < min_dist[0]:
                min_dist = (rms, rot, tran)

    return min_dist

def adjust_coordinates_periodic_boundary(c1, c2):
    if not (isinstance(c1, (list, np.ndarray)) and isinstance(c2, (list, np.ndarray))):
        print("Error: c1 and c2 should be lists or numpy arrays.")
        return None
    if len(c1) < 3 or len(c2) < 3:
        print("Error: c1 and c2 should have at least 3 dimensions.")
        return None

    c2 = np.copy(c2)

    diffa = c1[0] - c2[0]
    diffb = c1[1] - c2[1]
    diffc = c1[2] - c2[2]

    if diffa > 0.5:
        c2[0] = c2[0] + 1.0
    elif diffa < -0.5:
        c2[0] = c2[0] - 1.0

    if diffb > 0.5:
        c2[1] = c2[1] + 1.0
    elif diffb < -0.5:
        c2[1] = c2[1] - 1.0

    if diffc > 0.5:
        c2[2] = c2[2] + 1.0
    elif diffc < -0.5:
        c2[2] = c2[2] - 1.0

    return c2

def calculate_objective_value(V, ncra, ncca, Alpha, ne, nv, Bstar_inv, SBU_IP):
    def metric_tensor(l):
        a, b, c, alpha, beta, gamma = l
        Z = np.array([[a * a, a * b * np.cos(gamma * (pi / 180)), a * c * np.cos(beta * (pi / 180))],
                      [a * b * np.cos(gamma * (pi / 180)), b * b, b * c * np.cos(alpha * (pi / 180))],
                      [a * c * np.cos(beta * (pi / 180)), b * c * np.cos(alpha * (pi / 180)), c * c]])
        return Z

    P = V[0:6]
    Z = metric_tensor(P)
    X = np.array(V[6:])

    X = X.reshape(ncra, ncca)
    var_alpha = np.r_[Alpha[0:ne - nv + 1, :], X]

    omp = np.dot(Bstar_inv, var_alpha)
    g = np.dot(omp, np.dot(Z, omp.T))

    O1 = 0
    for sbu in SBU_IP:
        O2 = 0
        for i in sbu[1]:
            O2 = O2 + (i[2] - g[i[0] - 1][i[1] - 1]) ** 2
        O1 = O1 + O2

    return O1

def optimize_and_scale_unit_cell(all_SBU_coords, a, b, c, ang_alpha, ang_beta, ang_gamma, max_le, num_vertices, Bstar,
                                 alpha, num_edges):
    max_length = 0
    for line in all_SBU_coords:
        for length in [np.linalg.norm(s[1]) for s in line[1]]:
            if length > max_length:
                max_length = length

    scale_guess = (max_length / max_le) * 1.0

    all_SBU_ip = []
    all_SBU_ip_append = all_SBU_ip.append
    for sbu in all_SBU_coords:
        SBU_ip = []
        SBU_ip_append = SBU_ip.append
        w = len(sbu[1])
        for i in range(w):
            ivec = sbu[1][i][1]
            iind = sbu[1][i][0]
            for j in range(i, w):
                jvec = sbu[1][j][1]
                jind = sbu[1][j][0]
                dot = np.dot(ivec, jvec)
                SBU_ip_append([iind, jind, dot])
        all_SBU_ip_append((sbu[0], SBU_ip))
    ncra = num_vertices - 1
    ncca = 3

    covars_values = []
    covars_values_append = covars_values.append

    for i in range(ncra):
        for j in range(ncca):
            covars_values_append(0)

    ucvars = [a, b, c, ang_alpha, ang_beta, ang_gamma]

    max_a = 2 * scale_guess * a
    max_b = 2 * scale_guess * b
    max_c = 2 * scale_guess * c

    MIN_CELL_LENGTH = 5.0
    uc_bounds = ((MIN_CELL_LENGTH, max_a), (MIN_CELL_LENGTH, max_b), (MIN_CELL_LENGTH, max_c), (20, 160), (20, 160),
                 (20, 160))
    init_variables = [scale_guess * a, scale_guess * b, scale_guess * c, ang_alpha, ang_beta, ang_gamma] + covars_values
    x_bounds = tuple([(-10.0, 10.0) for x in covars_values])
    bounds = uc_bounds + x_bounds
    Bstar_inv = np.linalg.inv(Bstar)

    logging.debug('start scaling unit cell and vertex positions...')
    logging.debug('optimizing with local minimization algorithm [L-BFGS-B]')

    niter = 2
    uc_press = 0.05
    covars_perturb = 0.0001
    callbackresults = [
        [init_variables,
         calculate_objective_value(init_variables, ncra, ncca, alpha, num_edges, num_vertices, Bstar_inv, all_SBU_ip)]]

    def callbackF(X):
        var_string = ''
        for val in X:
            var_string += str(val) + '   '
        callbackresults.append(
            [X, calculate_objective_value(X, ncra, ncca, alpha, num_edges, num_vertices, Bstar_inv, all_SBU_ip)])

    for it in range(niter):
        res = minimize(calculate_objective_value, init_variables,
                       args=(ncra, ncca, alpha, num_edges, num_vertices, Bstar_inv, all_SBU_ip), method='L-BFGS-B',
                       bounds=bounds, options={'disp': False}, callback=callbackF)

        uc_params = res.x[0:6]
        uc_params = [i - (i * uc_press) for i in uc_params]
        mult = [random.choice([-1, 1]) for i in range(len(res.x[6:]))]
        covars = [i + j * (i * covars_perturb) for i, j in zip(res.x[6:], mult)]
        init_variables = uc_params + covars

    if niter != 0:
        sc_a, sc_b, sc_c, sc_alpha, sc_beta, sc_gamma = res.x[0:6]
        sc_covar = res.x[6:].reshape(ncra, ncca)
    else:
        init_variables = np.asarray(init_variables)
        sc_a, sc_b, sc_c, sc_alpha, sc_beta, sc_gamma = init_variables[0:6]
        sc_covar = init_variables[6:].reshape(ncra, ncca)

    ab = [a / b, sc_a / sc_b]
    ac = [a / c, sc_a / sc_c]
    bc = [b / c, sc_b / sc_c]
    alpha = [ang_alpha, sc_alpha]
    beta = [ang_beta, sc_beta]
    gamma = [ang_gamma, sc_gamma]
    covar = [np.average(abs(np.array(callbackresults[0][0][6:]) - np.array(callbackresults[-1][0][6:])))]
    final_obj = [callbackresults[-1][1]]

    scaling_data = [ab, ac, bc, alpha, beta, gamma, covar, final_obj]

    return (sc_a, sc_b, sc_c, sc_alpha, sc_beta, sc_gamma, sc_covar, Bstar_inv, max_length, callbackresults, ncra, ncca,
            scaling_data)