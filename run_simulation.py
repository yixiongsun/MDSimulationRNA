from clean_pdb import clean_pdb
import os.path
from pdb import PDB
from parameters import Parameters
import numpy as np
import taichi as ti
import itertools
import sys

arch = sys.argv[6]
if arch == "gpu":
    ti.init(arch=ti.gpu)

else:
    ti.init(arch=ti.cpu)

# Main file to run simulation

# load parameters
parameters = Parameters("data/top_all36_na.rtf", "data/par_all36_na.prm")

filename = sys.argv[1]

# check if file already cleaned
#filename = "1ivs.cif"
cleaned_filename = filename[:-4] + "_cleaned.cif"
if not os.path.exists(cleaned_filename):
    chain_id = sys.argv[2]

    if chain_id is None:
        raise Exception("Chain ID must be provided!")
    clean_pdb(filename, cleaned_filename, chain_id, parameters)

h = float(sys.argv[3])
steps = int(sys.argv[4])
save_step = int(sys.argv[5])

# Load PDB file
print("initializing")
pdb = PDB(cleaned_filename)
N = pdb.positions.shape[0]

def parse_name(name):
    if name == "OP1":
        name = "O1P"
    if name == "OP2":
        name = "O2P"
    return name

#get neighbours of i and then angles
all_angles = []
for i in range(0, N):
    neighbours = []
    for j in range(0, N):
        a_i = pdb.atoms[i]
        a_j = pdb.atoms[j]

        distance = 5
        # only if its in the same residue or specific atoms: P,

        if a_i.parent.id == a_j.parent.id:
            distance = parameters.get_distance(a_i.parent.resname, parse_name(a_i.name), parse_name(a_j.name))
        # otherwise check if i j are consecutive residues and thus in range
        elif a_j.parent.id[1] == a_i.parent.id[1] + 1:
            # check if i attached to j
            distance = parameters.get_distance(a_i.parent.resname, "+P", parse_name(a_i.name))
            # add distance of next residue
            distance += parameters.get_distance(a_j.parent.resname, "P", parse_name(a_j.name))
        elif a_j.parent.id[1] + 1 == a_i.parent.id[1]:
            distance = parameters.get_distance(a_j.parent.resname, "+P", parse_name(a_j.name))
            distance += parameters.get_distance(a_i.parent.resname, "P", parse_name(a_i.name))
        else:
            if a_j.parent.id[1] > a_i.parent.id[1] + 1:
                break
        if distance == 1:
            neighbours.append(j)

    neighbours_pairs = list(itertools.combinations(neighbours, 2))
    if len(neighbours_pairs) > 1:
        for neighbours in neighbours_pairs:
            a_i = pdb.atoms[i]
            a_j = pdb.atoms[neighbours[0]]
            a_k = pdb.atoms[neighbours[1]]

            # get angle parameters
            kthetha, theta0 = parameters.get_angle(a_i.parent.resname, parse_name(a_i.name), parse_name(a_j.name), parse_name(a_k.name))
            all_angles.append([neighbours[0], i, neighbours[1], kthetha, theta0])

np_angles = np.asarray(all_angles)
has_angles = False
N_angles = np_angles.shape[0]



angles = ti.Vector.field(5, dtype=ti.f32, shape=N_angles)
# get positions from pdb
pos = ti.Vector.field(3, dtype=ti.f32, shape=N, needs_grad=True)
p = ti.Vector.field(6, dtype=ti.f32, shape=(N, N))
force = ti.Vector.field(3, dtype=ti.f32, shape=N)

a = ti.Vector.field(3, dtype=ti.f32, shape=N)
a_next = ti.Vector.field(3, dtype=ti.f32, shape=N)
v = ti.Vector.field(3, dtype=ti.f32, shape=N)
mass = ti.field(ti.f32, shape=N)

U = ti.field(ti.f32, (), needs_grad=True)
intra = ti.Vector.field(2, dtype=ti.f32, shape=(N, N))

pi = np.pi
# angles store (i, k), j is the central atom, Ktheta, Theta0


pos.from_numpy(pdb.positions)
params = np.zeros((N, N, 6))
angles.from_numpy(np_angles)

# force_ij is skew symmetric because f_ij = - f_ji
#force_ij = np.zeros((N, N, 3))


vel = np.zeros((N, 3))
masses = 0
momentum = np.zeros(3)
avg_energy = 0
for i in range(0, N):  #
    mass[i] = parameters.get_mass(pdb.atoms[i].parent.resname, parse_name(pdb.atoms[i].name))

    # initialize velocities, random distribution
    vel[i] = np.random.rand(3) - 0.5
    masses += mass[i]
    momentum += mass[i] * vel[i]
    avg_energy += mass[i] * np.multiply(vel[i], vel[i])

    neighbours = []

    for j in range(i + 1, N):  # no interactions with itself
        # check if in range/interacts
        a_i = pdb.atoms[i]
        a_j = pdb.atoms[j]

        distance = 5
        # only if its in the same residue or specific atoms: P,

        if a_i.parent.id == a_j.parent.id:
            distance = parameters.get_distance(a_i.parent.resname, parse_name(a_i.name), parse_name(a_j.name))

        # otherwise check if i j are consecutive residues and thus in range
        if a_j.parent.id[1] == a_i.parent.id[1] + 1:
            # check if i attached to j
            distance = parameters.get_distance(a_i.parent.resname, "+P", parse_name(a_i.name))
            # add distance of next residue
            distance += parameters.get_distance(a_j.parent.resname, "P", parse_name(a_j.name))

        # check if bond exists between i and j, ie dist = 1
        if distance == 1:
            a2_name = "+P" if (a_j.name == "P" and a_j.parent.id[1] == a_i.parent.id[1] + 1) else parse_name(a_j.name)
            kb, b0 = parameters.get_bond_properties(a_i.parent.resname, parse_name(a_i.name), a2_name)
            intra[i, j][0] = kb
            intra[i, j][1] = b0
            intra[j, i][0] = kb
            intra[i, j][1] = b0

        # take into account UB for angle correction
        if distance == 2:
            a2_name = "+P" if (a_j.name == "P" and a_j.parent.id[1] == a_i.parent.id[1] + 1) else parse_name(a_j.name)
            kb, b0 = parameters.get_UB(a_i.parent.resname, parse_name(a_i.name), a2_name)
            intra[i, j][0] = kb
            intra[i, j][1] = b0
            intra[j, i][0] = kb
            intra[i, j][1] = b0

        if distance < 4:
            continue


        elif distance == 4:
            eps_i, rmin_i = parameters.get_LJ_parameters(a_i.parent.resname, parse_name(a_i.name), True)
            eps_j, rmin_j = parameters.get_LJ_parameters(a_j.parent.resname, parse_name(a_j.name), True)
            q_i = parameters.get_charge(a_i.parent.resname, parse_name(a_i.name))
            q_j = parameters.get_charge(a_j.parent.resname, parse_name(a_j.name))
            #p[i, j] = ti.Vector([rmin_i, rmin_j, eps_i, eps_j, q_i, q_j])
            params[i,j] = np.asarray([rmin_i, rmin_j, eps_i, eps_j, q_i, q_j])

        else:
            eps_i, rmin_i = parameters.get_LJ_parameters(a_i.parent.resname, parse_name(a_i.name))
            eps_j, rmin_j = parameters.get_LJ_parameters(a_j.parent.resname, parse_name(a_j.name))
            q_i = parameters.get_charge(a_i.parent.resname, parse_name(a_i.name))
            q_j = parameters.get_charge(a_j.parent.resname, parse_name(a_j.name))
            #p[i, j] = ti.Vector([rmin_i, rmin_j, eps_i, eps_j, q_i, q_j])
            params[i,j] = np.asarray([rmin_i, rmin_j, eps_i, eps_j, q_i, q_j])





vcm = momentum / masses
avg_energy /= N
kboltz = 0.0083144621
temp = 100
scale = np.full(3, 3 * kboltz * temp)
scale = np.sqrt(np.divide(scale, avg_energy))

for i in range(0, N):
    vel[i] -= vcm
    vel[i] = vel[i] * scale

v.from_numpy(vel)
p.from_numpy(params)

# calculate forces
# to calculate forces, have matrix then sum up after for efficiency, some values don't calculate twice
# return a vector of forces
# takes in pdb

# function to calculate 12-6 LJ potential (van der waals)
# direction of lennard jones force is radial, either towards or away
# modified LJ potential is U = eps_ij (Rmin_ij ^ 12 * r_ij ^ -12 - 2 * Rmin_ij^6 * r_ij ^ -6)
# LJ force is then: F = eps_ij( -12 Rmin_ij ^ 12 * r_ij *^-13 + 12 *  Rmin_ij^6 * r_ij ^ -7)
# to get direction of LJ force: r_ij scaled to unit length (vector pointing from j to i)
@ti.func
def lennard_jones_force(i, j):
    rmin_i = p[i, j][0]
    rmin_j = p[i, j][1]
    eps_i = p[i, j][2]
    eps_j = p[i, j][3]
    r_i = pos[i]
    r_j = pos[j]
    out = ti.Vector([0.0, 0.0, 0.0])
    if eps_i == 0.0:
        out = ti.Vector([0.0, 0.0, 0.0])
    else:
        eps_ij = ti.sqrt(eps_i * eps_j)
        Rmin_ij = rmin_i + rmin_j

        r_ij = r_j - r_i
        dist = r_ij.norm()
        force_magnitude = eps_ij * (-12 * (Rmin_ij ** 12) / (dist ** 13) + 12 * (Rmin_ij ** 6) / (dist ** 7))
        r_norm = r_ij / dist
        out = force_magnitude * r_norm

    return out



# function to calculate electrostatic force
# force_i = ke q_iq_j/r_ij^2, in the direction of r_ij vector pointing from j to i
# ke = 1/4pie0
@ti.func
def coulomb_force(i, j):
    q_i = p[i, j][4]
    q_j = p[i, j][5]
    r_i = pos[i]
    r_j = pos[j]
    # ke in KJ * A / mol e^2
    out = ti.Vector([0.0, 0.0, 0.0])

    if q_i == 0.0:
        out = ti.Vector([0.0, 0.0, 0.0])
    else:
        #ke = 13.893548
        r_ij = r_i - r_j
        dist = r_ij.norm()
        force_magnitude = q_i * q_j / (dist**2)
        r_norm = r_ij / dist
        out = force_magnitude * r_norm
    # convert distances from angstroms to meters
    # 1 angstrom = 1e-10 meters


    return out


@ti.kernel
def intermolecular_forces():

    for i in range(0, N):
        for j in range(i+1, N):
            f_ij = coulomb_force(i, j) + lennard_jones_force(i, j)
            force[i] = force[i] + f_ij
            force[j] = force[j] - f_ij

@ti.kernel
def intramolecular_forces():
    for i in range(0, N):
        for j in range(i + 1, N):
            kb = intra[i, j][0]
            b0 = intra[i, j][1]
            if kb != 0:
                r = pos[j] - pos[i]
                r_norm = r.norm()
                U[None] += kb * ((r_norm - b0) ** 2)
    # angle potential

    
    for i in range(0, N_angles):
        ri = pos[int(angles[i][0])]
        rj = pos[int(angles[i][1])]
        rk = pos[int(angles[i][2])]
        rij = ri - rj
        rkj = rk - rj
        theta = ti.acos(rij.dot(rkj) / (rij.norm() * rkj.norm()))
        kthetha = angles[i][3]
        theta0 = pi * angles[i][4] / 180
        U[None] += kthetha * ((theta - theta0) ** 2)


@ti.kernel
def velocity_verlet_positions(h: ti.f32):
    # integrate
    # p(t + h) = p(t) + h * v(t) + 0.5 * a(t) h ^ 2
    # v(t + h) = v(t) + 0.5 (a(t) + a(t + h)) * h
    for i in range(0, N):
        pos[i] = pos[i] + h * v[i] + 0.5 * h * h * a[i]

@ti.kernel
def velocity_verlet_velocities(h: ti.f32):
    for i in range(0, N):
        v[i] = v[i] + 0.5 * h * (a[i] + a_next[i])


@ti.kernel
def calculate_acceleration():
    for i in range(0, N):
        a_next[i][0] = (force[i][0] - pos.grad[i][0])/ mass[i]
        a_next[i][1] = (force[i][1] - pos.grad[i][0])/ mass[i]
        a_next[i][2] = (force[i][2] - pos.grad[i][0])/ mass[i]
        force[i][0] = 0
        force[i][1] = 0
        force[i][2] = 0

@ti.kernel
def update_acceleration():
    for i in range(0, N):
        a[i][0] = a_next[i][0]
        a[i][1] = a_next[i][1]
        a[i][2] = a_next[i][2]

print("init finished")

# Simulation loop
# 1. Set initial conditions from PDB
# 2. Calculate forces
# 3. Integrate over timestep
# 4. Save into new structure
# end of loop, save to PDB structures and render with pymol


intermolecular_forces()
with ti.Tape(U):
    intramolecular_forces()
for i in range(0, N):
    a[i][0] = (force[i][0] - pos.grad[i][0]) / mass[i]
    a[i][1] = (force[i][1] - pos.grad[i][0]) / mass[i]
    a[i][2] = (force[i][2] - pos.grad[i][0]) / mass[i]
    force[i][0] = 0
    force[i][1] = 0
    force[i][2] = 0

for i in range(0, steps):

    # h = 1 ~ 50 fs
    # h = 0.01 ~ 0.5 fs
    velocity_verlet_positions(h)
    # calculate a(t+h)
    intermolecular_forces()
    with ti.Tape(U):
        intramolecular_forces()
    calculate_acceleration()
    velocity_verlet_velocities(h)

    update_acceleration()
    if i % save_step == 0:
        print("Time step:", i)
        new_pos = pos.to_numpy()
        if np.isnan(np.min(new_pos)):
            break
        pdb.set_new_positions(new_pos)


output = filename.split("/")[-1]
pdb.save("output/" + output[:-4] + "_simulation.cif")

