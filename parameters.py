import networkx as nx

class Residue(object):

    atoms = {}
    bonds = []
    improper = []

    def __init__(self, name):
        self.name = name
        self.G = nx.Graph()

    def add_atom(self, name, node_type, node_charge):
        #self.atoms[name] = {"type": type, "charge": charge}

        # add nodes to graph
        self.G.add_node(name,type=node_type, charge=node_charge)

    def add_bond(self, a1, a2):


        if a2 == "+P":
            self.G.add_node(a2, type=self.G.nodes["P"]["type"], charge=self.G.nodes["P"]["charge"])
            self.G.add_edge(a1, a2)

           # self.bonds.append({"atoms": (a1, a2), "atom_types": (self.atoms[a1]["type"], self.atoms["P"]["type"])})
        else:
            self.G.add_edge(a1, a2)
            #self.bonds.append({"atoms": (a1, a2), "atom_types": (self.atoms[a1]["type"], self.atoms[a2]["type"])})

    def add_improper(self, a1, a2, a3, a4):
        self.improper.append((a1, a2, a3, a4))

    # Stored as bond type
    def add_bond_property(self, a1, a2, Kb, b0):
        for e in self.G.edges:
            if (self.G.nodes[e[0]]["type"] == a1 and self.G.nodes[e[1]]["type"] == a2) or (self.G.nodes[e[1]]["type"] == a1 and self.G.nodes[e[0]]["type"] == a2):
                self.G[e[0]][e[1]]["Kb"] = Kb
                self.G[e[0]][e[1]]["b0"] = b0

    def add_angle(self, a1, a2, a3, Ktheta, Theta0):
        # a2 is the central atom
        for node in self.G.nodes:
            if self.G.nodes[node]["type"] == a2:
                self.G.nodes[node].setdefault("angle", []).append({"branches":(a1, a3), "Ktheta": Ktheta, "Theta0": Theta0})

    def add_urey_bradley(self, a1, a3, Kub, S0):
        for node in self.G.nodes:
            if self.G.nodes[node]["type"] == a1:
                self.G.nodes[node].setdefault("ub",[]).append({"opposite": a3, "Kub":Kub, "S0": S0})
            elif self.G.nodes[node]["type"] == a3:
                self.G.nodes[node].setdefault("ub",[]).append({"opposite": a1, "Kub":Kub, "S0": S0})

    def add_LJ_parameters(self, a, eps, rmin, eps_14, rmin_14):
        for node in self.G.nodes:
            if self.G.nodes[node]["type"] == a:
                self.G.nodes[node]["LJ"] = {"eps": eps, "rmin": rmin, "eps_14": eps_14, "rmin_14": rmin_14}

    def add_mass(self, a_type, mass):
        for node in self.G.nodes:
            if self.G.nodes[node]["type"] == a_type:
                self.G.nodes[node]["mass"] = mass

    def get_atom(self, atom):
        if atom in self.G.nodes:
            return self.G.nodes[atom]
        return None

    def build_distances(self):
        self.distances = dict(nx.all_pairs_shortest_path_length(self.G))

    def get_distance(self, a1, a2):
        return self.distances[a1][a2]

    def get_LJ_parameters(self, atom, dist_14=False):
        if dist_14:
            return self.G.nodes[atom]["LJ"]["eps_14"], self.G.nodes[atom]["LJ"]["rmin_14"]
        return self.G.nodes[atom]["LJ"]["eps"], self.G.nodes[atom]["LJ"]["rmin"]

    def get_charge(self, atom):
        return self.G.nodes[atom]["charge"]

    def get_mass(self, atom):
        return self.G.nodes[atom]["mass"]

    def get_bond_properties(self, a1, a2):
        bond = self.G[a1][a2]
        return bond["Kb"], bond["b0"]

    def get_angle(self, central, a2, a3):
        t1 = self.G.nodes[a2]["type"]
        t3 = self.G.nodes[a3]["type"]
        angles = self.G.nodes[central]["angle"]
        for a in angles:
            if (t1, t3) == a['branches'] or (t3, t1) == a['branches']:
                return a['Ktheta'], a['Theta0']

        return 0, 0

    def get_UB(self, a1, a2):
        if "ub" in self.G.nodes[a1]:
            UB = self.G.nodes[a1]["ub"]
            for ub in UB:
                if self.G.nodes[a2]["type"] == ub['opposite']:
                    return ub["Kub"], ub["S0"]

        return 0, 0

class Parameters(object):

    residues = {"A": Residue("A"), "C": Residue("C"), "U": Residue("U"), "G": Residue("G")}
    r_names = {"ADE": "A", "CYT": "C", "URA": "U", "GUA": "G"}

    def __init__(self, topology_file, parameter_file):
        self.load_topology(topology_file)
        self.load_parameter(parameter_file)
        for k, r in self.residues.items():
            r.build_distances()


    def load_topology(self, topology_file):
        with open(topology_file, encoding="utf8") as f:
            current_residue = None
            while True:
                line = f.readline()
                if line == "end\n":
                    break

                if "RESI" in line or "PRES" in line:
                    r = line.split()[1]

                    if r in self.r_names:
                        current_residue = self.residues[self.r_names[r]]
                    else:
                        current_residue = None
                else:
                    if current_residue is None:
                        continue
                    if "ATOM" in line:
                        stripped = line.split("!", 1)[0]
                        parse = stripped.split()
                        if "H" in stripped:
                            continue
                        current_residue.add_atom(parse[1], parse[2], float(parse[3]))
                    elif "BOND" in line or "DOUBLE" in line:
                        stripped = line.split("!", 1)[0]
                        bonds = stripped.split()

                        for i in range(1, len(bonds), 2):
                            if "H" in bonds[i] or "H" in bonds[i + 1]:
                                continue
                            current_residue.add_bond(bonds[i], bonds[i + 1])
                    elif "IMPR" in line:
                        stripped = line.split("!", 1)[0]

                        impr = stripped.split()

                        for i in range(1, len(impr), 4):
                            if "H" in impr[i] + impr[i + 1] + impr[i + 2] + impr[i + 3]:
                                continue
                            current_residue.add_improper(impr[i], impr[i + 1], impr[i + 2], impr[i + 3])

    def load_parameter(self, parameter_file):
        with open(parameter_file, encoding="utf8") as f:
            current = None
            while True:
                line = f.readline()
                if line == "END\n":
                    break

                if line[0] == "!" or line == "\n":
                    continue

                if "MASS" in line:
                    stripped = line.split("!", 1)[0]
                    mass = stripped.split()
                    for key, res in self.residues.items():
                        res.add_mass(mass[2], float(mass[3]))
                    continue

                if "BONDS" in line:
                    current = "BOND"
                elif "ANGLES" in line:
                    current = "ANGLES"
                elif "DIHEDRALS" in line:
                    current = "DIHEDRALS"
                elif "NONBONDED" in line:
                    current = "LJ"
                elif "NBFIX" in line:
                    current = None
                else:
                    stripped = line.split("!", 1)[0]

                    if current == "BOND":
                        bond = stripped.split()
                        for key, res in self.residues.items():
                            res.add_bond_property(bond[0], bond[1], float(bond[2]), float(bond[3]))

                    elif current == "ANGLES":
                        angle = stripped.split()
                        # angle is a1-a2-a3, a2 is the center atom
                        # urey bradley term only for a1-a3 properties

                        for key, res in self.residues.items():
                            res.add_angle(angle[0], angle[1], angle[2], float(angle[3]), float(angle[4]))
                            if len(angle) > 5:
                                res.add_urey_bradley(angle[0], angle[1], float(angle[5]), float(angle[6]))

                    elif current == "LJ":
                        angle = stripped.split()
                        for key, res in self.residues.items():
                            if len(angle) > 4:
                                res.add_LJ_parameters(angle[0], float(angle[2]), float(angle[3]), float(angle[5]), float(angle[6]))
                            else:
                                res.add_LJ_parameters(angle[0], float(angle[2]), float(angle[3]), float(angle[2]), float(angle[3]))


    def has_atom(self, res, atom):
        if self.residues[res].get_atom(atom) is not None:
            return True
        return False

    def has_residue(self, res):
        if res in self.residues:
            return True
        return False

    def get_distance(self, res, a1, a2):
        return self.residues[res].get_distance(a1, a2)

    def get_LJ_parameters(self, res, atom, dist_14=False):
        return self.residues[res].get_LJ_parameters(atom, dist_14)

    def get_charge(self, res, atom):
        return self.residues[res].get_charge(atom)

    def get_mass(self, res, atom):
        return self.residues[res].get_mass(atom)

    def get_bond_properties(self, res, a1, a2):
        return self.residues[res].get_bond_properties(a1, a2)

    def get_angle(self, res, central, a2, a3):
        return self.residues[res].get_angle(central, a2, a3)

    def get_UB(self, res, a1, a2):
        return self.residues[res].get_UB(a1, a2)

