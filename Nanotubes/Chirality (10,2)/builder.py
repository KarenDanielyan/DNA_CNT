import argparse
import numpy as np
import sys
import subprocess as sp
from math import cos, sin, sqrt, radians, pi
from functools import partial
from os.path import basename

from moleculex import Molecule, Atom

_tubegen_input = """
set format pdb
set chirality {chirality[0]},{chirality[1]}
set units angstrom
set gutter 1.6735,1.6735,0
set shape cubic
set relax_tube yes
set cell_count 1,1,{cell_count_z}
generate
save {prefix}.pdb
exit
"""
_tube_user_input = """!
! User settings
!
set prefix   = %s
set material = %s
set pbc      = %s
set rtube    = %.3f
set ltube    = %.3f
set imtube   = %s
"""

def cGCD(x, y):
    if x == 0 or y == 0:
        if x > y:
            gcd = x
        else:
            gcd = y
    if x > y and x * y != 0:
        small = y
    else:
        small = x
    for i in range(1, small+1):
        if((x % i == 0) and (y % i == 0)):
            gcd = i
    return gcd

def det_thickness(n,m):
    d = cGCD(n,m)
    if (n-m) % ( 3.0 * d) == 0 :
        dr = 3.0 * d
    else:
        dr =  d
    a1 = np.array([2.1315 , 1.23062])
    a2 = np.array([2.1315 ,-1.23062])
    Tv = ((2*m + n) * a1 - (2*n + m ) * a2) / dr
    T = np.linalg.norm(Tv)
    return T

def det_radius(n,m):
    return sqrt(3* (n**2 + n*m + m**2)) * 1.421 / (2 * pi)

def range_type(astr, min=0, max=101):
    value = float(astr)
    if min<= value <= max:
        return value
    else:
        raise argparse.ArgumentTypeError('value not in range %s-%s'%(min,max))

def csv_list_type(string):
    try:
        return [int(_) for _ in string.split(',')]
    except:
        raise ValueError('comma separated integers expected')

def tubegen(args):
    pid = sp.Popen(args.tubegen, stdin=sp.PIPE, stdout=sp.PIPE, close_fds=True, universal_newlines=True)
    (stdin, stdout) = (pid.stdin, pid.stdout)
    stdin.write(_tubegen_input.format(**args.__dict__))
    cap = (args.pbc == 'no')
    polar_carbon = args.pcarbon
    n, m = args.chirality
    cell_copy = args.cell_count_z
    prefix = args.prefix
    pbc = args.pbc

    stdin.flush()
    ztrans = det_thickness(n, m) * int(cell_copy)
    ltube = ztrans
    rtube = det_radius(n , m)
    print 'length_of_tube : %5.3f' % ( det_thickness(n, m) * int(cell_copy))
    print 'radius_of_tube : %5.3f' % rtube
    if det_thickness(n, m) * int(cell_copy) > 400:
        print("this n,m combination will generate too long CNT")
        exit()
    print stdout.read()
    pdbfile = '{}.pdb'.format(args.prefix)
    rtffile = '{}.rtf'.format(args.prefix)
    crdfile = '{}.crd'.format(args.prefix)
    usrfile = 'step1.1_user_input.str'
    mol = Molecule.from_pdb(open(pdbfile))

    atoms = list(mol.atoms())
    if cap:
        terminal_carbons = []
        for atom in atoms:
            neighbors = [_ for _ in mol.graph.neighbors(atom) if _.element.symbol == 'C']
            if atom.element.symbol == 'C':
                num_neighbors = len(neighbors)
                if num_neighbors == 1:
                    neighbor = neighbors[0]
                    mol.remove_bond(atom, neighbor)
                    mol.remove_atom(atom)
                    terminal_carbons.append(neighbor)
                elif num_neighbors == 2:
                    if atom in terminal_carbons:
                        continue
                    terminal_carbons.append(atom)
        for atom in terminal_carbons:
            neighbors = [_ for _ in mol.graph.neighbors(atom) if _.element.symbol == 'C']
            r0 = np.array([atom.x, atom.y, atom.z])
            r1 = np.array([neighbors[0].x, neighbors[0].y, neighbors[0].z])
            r2 = np.array([neighbors[1].x, neighbors[1].y, neighbors[1].z])
            rc = (r1 + r2)/2
            vc = (r0 - rc) / np.linalg.norm(r0 - rc)
            rh = vc * 1.2 + r0

            cap = Atom(name='HT', resname=atom.resname, resnr=atom.resnr, x=rh[0], y=rh[1], z=rh[2])
            mol.add_atom(cap)
            mol.add_bond(atom, cap)
        imtube = 'no'
    else:
        bt_carbons = []
        up_carbons = []
        for atom in atoms:
            neighbors = [_ for _ in mol.graph.neighbors(atom) if _.element.symbol == 'C']
            if atom.element.symbol == 'C':
                num_neighbors = len(neighbors)
                if num_neighbors <= 2:
                    if atom.z < 0:
                        bt_carbons.append(atom)
                    else:
                        up_carbons.append(atom)

        up_pos = np.zeros([len(up_carbons), 3])
        for i, upc in enumerate(up_carbons):
            up_pos[i,0] = upc.x
            up_pos[i,1] = upc.y
            up_pos[i,2] = upc.z
        imp_dic = {}
        for btc in bt_carbons:
            impatch = []
            btpos = np.array([btc.x, btc.y, btc.z + ztrans])
            distmat = np.linalg.norm(up_pos - btpos, axis = 1)
            neighbors = [_ for _ in mol.graph.neighbors(btc) if _.element.symbol == 'C']
            num_neighbors = len(neighbors)
            if num_neighbors == 2:
               dist = np.partition(distmat, 1)[0:1]
               idx = int(np.where(distmat == dist)[0])
               impatch.append(up_carbons[idx])
            elif num_neighbors == 1:
               dist1, dist2 = np.partition(distmat, 1)[0:2]
               idxes = np.where(distmat == dist1)[0]
               if len(idxes) == 2:
                   idx1, idx2 = idxes
               else:
                   idx1 = int(np.where(distmat == dist1)[0])
                   idx2 = int(np.where(distmat == dist2)[0])
               impatch.append(up_carbons[idx1])
               impatch.append(up_carbons[idx2])
            imp_dic[btc] = impatch
        imtube = 'yes'
    if polar_carbon:
        polar_h = []
        for atom in atoms:
            neighbors = [_ for _ in mol.graph.neighbors(atom) if _.element.symbol == 'C']
            if atom.element.symbol == 'C':
                cpos = np.array([atom.x, atom.y])
                cpos_norm = np.linalg.norm(cpos)
                hvec = cpos / cpos_norm
                hvec_p = hvec * 0.65 + cpos
                hvec_p = [hvec_p[0], hvec_p[1], atom.z]
                hp1 = Atom(name='HP', resname=atom.resname, resnr=atom.resnr, x=hvec_p[0], y=hvec_p[1], z=hvec_p[2])
                mol.add_atom(hp1)
                mol.add_bond(atom, hp1)
                hvec_n = hvec * -0.65 + cpos
                hvec_n = [hvec_n[0], hvec_n[1], atom.z]
                hp2 = Atom(name='HP', resname=atom.resname, resnr=atom.resnr, x=hvec_n[0], y=hvec_n[1], z=hvec_n[2])
                mol.add_atom(hp2)
                mol.add_bond(atom, hp2)

    # set atom names
    mol.name_atoms()

    # set resname and segid
    for atom in mol.atoms():
        atom.resname = "TUBE"
        atom.segid = "TUBE"

    mol.to_rtf(open(rtffile, 'w'), polar_carbon)
    if not cap: mol.to_impatch(open(rtffile, 'a'), imp_dic)
    mol.to_pdb(open(pdbfile, 'w'))
    mol.to_crd(open(crdfile, 'w'))
    f = open(usrfile, 'w')

    f.write(_tube_user_input %(basename(args.prefix), "tube", pbc, float(rtube), float(ltube), imtube))
def graphene(args):
    import graphene
    import networkx as nx

    nrings = args.nrings
    defect = args.defect
    dense_defect = args.dense_defect
    cap_hydrogen = args.cap

    g = nx.Graph()
    g.defect_level = defect
    g.dense_defect = dense_defect
    graphene.add_unit(g)
    closed_node = []
    for i in range(nrings):
        node = graphene.find_neighbor(g)
        nodetype = g.node[node]['vertices']
        for j in range(nodetype):
            graphene.add_unit_neighbor(g, node)

        for n in g.nodes():
            if n not in closed_node and g.node[n]['vertices'] == g.degree(n):
                graphene.check_closure(g, n)
                closed_node.append(n)

    h = graphene.atom_graph(g)

    # cap hydrogen
    if cap_hydrogen:
        atoms = list(h.nodes())
        count = len(h.nodes())
        for n in atoms:
            if h.degree(n) != 3:
                nn = 3 - h.degree(n)
                for _ in range(nn):
                    h.add_node(count, type='HT')
                    h.add_edge(n, count)
                    count += 1

    pos=nx.nx_pydot.graphviz_layout(h, prog='neato')
    pdbname = '{}.pdb'.format(args.prefix)
    rtfname = '{}.rtf'.format(args.prefix)
    atom_names = graphene.build_initial_pdb(g, pos, pdbname)
    graphene.build_topology(h, atom_names, rtfname)

def run(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers()
    subparsers.required = True
    subparsers.dest = 'command'

    ptubegen = subparsers.add_parser('tubegen')
    ptubegen.add_argument('--tubegen', default="./tubegen-3.4/src/tubegen", help="path to tubegen executable")
    ptubegen.add_argument('--chirality', default=[3, 3], type=csv_list_type, metavar='[N,M]', help="chirality parameter n, m")
    ptubegen.add_argument('--cell_count_z', default=1, type=int, metavar='N', help="cell count along Z")
    ptubegen.add_argument('--prefix', default="tube", help="prefix for output file")
    ptubegen.add_argument('--pbc', default="no", help="PBC direction (x, y, z, or 'no': default: no)")
    ptubegen.add_argument('--pcarbon', action="store_true", help="polar carbon  direction (default: no)")
    ptubegen.set_defaults(func=tubegen)

    pgraphene = subparsers.add_parser('graphene')
    pgraphene.add_argument('--charmm', default="charmm", help="path to CHARMM executable")
    pgraphene.add_argument('--nrings', default=100, type=int, help="number of rings")
    pgraphene.add_argument('--cap', default=False, action='store_true', help="cap hydrogen if turned on")
    pgraphene.add_argument('--prefix', default="graphene", help="prefix for output file")
    pgraphene.add_argument('--defect', default=0, type=partial(range_type, min=0, max=1), metavar='[0-1]', help="fraction of rings to have defect")
    pgraphene.add_argument('--dense_defect', default=False, action='store_true', help="favors defects to be near each other when turned on")
    pgraphene.set_defaults(func=graphene)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    run()
