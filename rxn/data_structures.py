import networkx as nx
from networkx import Graph, DiGraph
import networkx.algorithms.isomorphism as iso
from networkx.readwrite import json_graph

from imolecule.notebook import generate


class MolGraph(Graph):
    def __init__(self, data_dict=None):
        Graph.__init__(self)
        if data_dict is not None:
            self.from_dict(data_dict)

    ### Read/Write functions ###
    def to_dict(self):
        """Construct molecule dict with atoms/bonds from molecule graph"""
        mol = {}
        atoms = []
        bonds = []
        nodes = []
        for n in self.nodes():
            if self.node[n]:
                atoms.append(self.node[n].copy())
                nodes.append(n)

        for e in self.edges():
            idxs = [nodes.index(e[0]), nodes.index(e[1])]
            bonds.append({'atoms': idxs, 'order': 1})
        mol['atoms'] = atoms
        mol['bonds'] = bonds
        return mol

    def from_dict(self, moldict):
        """Construct molecule graph from molecule dict with atoms/bonds"""
        bonds = moldict['bonds']
        atoms = moldict['atoms']
        nodes = range(len(moldict['atoms']))
        node_attr = []
        for ai in atoms:
            node = ai['element'] + str(ai['location'])
            attr = ai
            node_attr.append((node, attr))

        self.add_nodes_from(node_attr)

        for b in bonds:
            idxs = b['atoms']
            node_1 = node_attr[idxs[0]][0]
            node_2 = node_attr[idxs[1]][0]
            self.add_edge(node_1, node_2)

    def generate(self, data, data_format):
        """Wrapper around imolecule generate function"""
        moldict = eval(generate(data, data_format))
        self.from_dict(moldict)
        return self

    ### Data structure functions ###
    def __str__(self):
        """Shorthand linear text representation of the molecule"""
        # possible to use SMILES notation instead?
        txt = ''
        for i in nx.dfs_postorder_nodes(self, self.nodes()[0]):
            txt += self.node[i]['element']
        return txt

    def __repr__(self):
        """Comprehensive representation of object. eval(repr(x)) == x."""
        return str('MolGraph({})'.format(self.to_dict()))

    def __eq__(self, other):
        """Check equivalencey between two molecules. Note that this is
        topological equivalency (e.g. atom/bonds) as defined by graph
        isomorphism. Atoms are considered equivalent if the element and
        charge are the same."""
        nm = iso.categorical_node_match('element', 'charge')
        return nx.is_isomorphic(self, other, node_match=nm)


class RxnGraph(DiGraph):

    ### These read/write functions need cleanup ###
    def from_rxn_list(self, rxn_list):
        for rxn in rxn_list:
            reactants, products = rxn
            rxn_name = '+'.join([str(r) for r in reactants])
            rxn_name += '->'
            rxn_name += '+'.join([str(p) for p in products])

            name_0 = rxn_name
            i = 1
            while rxn_name in self.nodes():
                rxn_name = name_0 + '(' + str(i) + ')'
                i += 1

            self.add_node(rxn_name, attr_dict={'type': 'reaction'})

            for r in reactants + products:
                assert hasattr(r, 'is_directed')  # ensure that r is MolGraph
                attrs = {'type': 'molecule', 'graph': r}
                if r in reactants:
                    attrs['molecule_type'] = 'reactant'
                    self.add_node(str(r), attrs)
                    self.add_edge(str(r), rxn_name)
                if r in products:
                    attrs['molecule_type'] = 'product'
                    self.add_node(str(r), attrs)
                    self.add_edge(rxn_name, str(r))

    def to_rxn_list(self):
        rxns = [r for r in self.nodes() if self.node[r]['type'] == 'reaction']
        all_rxns = []
        for rxn in rxns:
            prods = [self.node[p]['graph'] for p in self.successors(rxn)]
            reacts = [self.node[p]['graph'] for p in self.predecessors(rxn)]
            all_rxns.append([reacts, prods])
        return all_rxns

    @staticmethod
    def node_matcher(n1, n2):
        """Define whether or not nodes are equivalent. Static because
        it does not depend on the properties of the entire graph"""
        if n1['type'] != n2['type']:
            # reactions are not the same as molecules
            return False
        elif n1['type'] == n2['type'] == 'molecule':
            # molecules are only equivalent if their graphs are.
            # see MolGraph.__eq__
            return n1['graph'] == n2['graph']
        elif n1['type'] == n2['type'] == 'reaction':
            # treat all reaction nodes as equivalent
            return True

    def __str__(self):
        """Shorthand text representation of the reactions"""
        rxns = []
        for n in self.nodes():
            if self.node[n]['type'] == 'reaction':
                rxns.append(str(n))
        txt = '\n'.join(rxns)
        return txt

    def __repr__(self):
        data = json_graph.node_link_data(self)
        return repr(data)

    def __eq__(self, other):
        """Check equivalencey between two reaction networks.
        Molecules are considered equivalent based on topology."""
        return nx.is_isomorphic(self, other, node_match=self.node_matcher)
