from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dag
from matplotlib import pyplot as plt
from qiskit.dagcircuit import DAGOpNode
from collections import defaultdict 
from qiskit.circuit.quantumregister import Qubit
from qiskit import QuantumRegister
from qiskit.transpiler import CouplingMap, Layout
from copy import copy 
from qiskit.circuit.library.standard_gates import SwapGate, CXGate
import networkx as nx
import numpy as np
import pydot

from tree import Tree, Node



path = 'E:/files/quant/code/quantum/qasm_circuits/sqn_258.qasm'

width = 8
depth = 17
score_layer = 8 # 1020 rd84_253.qasm # decay = 0.7 7

# rd84_253: 923 score_layer==7
# sym6_145: 169
# cycle10_2_110: 358

# rd84_253: 969
# sym6_145: 235
# cycle10_2_110: 469
# radd_250: 264
# rd84_142: 42
# adr4_197: 212
# sqn_258: 



# width = 7
# depth = 17
# score_layer = 6 # 1020 rd84_253.qasm # decay = 0.7

# rd84_142: 40 +2
# adr4_197: 228 +15
# sqn_258: 686 +57
# square_root_7: 417 +80
# rd73_252: 432 +93
# radd_250: 271 -31
# rd84_253: 991 +49
# z4_268: 222 +36
# sym6_145: 422 -96
# cycle10_2_110: 463 -93
# misex1_241: 252 +1

# width = 6
# depth = 40 ~ 17
# score_layer = 3 # 1020 rd84_253.qasm # decay = 0.7

# rd84_142: 41 
# adr4_197: 294
# sqn_258: 764 
# square_root_7: 429
# rd73_252: 379
# radd_250: 257
# rd84_253: 1020
# z4_268: 281  
# sym6_145: 311
# cycle10_2_110: 448
# misex1_241: 266

# width = 2 
# depth = 1 
# score_layer = 10 # 1776 rd84_253.qasm # decay = 0.6





if __name__ == '__main__':



    edges = [(0,1),(1,0),(1,2),(2,1),(2,3),(3,2),(3,4),(4,3),(0,5),(5,0),(1,6),(6,1),(2,7),(7,2),(3,8),(8,3),(4,9),(9,4),(1,7),(7,1),(2,6),(6,2),(3,9),(9,3),(4,8),(8,4),(5,6),(6,5),(6,7),(7,6),(7,8),(8,7),(8,9),(9,8),(5,10),(10,5),(6,11),(11,6),(7,12),(12,7),(8,13),(13,8),(9,14),(14,9),(5,11),(11,5),(6,10),(10,6),(7,13),(13,7),(8,12),(12,8),(10,11),(11,10),(11,12),(12,11),(12,13),(13,12),(13,14),(14,13),(10,15),(15,10),(11,16),(16,11),(12,17),(17,12),(13,18),(18,13),(14,19),(19,14),(11,17),(17,11),(12,16),(16,12),(13,19),(19,13),(14,18),(18,14),(15,16),(16,15),(16,17),(17,16),(17,18),(18,17),(18,19),(19,18)]

    coupling_map = CouplingMap(couplinglist=edges)
    print('width:', width)
    print('depth:', depth)
    print('score_layer:', score_layer)
    
    # if coupling_map.graph.has_edge(0, 1):
    #     print(coupling_map.neighbors(1))

    # print(coupling_map.neighbors(7))
    num_q = 20
    mapped_cir = QuantumCircuit(num_q)
    canonical_register = QuantumRegister(num_q, "q")
    circuit = QuantumCircuit.from_qasm_file(path)
    dag = circuit_to_dag(circuit)
    count_cx = dag._op_names['cx']
    _bit_indices = {bit: idx for idx, bit in enumerate(canonical_register)}

    tree = Tree(dag, width, depth, coupling_map, canonical_register, mapped_cir)

    count = 0

    x = []
    y = []

    while tree.nodes[tree.root_node].front_layer: # and not Tree.find_fast_path:
        tree.expansion(coupling_map, _bit_indices, score_layer, dag, canonical_register)

        if not tree.find_fast_path:
            tree.selection()

            tree.decision()
        count += 1
        print(str(count) + ': ' + str(tree.nodes[tree.root_node].add_gates) + ', ' + str(count_cx - tree.nodes[tree.root_node].exe_gates), end='\r', flush=True)
        # print(tree.nodes[tree.root_node].front_layer)
        if tree.nodes[tree.root_node].add_gates not in x:
            x.append(tree.nodes[tree.root_node].add_gates)
            y.append(count_cx - tree.nodes[tree.root_node].exe_gates) 
    print()

    plt.plot(x, y)
    plt.show()

    # tree.nodes[tree.root_node].mapped_cir.draw(output='mpl', filename='result3.pdf')
    print(tree.nodes[tree.root_node].add_gates)# 

    

    