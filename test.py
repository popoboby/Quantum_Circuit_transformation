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
import time as tm

from tree import Tree, Node

test_files = ['3_17_13.qasm', 'rd84_142.qasm', 'sym6_145.qasm', 'cycle10_2_110.qasm', 'radd_250.qasm', 'rd73_252.qasm', 'adr4_197.qasm', 'z4_268.qasm', 'misex1_241.qasm', 'sqn_258.qasm', 'square_root_7.qasm', 'rd84_253.qasm']

# 你好



path = 'E:/files/quant/code/quantum/qasm_circuits/'

width = 7
depth = 17
score_layer = 7 # 1020 rd84_253.qasm # decay = 0.7 7
is_draw = True

# rd84_253: 923 # score_layer==7
# sym6_145: 169
# cycle10_2_110: 358

# rd84_253: 969 score_layer==8
# sym6_145: 235
# cycle10_2_110: 469
# radd_250: 264
# rd84_142: 42
# adr4_197: 212
# sqn_258: //

# rd84_142: 39 
# adr4_197: 250
# sqn_258: 682 
# square_root_7: 422
# rd73_252: 398 
# radd_250: 263 
# rd84_253: 1031 
# z4_268: //
# sym6_145: 211 
# cycle10_2_110: 472
# misex1_241: 266



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

# rd84_142: 38  +4 # 给加入门数赋上2的权重
# adr4_197: 282 -29
# sqn_258: 638 +105
# square_root_7: 432 +65
# rd73_252: 448 +75
# radd_250: 282 -42
# rd84_253: 976 +64
# z4_268: 238 +20
# sym6_145: 203 +123
# cycle10_2_110: 475 -105
# misex1_241: 256 -3

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





if __name__ == '__main__':



    edges = [(0,1),(1,0),(1,2),(2,1),(2,3),(3,2),(3,4),(4,3),(0,5),(5,0),(1,6),(6,1),(2,7),(7,2),(3,8),(8,3),(4,9),(9,4),(1,7),(7,1),(2,6),(6,2),(3,9),(9,3),(4,8),(8,4),(5,6),(6,5),(6,7),(7,6),(7,8),(8,7),(8,9),(9,8),(5,10),(10,5),(6,11),(11,6),(7,12),(12,7),(8,13),(13,8),(9,14),(14,9),(5,11),(11,5),(6,10),(10,6),(7,13),(13,7),(8,12),(12,8),(10,11),(11,10),(11,12),(12,11),(12,13),(13,12),(13,14),(14,13),(10,15),(15,10),(11,16),(16,11),(12,17),(17,12),(13,18),(18,13),(14,19),(19,14),(11,17),(17,11),(12,16),(16,12),(13,19),(19,13),(14,18),(18,14),(15,16),(16,15),(16,17),(17,16),(17,18),(18,17),(18,19),(19,18)]

    result_gates = []
    result_time = []
    num_q = 20
    mapped_cir = QuantumCircuit(num_q)
    canonical_register = QuantumRegister(num_q, "q")
    coupling_map = CouplingMap(couplinglist=edges)

    for test in test_files:
        print(test)
        print('width:', width)
        print('depth:', depth)
        print('score_layer:', score_layer)
        circuit = QuantumCircuit.from_qasm_file(path + test)
        dag = circuit_to_dag(circuit)
        count_cx = dag._op_names['cx']
        _bit_indices = {bit: idx for idx, bit in enumerate(canonical_register)}

        tree = Tree(dag, width, depth, coupling_map, canonical_register, mapped_cir, is_draw)

        count = 0

        x = []
        y = []

        t_start = tm.time()

        t_exp = 0.0
        t_exp_1 = 0.0
        t_exp_2 = 0.0
        t_sel = 0.0
        t_dec = 0.0
        total = 0.0
        while tree.nodes[tree.root_node].front_layer: # and not Tree.find_fast_path:
            total_start = tm.time()
            t_exp_start = tm.time()
            t_1, t_2 = tree.expansion(coupling_map, _bit_indices, score_layer, dag, canonical_register)
            t_exp_1 += t_1
            t_exp_2 += t_2
            t_exp_end = tm.time()

            t_exp += (t_exp_end - t_exp_start)

            if not tree.find_fast_path:
                t_sel_start = tm.time()
                tree.selection()
                t_sel_end = tm.time()
                t_sel += (t_sel_end - t_sel_start)


                t_dec_start = tm.time()
                tree.decision()
                t_dec_end = tm.time()
                t_dec += (t_dec_end - t_dec_start)

            count += 1
            total_end = tm.time()
            total += (total_end - total_start)
            print(str(count) + ': ' + str(tree.nodes[tree.root_node].add_gates) + ', ' + str(count_cx - tree.nodes[tree.root_node].exe_gates) + ', ' + str(tree.nodes[0].swap_gates), end='\r', flush=True)
        t_end = tm.time()

        result_gates.append(tree.nodes[tree.root_node].add_gates)
        result_time.append(t_end - t_start)

        print()
        print(t_end - t_start)
        print(total)
        print('expand time spend: ' + str(t_exp))
        print('expand_1 time spend: ' + str(t_exp_1))
        print('expand_2 time spend: ' + str(t_exp_2))
        print('select time spend: ' + str(t_sel))
        print('decision time spend: ' + str(t_dec))
        # plt.plot(x, y)
        # plt.show()

        if is_draw:
            tree.nodes[tree.root_node].mapped_cir.draw(output='mpl', filename='result3.pdf')
            tree.nodes[tree.root_node].mapped_cir.qasm(filename='test1.qasm')
        print(tree.nodes[tree.root_node].add_gates)

    i = 0
    for res in test_files:
        print(res)
        print('add_gates: ' + str(result_gates[i]))
        print('spend_time: ' + str(result_time[i]))
        i += 1