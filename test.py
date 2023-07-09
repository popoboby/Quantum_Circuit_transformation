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

test_files0 = ['cycle10_2_110.qasm']

test_files1 = ['rd84_142.qasm',
        'adr4_197.qasm',
        'radd_250.qasm',
        'z4_268.qasm',
        'sym6_145.qasm',
        'misex1_241.qasm',
        'rd73_252.qasm',
        'cycle10_2_110.qasm',
        'square_root_7.qasm',
        'sqn_258.qasm',
        'rd84_253.qasm']

test_files2 = ['graycode6_47.qasm',
        'xor5_254.qasm',
        'ex1_226.qasm',
        '4gt11_84.qasm',
        'ex-1_166.qasm',
        'ham3_102.qasm',
        '4mod5-v0_20.qasm',
        '4mod5-v1_22.qasm',
        'mod5d1_63.qasm',
        '4gt11_83.qasm',
        '4gt11_82.qasm',
        'rd32-v0_66.qasm',
        'mod5mils_65.qasm',
        '4mod5-v0_19.qasm',
        'rd32-v1_68.qasm',
        'alu-v0_27.qasm',
        '3_17_13.qasm',
        '4mod5-v1_24.qasm',
        'alu-v1_29.qasm',
        'alu-v1_28.qasm',
        'alu-v3_35.qasm',
        'alu-v2_33.qasm',
        'alu-v4_37.qasm',
        'miller_11.qasm',
        'decod24-v0_38.qasm',
        'alu-v3_34.qasm',
        'decod24-v2_43.qasm',
        'mod5d2_64.qasm',
        '4gt13_92.qasm',
        '4gt13-v1_93.qasm',
        'one-two-three-v2_100.qasm',
        '4mod5-v1_23.qasm',
        '4mod5-v0_18.qasm',
        'one-two-three-v3_101.qasm',
        '4mod5-bdd_287.qasm',
        'decod24-bdd_294.qasm',
        '4gt5_75.qasm',
        'alu-v0_26.qasm',
        'rd32_270.qasm',
        'alu-bdd_288.qasm',
        'decod24-v1_41.qasm',
        '4gt5_76.qasm',
        '4gt13_91.qasm',
        '4gt13_90.qasm',
        'alu-v4_36.qasm',
        '4gt5_77.qasm',
        'one-two-three-v1_99.qasm',
        'rd53_138.qasm',
        'one-two-three-v0_98.qasm',
        '4gt10-v1_81.qasm',
        'decod24-v3_45.qasm',
        'aj-e11_165.qasm',
        '4mod7-v0_94.qasm',
        'alu-v2_32.qasm',
        '4mod7-v1_96.qasm',
        'cnt3-5_179.qasm',
        'mod10_176.qasm',
        '4gt4-v0_80.qasm',
        '4gt12-v0_88.qasm',
        '0410184_169.qasm',
        '4_49_16.qasm',
        '4gt12-v1_89.qasm',
        '4gt4-v0_79.qasm',
        'hwb4_49.qasm',
        '4gt4-v0_78.qasm',
        'mod10_171.qasm',
        '4gt12-v0_87.qasm',
        '4gt12-v0_86.qasm',
        '4gt4-v0_72.qasm',
        '4gt4-v1_74.qasm',
        'mini-alu_167.qasm',
        'one-two-three-v0_97.qasm',
        'rd53_135.qasm',
        'ham7_104.qasm',
        'decod24-enable_126.qasm',
        'mod8-10_178.qasm',
        '4gt4-v0_73.qasm',
        'ex3_229.qasm',
        'mod8-10_177.qasm',
        'alu-v2_31.qasm',
        'C17_204.qasm',
        'rd53_131.qasm',
        'alu-v2_30.qasm',
        'mod5adder_127.qasm',
        'rd53_133.qasm',
        'majority_239.qasm',
        'ex2_227.qasm',
        'cm82a_208.qasm',
        'sf_276.qasm',
        'sf_274.qasm',
        'con1_216.qasm',
        'rd53_130.qasm',
        'f2_232.qasm',
        'rd53_251.qasm',
        'hwb5_53.qasm',
        'radd_250.qasm',
        'rd73_252.qasm',
        'cycle10_2_110.qasm',
        'hwb6_56.qasm',
        'cm85a_209.qasm',
        'rd84_253.qasm',
        'root_255.qasm',
        'mlp4_245.qasm',
        'urf2_277.qasm',
        'sym9_148.qasm',
        'hwb7_59.qasm',
        'clip_206.qasm',
        'sym9_193.qasm',
        'dist_223.qasm',
        'sao2_257.qasm',
        'urf5_280.qasm',
        'urf1_278.qasm',
        'sym10_262.qasm',
        'hwb8_113.qasm',
        ]



benchmark = 'B11'

if benchmark == 'B11':
    path = 'E:/files/quant/code/quantum/qasm_circuits/'
    test_files = test_files0
elif benchmark == 'B114':
    path = 'E:/files/quant/code/quantum/qasm_circuits_new/'
    test_files = test_files2


width = 8 # 8
depth = 15 # 15
score_layer = 7 # 1020 rd84_253.qasm # decay = 0.7 7
is_draw = False


if __name__ == '__main__':



    edges = [(0,1),(1,0),(1,2),(2,1),(2,3),(3,2),(3,4),(4,3),(0,5),(5,0),(1,6),(6,1),(2,7),(7,2),(3,8),(8,3),(4,9),(9,4),(1,7),(7,1),(2,6),(6,2),(3,9),(9,3),(4,8),(8,4),(5,6),(6,5),(6,7),(7,6),(7,8),(8,7),(8,9),(9,8),(5,10),(10,5),(6,11),(11,6),(7,12),(12,7),(8,13),(13,8),(9,14),(14,9),(5,11),(11,5),(6,10),(10,6),(7,13),(13,7),(8,12),(12,8),(10,11),(11,10),(11,12),(12,11),(12,13),(13,12),(13,14),(14,13),(10,15),(15,10),(11,16),(16,11),(12,17),(17,12),(13,18),(18,13),(14,19),(19,14),(11,17),(17,11),(12,16),(16,12),(13,19),(19,13),(14,18),(18,14),(15,16),(16,15),(16,17),(17,16),(17,18),(18,17),(18,19),(19,18)]

    result_gates = []
    result_time = []
    num_q = 20
    mapped_cir = QuantumCircuit(num_q)
    canonical_register = QuantumRegister(num_q, "q")
    coupling_map = CouplingMap(couplinglist=edges)


    for i in range(2):
        if i == 0:
            width = 8 # 8
            depth = 15 # 15
            score_layer = 7 # 1020 rd84_253.qasm # decay = 0.7 7

        elif i == 1:
            width = 7 # 8
            depth = 17 # 15
            score_layer = 7 # 1020 rd84_253.qasm # decay = 0.7 7

        print('width:', width)
        print('depth:', depth)
        print('score_layer:', score_layer)
        print()

        for test in test_files:
            print(test)
        
            circuit = QuantumCircuit.from_qasm_file(path + test)
            dag = circuit_to_dag(circuit)
            count_cx = dag._op_names['cx']
            _bit_indices = {bit: idx for idx, bit in enumerate(canonical_register)}

            tree = Tree(dag, width, depth, coupling_map, canonical_register, mapped_cir, is_draw)

            count = 0

            x = []
            y = []

            t_start = tm.time()

            while tree.nodes[tree.root_node].front_layer: # and not Tree.find_fast_path:
                tree.expansion(coupling_map, _bit_indices, score_layer, dag, canonical_register)

                if not tree.find_fast_path:
                    tree.selection()

                    tree.decision()
                count += 1
                print(str(count) + ': ' + str(tree.nodes[tree.root_node].add_gates) + ', ' + str(count_cx - tree.nodes[tree.root_node].exe_gates), end='\r', flush=True)
                if tree.nodes[tree.root_node].add_gates not in x:
                    x.append(tree.nodes[tree.root_node].add_gates)
                    y.append(count_cx - tree.nodes[tree.root_node].exe_gates) 
            t_end = tm.time()
            result_gates.append(tree.nodes[tree.root_node].add_gates)
            result_time.append(t_end - t_start)

            print()
            print(t_end - t_start)
            print('total_gates: ' + str(sum(result_gates)))
            if i == 0:
                plt.plot(x, y, color='green', label='Beam_search')
            elif i == 1:
                plt.plot(x, y, color='orange', label='Monte_Carlo')
            if i == 1:
                plt.ylabel('rest_cnots')
                plt.xlabel('added_swaps')
                plt.legend()
                plt.show()

            if is_draw:
                tree.nodes[tree.root_node].mapped_cir.draw(output='mpl', filename='result3.pdf')
            print(tree.nodes[tree.root_node].add_gates)
            print()

    i = 0

    print('width:', width)
    print('depth:', depth)
    print('score_layer:', score_layer)
    print()
    for res in test_files:
        print(res)
        print('add_gates: ' + str(result_gates[i]))
        print('spend_time: ' + str(result_time[i]))
        print()
        i += 1

    print()
    print('add_gates: ' + str(sum(result_gates)))
    print('total_time: ' + str(sum(result_time)))


