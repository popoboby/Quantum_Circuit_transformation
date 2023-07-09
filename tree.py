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

# dag
# front_layer
# initial_layout 
# applied_predecessors
# canonical_register
# coupling_map
# mapped_cir


class Node:
    def __init__(self, front_layer, initial_layout, mapped_cir, applied_predecessors):
        # self.id = 0
        self.score = 0
        self.front_layer = front_layer.copy()
        self.initial_layout = initial_layout.copy()
        self.applied_predecessors = applied_predecessors.copy()
        self.successors = [] # 记录该节点的子孙系欸但
        self.mapped_cir = mapped_cir
        self.father_node = -1
        self.add_gates = 0
        self.exe_gates = 0
        self.swap_gates = None
        self.exe_score = 0
        self.exe_num = 0

    # def _apply_gate(mapped_cir, node, current_layout, canonical_register): # 应用门
    #     new_node = _transform_gate_for_layout(node, current_layout, canonical_register)
    #     mapped_cir.append(new_node.op, new_node.qargs, new_node.cargs)

    def _successors(self, node, dag):
        for _, successor, edge_data in dag.edges(node):
            if type(successor) != DAGOpNode:
                continue
            if isinstance(edge_data, Qubit):
                yield successor

    def _is_resolved(self, node):
        return self.applied_predecessors[node] == len(node.qargs)

    def _transform_gate_for_layout(self, op_node, layout, device_qreg): # 未知
        """Return node implementing a virtual op on given layout."""
        mapped_op_node = copy(op_node)
        premap_qargs = op_node.qargs
        mapped_qargs = map(lambda x: device_qreg[layout[x]], premap_qargs)
        mapped_op_node.qargs = list(mapped_qargs)
        return mapped_op_node

    def _apply_gate(self, node, canonical_register):
        new_node = self._transform_gate_for_layout(node, self.initial_layout, canonical_register)
        self.mapped_cir.append(new_node.op, new_node.qargs, new_node.cargs)

    def is_exp(self): # 判断当前节点是否能够扩展
        if self.front_layer:
            return True
        else:
            return False

    def copy(self, node):
        node.score = 0
        node.front_layer = self.front_layer.copy()
        node.initial_layout = self.initial_layout.copy()
        node.applied_predecessors = self.applied_predecessors.copy()
        # node.successors = [] # 记录该节点的子孙节点
        node.mapped_cir = self.mapped_cir
        node.father_node = -1
        node.add_gates = 0


    def get_pertinent_swaps(self, coupling_map, _bit_indices, score_layer, dag): # 得到当前选择节点的最前层的相关交换门
        candidate_swaps = [] # 一个集合变量
        scores = []
        add_gates = []
        for node in self.front_layer: # 遍历最前层的节点
            for virtual in node.qargs: # 节点的门操纵 
                physical = self.initial_layout[virtual] # 物理比特 （关注一下）
                K = []
                for neighbor in coupling_map.neighbors(physical):
                    K.append(neighbor)

                K.sort()
                for neighbor in K: # 只有相邻的物理比特才可以执行
                    virtual_neighbor = self.initial_layout[neighbor]
                    swap = sorted([virtual, virtual_neighbor], key=lambda q: _bit_indices[q])
                    candidate_swaps.append(tuple(swap))

        valid_swaps = []
        for swap in candidate_swaps:
            trial_layout = self.initial_layout.copy()
            trial_layout.swap(*swap)
            cost = 0
            for node in self.front_layer:
                v0, v1 = node.qargs
                cost += (coupling_map.distance(self.initial_layout[v0], self.initial_layout[v1]) - coupling_map.distance(trial_layout[v0], trial_layout[v1]))

            if cost < 0:
                # candidate_swaps.remove(swap)
                continue

            valid_swaps.append(swap)

        
            score, add_gate = self.get_scores(coupling_map, trial_layout, score_layer, dag, 3)
            scores.append(score)
            add_gates.append(add_gate)

        return valid_swaps, scores, add_gates

    def execute_gates(self, coupling_map, dag, canonical_register, is_draw): # 该节点一口气执行到front_layer不再更新为止
        # canonical_register = QuantumRegister(num_q, 'q')
        # _bit_indices = {bit: idx for idx, bit in enumerate(canonical_register)}
        flag = True
        score = 0
        while flag:
            flag = False
            # applied_predecessors = defaultdict(int)
            execute_gate_list = []
            for node in self.front_layer:
                if node.name == 'cx':
                    v0, v1 = node.qargs
                    if coupling_map.graph.has_edge(self.initial_layout[v0], self.initial_layout[v1]):
                        score += 1
                        execute_gate_list.append(node)
                else:
                    execute_gate_list.append(node)

            if execute_gate_list:
                flag = True
                for node in execute_gate_list:
                    if is_draw:
                        self._apply_gate(node, canonical_register)
                    self.front_layer.remove(node)
                    
                    for successor in self._successors(node, dag):
                        self.applied_predecessors[successor] += 1

                        if self._is_resolved(successor):
                            self.front_layer.append(successor)

        # print(self.front_layer)
        self.score = score
        self.exe_gates = score
        self.exe_num = score
        # return score

    def get_scores(self, coupling_map, trial_layout, score_layer, dag, EXTENDED_SET_SIZE):
        cost = 0.0

        involved_nodes = self.front_layer.copy()
        # for node in self.front_layer:
        #     v0, v1 = node.qargs
        #     cost += (coupling_map.distance(current_layout[v0], current_layout[v1]) - coupling_map.distance(trial_layout[v0], trial_layout[v1]))

        # if cost < 0: # 当前操作无法使第一层节点的距离之差之和 > 0
        #     # 不再执行后续操作，直接返回值
        #     return cost
        
        cost = 0.0
        # 如果该交换可以使得第一层之和值大于0，则继续后面五层的操作
        score_add = 0.0

        added_swaps = 1
        isPos = True
        decay = 1
        i = 0
        for layer_i in range(score_layer): # 往后看score_layer层
            i_max = len(involved_nodes)
            flag = False
            while i <= i_max - 1:
                node = involved_nodes[i]
                i += 1
                if len(node.qargs) == 1:
                    for successor in self._successors(node, dag):
                        if not successor in involved_nodes:
                            involved_nodes.append(successor)
                    continue
                flag = True
                v0, v1 = node.qargs
                score_add = coupling_map.distance(self.initial_layout[v0], self.initial_layout[v1]) - coupling_map.distance(trial_layout[v0], trial_layout[v1])
                if score_add > 0 and isPos:
                    isPos = True

                elif score_add > 0 and not isPos:
                    isPos = True
                    added_swaps += 1

                elif score_add < 0 and isPos:
                    added_swaps += 1
                    score_add = 0
                    isPos = False

                elif isPos == False and score_add < 0:
                    score_add = 0

                cost += score_add * decay
                for successor in self._successors(node, dag):
                    if not successor in involved_nodes:
                        involved_nodes.append(successor)
            if flag:
                decay *= 0.7
            else:
                layer_i -= 1
            
        return cost, added_swaps


    def swap(self, swap, coupling_map, canonical_register, dag, is_draw): # 执行完swap操作，并把front_layer推到最后一层，并计算其得分
        swap_node = DAGOpNode(op=SwapGate(), qargs=swap)
        if is_draw:
            self._apply_gate(swap_node, canonical_register)
        self.initial_layout.swap(*swap)
        self.execute_gates(coupling_map, dag, canonical_register, is_draw)


            
        

class Tree:
    def __init__(self, dag, width, depth, coupling_map, canonical_register, mapped_cir, is_draw):

        self.node_count = 0
        self.width = width
        self.depth = depth
        self.selec_count = 0
        self.nodes = []
        applied_predecessors = defaultdict(int)
        self.root_node = 0
        self.is_draw = is_draw
        # self.root_node.execute_gates()
        initial_mapping = [i for i in range(len(dag.qubits))]
        initial_layout = Layout({q: dag.qubits[i] for i, q in enumerate(initial_mapping)})

        front_layer = dag.front_layer()
        applied_predecessors = defaultdict(int)
        for _, input_node in dag.input_map.items():
            for successor in self._successors(input_node, dag):
                applied_predecessors[successor] += 1

        new_node = Node(front_layer, initial_layout, mapped_cir, applied_predecessors)
        new_node.execute_gates(coupling_map, dag, canonical_register, self.is_draw)
        self.add_node(new_node)
        self.selec_nodes = [] # 待选择节点
        self.selec_nodes
        self.exp_nodes = [] # 待拓展节点
        self.exp_nodes.append(self.root_node)
        # self.dec_nodes = [] # 待决策节点
        # self.first_exp_layer_num = 0
        # self.add_gates = 0
        self.find_fast_path = False

    def _successors(self, node, dag):
        for _, successor, edge_data in dag.edges(node):
            if type(successor) != DAGOpNode:
                continue
            if isinstance(edge_data, Qubit):
                yield successor

    def get_best_nodes(self, nodes_id1, count): # 根据nodes_id和count返回一个或多个最优节点
        nodes_id = nodes_id1.copy()

        # 考虑两种情况
        if count >= len(nodes_id):
            return nodes_id
        else:

            seed = np.random.randint(0, np.iinfo(np.int32).max)
            rng = np.random.default_rng(seed)

            target_nodes = []

            
            node_scores = dict.fromkeys(nodes_id, 0)
            for node in node_scores:
                node_scores[node] = self.nodes[node].score
            i = 0
            while i < count:
                i += 1
                max_score = max(node_scores.values())
                best_nodes = [k for k, v in node_scores.items() if v == max_score]
                best_node = min(best_nodes) # rng.choice(best_nodes)
                target_nodes.append(best_node)
                node_scores.pop(best_node)

            return target_nodes






    def add_node(self, node):
        node.id = self.node_count
        self.node_count += 1
        self.nodes.append(node)


    def expansion(self, coupling_map, _bit_indices, score_layer, dag, canonical_register):

        self.selec_nodes = []

        end_nodes = []
        
        for node_id in self.exp_nodes: # 从遍历待拓展节点

        
            swaps, scores, add_gates = self.nodes[node_id].get_pertinent_swaps(coupling_map, _bit_indices, score_layer, dag)

        
        
            for swap, score, add_gate in zip(swaps, scores, add_gates):
                new_node = Node(self.nodes[node_id].front_layer.copy(), self.nodes[node_id].initial_layout.copy(), self.nodes[node_id].mapped_cir.copy(), self.nodes[node_id].applied_predecessors.copy())
                
                new_node.swap(swap, coupling_map, canonical_register, dag, self.is_draw) # 该节点执行交换操作
                if self.nodes[node_id].exe_num == 0 and swap == self.nodes[node_id].swap_gates:
                    continue
                new_node.exe_gates += self.nodes[node_id].exe_gates
                new_node.score = (new_node.score / add_gate) + score + self.nodes[node_id].score
                new_node.add_gates = self.nodes[node_id].add_gates + 1
                new_node.swap_gates = swap
                new_node.exe_score = new_node.exe_gates - self.nodes[self.root_node].exe_gates
                
                

                if self.selec_count == 0:
                    new_node.father_node = self.node_count
                if self.selec_count > 0:
                    new_node.father_node = self.nodes[node_id].father_node
                self.add_node(new_node)
                if not new_node.is_exp(): # 判断是否扩展出的节点已经到底
                    end_nodes.append(self.node_count - 1)
                # self.nodes[node_id].successors.append(self.node_count - 1)
                self.selec_nodes.append(self.node_count - 1)

        if end_nodes:
            self.find_fast_path = True
            root_node = self.get_best_nodes(end_nodes, 1)
            self.root_node = root_node[0]
            return 

    def get_best_nodes1(self, nodes_id1, count): # 根据nodes_id和count返回一个或多个最优节点
        nodes_id = nodes_id1.copy()

        # 考虑两种情况
        if count >= len(nodes_id):
            return nodes_id
        else:

            seed = np.random.randint(0, np.iinfo(np.int32).max)
            rng = np.random.default_rng(seed)

            target_nodes = []

            
            node_scores = dict.fromkeys(nodes_id, 0)
            for node in node_scores:
                node_scores[node] = self.nodes[node].exe_score
            i = 0
            while i < count:
                i += 1
                max_score = max(node_scores.values())
                best_nodes = [k for k, v in node_scores.items() if v == max_score]
                best_node = min(best_nodes) # rng.choice(best_nodes)
                target_nodes.append(best_node)
                node_scores.pop(best_node)

            return target_nodes

    def selection(self):
        self.exp_nodes = self.get_best_nodes(self.selec_nodes, self.width)
        self.selec_count += 1
        # 剪枝策略：
        father_node_id = self.nodes[self.exp_nodes[0]].father_node
        if father_node_id == -1:
            return
        count = 0
        for node_id in self.exp_nodes:
            if father_node_id == self.nodes[node_id].father_node:
                count += 1

        if count == len(self.exp_nodes):
            self.selec_count = self.depth

    def decision(self):
        if self.selec_count == self.depth:
            dec_node_id = self.get_best_nodes1(self.exp_nodes, 1)
            root_node_id = self.get_father_node(dec_node_id[0])
            self.root_node = root_node_id
            self.exp_nodes = []
            self.exp_nodes.append(0)
            self.selec_count = 0
            self.delete_nodes()

    def get_father_node(self, dec_node):
        for node in self.nodes:
            if node.id == dec_node:
                return node.father_node


    def delete_nodes(self): # 删除除根节点以外的全部节点

        node = Node(self.nodes[self.root_node].front_layer, self.nodes[self.root_node].initial_layout, self.nodes[self.root_node].mapped_cir, self.nodes[self.root_node].applied_predecessors)
        node.id = 0
        node.add_gates = self.nodes[self.root_node].add_gates
        node.score = self.nodes[self.root_node].score
        node.exe_gates = self.nodes[self.root_node].exe_gates
        node.swap_gates = self.nodes[self.root_node].swap_gates
        node.exe_num = self.nodes[self.root_node].exe_num
        node.exe_score = 0
        self.nodes.clear()
        self.nodes.append(node)
        self.root_node = 0
        self.node_count = 1
            



        

            
            
            
            





