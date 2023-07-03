from Ui_A import Ui_MainWindow as A_Ui
from Ui_B import Ui_Dialog as B_Ui
from Ui_C import Ui_Dialog as C_Ui

from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QFileDialog
from functools import partial
import sys

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
import os
from decimal import Decimal

from tree import Tree, Node


edges = [(0,1),(1,0),(1,2),(2,1),(2,3),(3,2),(3,4),(4,3),(0,5),(5,0),(1,6),(6,1),(2,7),(7,2),(3,8),(8,3),(4,9),(9,4),(1,7),(7,1),(2,6),(6,2),(3,9),(9,3),(4,8),(8,4),(5,6),(6,5),(6,7),(7,6),(7,8),(8,7),(8,9),(9,8),(5,10),(10,5),(6,11),(11,6),(7,12),(12,7),(8,13),(13,8),(9,14),(14,9),(5,11),(11,5),(6,10),(10,6),(7,13),(13,7),(8,12),(12,8),(10,11),(11,10),(11,12),(12,11),(12,13),(13,12),(13,14),(14,13),(10,15),(15,10),(11,16),(16,11),(12,17),(17,12),(13,18),(18,13),(14,19),(19,14),(11,17),(17,11),(12,16),(16,12),(13,19),(19,13),(14,18),(18,14),(15,16),(16,15),(16,17),(17,16),(17,18),(18,17),(18,19),(19,18)]


class AUi(QtWidgets.QMainWindow, A_Ui):
    def __init__(self):
        super(AUi, self).__init__()
        self.setupUi(self)


class BUi(QtWidgets.QDialog, B_Ui):
    def __init__(self):
        super(BUi, self).__init__()
        self.setupUi(self)
        self.label_3.setVisible(True)
        self.label_2.setVisible(False)
        self.label_7.setVisible(False)

    def openFile(self):
        #其中self指向自身，"读取文件夹"为标题名，"./"为打开时候的当前路径
        directory1, _ = QFileDialog.getOpenFileName(self,
                                                "选取文件",
                                                "./",
                                                "*.txt",)  # 起始路径
        self.file_path = directory1
        self.file_dic = os.path.dirname(os.path.abspath(directory1))

class CUi(QtWidgets.QDialog, C_Ui):
    def __init__(self):
        super(CUi, self).__init__()
        self.setupUi(self)
        self.file_path = None
        self.file_dic = None
        self.width = 8
        self.depth = 17
        self.score_layer = 7
        self.decay = 1

    def openFile(self):
        #其中self指向自身，"读取文件夹"为标题名，"./"为打开时候的当前路径
        directory1, _ = QFileDialog.getOpenFileName(self,
                                                "选取文件",
                                                "./")  # 起始路径
        self.file_path = directory1
        self.file_dic = os.path.dirname(os.path.abspath(directory1))
        self.textBrowser_2.setMarkdown(self.file_path)
        print(self.file_path)


    def pluswidth(self):
        if self.width < 10:
            self.width += 1
            self.lineEdit.setText(str(self.width))

    def minuswidth(self):
        if self.width > 1:
            self.width -= 1
            self.lineEdit.setText(str(self.width))
    def plusdepth(self):
        if self.depth < 30:
            self.depth += 1
            self.lineEdit_2.setText(str(self.depth))

    def minusdepth(self):
        if self.depth > 1:
            self.depth -= 1
            self.lineEdit_2.setText(str(self.depth))

    def plusscorelayer(self):
        if self.score_layer < 10:
            self.score_layer += 1
            self.lineEdit_3.setText(str(self.score_layer))

    def minusscorelayer(self):
        if self.score_layer > 1:
            self.score_layer -= 1
            self.lineEdit_3.setText(str(self.score_layer))
    def plusdecay(self):
        decay = float(self.lineEdit_4.text())
        if decay < 1:
            self.decay = float(Decimal(decay) + Decimal(0.1))
            self.lineEdit_4.setText(str(round(self.decay, 1)))

    def minusdecay(self):
        decay = float(self.lineEdit_4.text())
        if decay > 0.1:
            self.decay = float(Decimal(decay) - Decimal(0.1))
            self.lineEdit_4.setText(str(round(self.decay, 1)))

    def clearAll(self):
        self.lineEdit_5.setText('0')
        self.lineEdit_6.setText('0')
        self.textBrowser.setMarkdown('')
        self.progressBar.setValue(0)
        # self.lineEdit_7.setText('')

    def QCtransformation(self, edges, num_q, is_draw):
        mapped_cir = QuantumCircuit(num_q)
        canonical_register = QuantumRegister(num_q, "q")
        coupling_map = CouplingMap(couplinglist=edges)
        circuit = QuantumCircuit.from_qasm_file(self.file_path)
        dag = circuit_to_dag(circuit)
        count_cx = dag._op_names['cx']
        _bit_indices = {bit: idx for idx, bit in enumerate(canonical_register)}

        tree = Tree(dag, self.width, self.depth, coupling_map, canonical_register, mapped_cir, is_draw)
        count = 0
        t_start = tm.time()


        while tree.nodes[tree.root_node].front_layer: # and not Tree.find_fast_path:
            tree.expansion(coupling_map, _bit_indices, self.score_layer, dag, canonical_register, self.decay)
            if not tree.find_fast_path:
                tree.selection()
                tree.decision()
            count += 1
            print(str(count) + ': ' + str(tree.nodes[tree.root_node].add_gates) + ', ' + str(count_cx - tree.nodes[tree.root_node].exe_gates) + ', ' + str(tree.nodes[0].swap_gates), end='\r', flush=True)
            value = tree.nodes[tree.root_node].exe_gates / count_cx
            self.progressBar.setValue(int(value * 100))
        t_end = tm.time()

        result_gate = tree.nodes[tree.root_node].add_gates
        result_time = t_end - t_start
        if not self.checkBox.isChecked():
            self.textBrowser.setMarkdown(self.file_dic + '\output.qasm')
        self.lineEdit_6.setText(str(result_gate))
        self.lineEdit_5.setText(str(result_time) + 's')
        tree.nodes[tree.root_node].mapped_cir.qasm(filename=(self.file_dic + '\output.qasm'))

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    a = AUi()
    a.show()
    b = BUi()
    c = CUi()
    # button是你定义的按钮
    a.pushButton.clicked.connect(
    	lambda:{a.close(), b.show()}
   	)
    a.pushButton_2.clicked.connect(
        lambda:{a.close()}
    )
    b.pushButton.clicked.connect(
        lambda:{b.close(), c.show()}
    )
    b.pushButton_2.clicked.connect(
        lambda:{b.openFile()}
    )
    b.pushButton_4.clicked.connect(
        lambda:{b.label_2.setVisible(True), b.label_3.setVisible(False), b.label_7.setVisible(False)}
    )
    b.pushButton_3.clicked.connect(
        lambda:{b.label_2.setVisible(False), b.label_3.setVisible(True), b.label_7.setVisible(False)}
    )
    b.pushButton_5.clicked.connect(
        lambda:{b.label_2.setVisible(False), b.label_3.setVisible(False), b.label_7.setVisible(True)}
    )
    c.pushButton_2.clicked.connect(
        lambda:{c.close(), b.show()}
    )
    c.pushButton.clicked.connect(
        lambda:{c.close()}
    )
    c.pushButton_8.clicked.connect(
        lambda:{c.openFile(), c.clearAll()}
    )
    c.pushButton_7.clicked.connect(
        lambda:{c.QCtransformation(edges, 20, not c.checkBox.isChecked())}
    )
    c.pushButton_3.clicked.connect(
        lambda:{c.pluswidth()}
    )
    c.pushButton_4.clicked.connect(
        lambda:{c.minuswidth()}
    )
    c.pushButton_5.clicked.connect(
        lambda:{c.plusdepth()}
    )
    c.pushButton_6.clicked.connect(
        lambda:{c.minusdepth()}
    )
    c.pushButton_12.clicked.connect(
        lambda:{c.plusscorelayer()}
    )
    c.pushButton_11.clicked.connect(
        lambda:{c.minusscorelayer()}
    )
    c.pushButton_10.clicked.connect(
        lambda:{c.plusdecay()}
    )
    c.pushButton_9.clicked.connect(
        lambda:{c.minusdecay()}
    )

    sys.exit(app.exec_())
