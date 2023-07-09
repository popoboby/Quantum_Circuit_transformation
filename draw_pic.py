# 画出B11测试集上，MTC算法和集束搜索算法的各个线路加入swap门数的直方图
# 再画出rd84_253线路上两种方法的折现图


MTC_results = []
BS_results = []

import matplotlib.pyplot as plt
import numpy as np
x = np.arange(11)
y = [41, 255, 227, 196, 169, 240, 401, 358, 381, 635, 966]
y1 = [39, 268, 228, 234, 266, 222, 356, 455, 448, 702, 949]
# 多数据并列柱状图
bar_width = 0.25  
# tick_label=['rd84_142',
#         'adr4_197',
#         'radd_250',
#         'z4_268',
#         'sym6_145',
#         'misex1_241',
#         'rd73_252',
#         'cycle10_2_110',
#         'square_root_7',
#         'sqn_258',
#         'rd84_253']
tick_label=['c1',
        'c2',
        'c3',
        'c4',
        'c5',
        'c6',
        'c7',
        'c8',
        'c9',
        'c10',
        'c11']
# plt.bar(x, y, bar_width, align='center', color='#66c2a5', label='Monte_Carlo')
plt.bar(x, y, bar_width, align='center', color='green', label='Beam_search')
# plt.bar(x+bar_width, y1, bar_width, align='center', color='#8da0cb', label='Beam_search')
plt.bar(x+bar_width, y1, bar_width, align='center', color='orange', label='Monte_Carlo')
plt.xlabel('B11_circuits')
plt.ylabel('added_swaps')
plt.xticks(x, tick_label)
plt.legend()
plt.show()