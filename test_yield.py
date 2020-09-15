# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 23:31:07 2020

@author: shiba
"""

import numpy as np       ### 今回はnumpy配列使うので必要
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

hold_yield = pd.read_csv('C:/Users/shiba/Desktop/Python/Yield_test.csv')

### 3Dグラフの箱を用意(多分)
fig = plt.figure()
ax = Axes3D(fig)

### y軸の数値を入れる配列
y_num = np.array([])
### y軸に設定したい文字列を入れる配列
y_axis = np.array([])

for i in range(4):
    x = hold_yield['Maturity'][(i*10):(i*10+10)]
    y = hold_yield['Date2'][(i*10):(i*10+10)]
    z = hold_yield['Yield'][(i*10):(i*10+10)]
    # 3Dでプロット
    ax.plot(x, y, z)
    y_axis = np.append(y_axis,hold_yield['Date'][i*10])
    y_num = np.append(y_num,hold_yield['Date2'][i*10])

plt.yticks(y_num,y_axis)
# 軸ラベル
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
# 表示
plt.show()


### テストⅡ
hold_yield2 = pd.read_csv('C:/Users/shiba/Desktop/Python/spot-rate(JGB).csv')
    
fig = plt.figure() #プロット領域の作成
ax = fig.gca(projection='3d') #プロット中の軸の取得。gca は"Get Current Axes" の略。

y = hold_yield2['date2'] # x点として[-2, 2]まで0.05刻みでサンプル
x = np.arange(0.5, 20.5, 0.5)  # y点として[-2, 2]まで0.05刻みでサンプル
#x, y = np.meshgrid(x, y)  # 上述のサンプリング点(x,y)を使ったメッシュ生成
X, Y = np.meshgrid(x, y)
Z =  hold_yield2.loc[:, '0.5':'20']

y_axis = hold_yield2['date'].fillna(' ')
y_num = hold_yield2['date2']

plt.yticks(y_num,y_axis)
ax.set_xticks(np.arange(5,25 , step=5))
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.plot_wireframe(X, Y, Z, color='blue',linewidth=1) # ワイヤーフレームのプロット。linewidthは曲面のメッシュの線の太さ，をそれぞれ表す。

plt.show() # 絵の出力




