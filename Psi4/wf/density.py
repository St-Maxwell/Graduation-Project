import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in' 
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams["xtick.top"] = True
mpl.rcParams["ytick.right"] = True
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['mathtext.fontset'] = 'custom'
import numpy as np
fig, ax = plt.subplots(tight_layout=True)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
ax.set_xlabel(r'$y/\mathrm{Angstroms}$', fontsize=16)
ax.set_ylabel(r'$z/\mathrm{Angstroms}$', fontsize=16)
ax.set_xlim(-2.2, 2.2)
ax.set_ylim(-3.3, 3.3)
ax.grid(linestyle=':')
ax.set_facecolor('whitesmoke')
data = np.loadtxt('NCO-density.txt')
val = data[:,3]
X = np.linspace(-2.38130, 2.35748, 200)
Y = np.linspace(-3.55675, 3.57744, 200)
z = np.zeros((200, 200))
k = 0
for i in range(200):
  for j in range(200):
    z[j,i] = val[k]
    k += 1

line_lv = [0.001, 0.002, 0.004, 0.008, 0.020, 0.040, 0.080, 0.2, 0.4]
lb_lv = [0.001, 0.002, 0.004, 0.008, 0.020, 0.040, 0.080, 0.2, 0.4]
lv_fmt = {0.001:'0.001', 0.002:'0.002', 0.004:'0.004', 0.008:'0.008',
          0.020:'0.02', 0.040:'0.04', 0.080:'0.08', 0.2:'0.2', 0.4:'0.4'}
C = plt.contour(X, Y, z, levels=line_lv, colors='black')
plt.clabel(C, lb_lv, inline=True, fontsize=12, inline_spacing=1., fmt=lv_fmt)
plt.text(0.0, 1.22, 'N', fontsize=15, horizontalalignment='center', verticalalignment='center')
plt.text(0.0, -0.01, 'C', fontsize=15, horizontalalignment='center', verticalalignment='center')
plt.text(0.0, -1.19, 'O', fontsize=15, horizontalalignment='center', verticalalignment='center')
circle_N = mpl.patches.Circle((0.0, 1.231991), radius=0.2, facecolor=(.8,.8,1.,1.), edgecolor='black', linewidth=.8)
ax.add_artist(circle_N)
circle_C = mpl.patches.Circle((0.0, 0.0), radius=0.2, facecolor=(.8,.8,1.,1.), edgecolor='black', linewidth=.8)
ax.add_artist(circle_C)
circle_O = mpl.patches.Circle((0.0, -1.175455), radius=0.2, facecolor=(.8,.8,1.,1.), edgecolor='black', linewidth=.8)
ax.add_artist(circle_O)
fig.set_size_inches(4, 6)
figure_fig = plt.gcf()
figure_fig.savefig('NCO-density.pdf', format='pdf', dpi=1000)