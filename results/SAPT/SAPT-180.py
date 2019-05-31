import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
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
data = np.loadtxt('180.txt')
r = data[:,0]
elst = data[:,1]
exch = data[:,2]
ind = data[:,3]
disp = data[:,4]
sapt = data[:,5]
fig, ax = plt.subplots(tight_layout=True)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
ax.set_xlabel(r'$R/\mathrm{Angstroms}$', fontsize=16)
ax.set_ylabel(r'$E/\mathrm{cm}^{-1}$', fontsize=16)
line_elst, = ax.plot(r,elst,linestyle=(0, (3, 1, 1, 1, 1, 1)))
line_exch, = ax.plot(r,exch,linestyle=(0, (5, 1)))
line_ind, = ax.plot(r,ind,linestyle=(0, (3, 1, 1, 1)))
line_disp, = ax.plot(r,disp,linestyle=(0, (1, 1)))
line_sapt, = ax.plot(r,sapt,'-')
axins = inset_axes(ax,width='45%',height='60%',loc=1)
axins.set_xlim(4.75,5.15)
axins.set_ylim(-160,110)
line_elst_sub, = axins.plot(r,elst,linestyle=(0, (3, 1, 1, 1, 1, 1)),label='Elst')
line_exch_sub, = axins.plot(r,exch,linestyle=(0, (5, 1)),label='Exch')
line_ind_sub, = axins.plot(r,ind,linestyle=(0, (3, 1, 1, 1)),label='Ind')
line_disp_sub, = axins.plot(r,disp,linestyle=(0, (1, 1)),label='Disp')
line_sapt_sub, = axins.plot(r,sapt,'-',label='SAPT0')
axins.legend(loc=2,bbox_to_anchor=(-1.19,1.03),fontsize='large')
figure_fig = plt.gcf()
figure_fig.savefig('180-SAPT.pdf', format='pdf', dpi=1000)