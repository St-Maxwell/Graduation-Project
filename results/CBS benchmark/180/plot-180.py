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
b3lyp_sp = np.loadtxt('linear-180-B3LYPD3BJ-spline.txt')
cbs_sp = np.loadtxt('linear-180-cbs-spline.txt')
ccsdt_sp = np.loadtxt('linear-180-ccsd(t)-tz-spline.txt')
dsd_sp = np.loadtxt('linear-180-DSDPBEP86D3BJ-spline.txt')
sapt0_sp = np.loadtxt('linear-SAPT0-180-spline.txt')
b3lyp_sp_r = b3lyp_sp[:,0]
b3lyp_sp_e = b3lyp_sp[:,1]
cbs_sp_r = cbs_sp[:,0]
cbs_sp_e = cbs_sp[:,1]
ccsdt_sp_r = ccsdt_sp[:,0]
ccsdt_sp_e = ccsdt_sp[:,1]
dsd_sp_r = dsd_sp[:,0]
dsd_sp_e = dsd_sp[:,1]
sapt0_sp_r = sapt0_sp[:,0]
sapt0_sp_e = sapt0_sp[:,1]
fig, ax = plt.subplots(tight_layout=True)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
ax.set_xlabel(r'$R/\mathrm{Angstroms}$', fontsize=16)
ax.set_ylabel(r'$E/\mathrm{cm}^{-1}$', fontsize=16)
line_cbs, = ax.plot(cbs_sp_r,cbs_sp_e,'-')
line_ccsdt, = ax.plot(ccsdt_sp_r,ccsdt_sp_e,linestyle=(0, (5, 1)))
line_sapt0, = ax.plot(sapt0_sp_r,sapt0_sp_e,linestyle=(0, (3, 1, 1, 1)))
line_b3lyp, = ax.plot(b3lyp_sp_r,b3lyp_sp_e,linestyle=(0, (1, 1)))
line_dsd, = ax.plot(dsd_sp_r,dsd_sp_e,linestyle=(0, (3, 1, 1, 1, 1, 1)))
axins = inset_axes(ax,width='49%',height='60%',loc=1)
axins.set_xlim(4.55,5.25)
axins.set_ylim(-115,-10)
line_cbs_sub, = axins.plot(cbs_sp_r,cbs_sp_e,'-',label='CCSD(T)/CBS')
line_ccsdt_sub, = axins.plot(ccsdt_sp_r,ccsdt_sp_e,linestyle=(0, (5, 1)),label='CCSD(T)/avtz')
line_sapt0_sub, = axins.plot(sapt0_sp_r,sapt0_sp_e,linestyle=(0, (3, 1, 1, 1)),label='SAPT0')
line_b3lyp_sub, = axins.plot(b3lyp_sp_r,b3lyp_sp_e,linestyle=(0, (1, 1)),label='B3LYP-D3(BJ)')
line_dsd_sub, = axins.plot(dsd_sp_r,dsd_sp_e,linestyle=(0, (3, 1, 1, 1, 1, 1)),label='DSD-PBEP86-D3(BJ)')
axins.legend(loc=2,bbox_to_anchor=(-1.02,1.03),fontsize='large')
figure_fig = plt.gcf()
figure_fig.savefig('180-benchmark.pdf', format='pdf', dpi=1000)