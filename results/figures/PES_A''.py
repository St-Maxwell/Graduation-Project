import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in' 
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams["xtick.top"] = True
mpl.rcParams["ytick.right"] = True
mpl.rcParams['font.family'] = 'Arial'
import numpy as np
fig, ax = plt.subplots(tight_layout=True)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
ax.set_xlabel(r'$\theta/\mathrm{degrees}$', fontsize=14)
ax.set_ylabel(r'$R/\mathrm{Angstroms}$', fontsize=14)
z = np.loadtxt('A1.txt')
x = np.linspace(0,180,37)
y = np.linspace(2.7,7,44)
X,Y = np.meshgrid(x,y)
lines = [-150,-110,-80,-40,-15,-5,100,500,2000,7000,20000]
l_color = ['blue','blue','blue','blue','blue','blue','red','red','red','red','red']
l_style = ['dashed','dashed','dashed','dashed','dashed','dashed','solid','solid','solid','solid','solid']
lbs = [-150,-110,-80,-40,-15,-5,100,500,2000,7000,20000]
lb_pos = [(100,3.75),(100,4.2),(90,4.3),(90,4.9),(90,6), \
(90,6.8),(140,3.9),(90,2.9),(40,3.4),(30,3.2), \
(20,3),(140,3.4),(150,3.2),(160,3),(10,4.5),(170,4.7),(40,3.8)]
C = plt.contour(X, Y, z, levels=lines, colors=l_color, linestyles=l_style, linewidth=.5)
plt.clabel(C, lbs, inline=True, fontsize=12, inline_spacing=3.,manual=lb_pos,fmt='%d')
figure_fig = plt.gcf()
figure_fig.savefig('figure.pdf', format='pdf', dpi=1000)

# mpl.rcParams["contour.negative_linestyle"] = 'dashed'
# mpl.rcParams['text.usetex'] = True
# ax.set_xlabel(r'$\displaystyle \theta/\mathrm{deg}$', fontsize=14)
# ax.set_ylabel(r'$\displaystyle R/\mathrm{Angstrom}$', fontsize=14)
# C = plt.contourf(X, Y, z, levels=lines, alpha=.75, cmap=plt.cm.hot)


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
ax.set_xlabel(r'$\theta/\mathrm{degrees}$', fontsize=16)
ax.set_ylabel(r'$R/\mathrm{Angstroms}$', fontsize=16)
e_pro = np.loadtxt('A2-energy_profile.txt')
ept = e_pro[:,0]
epr = e_pro[:,1]
ax.plot(ept,epr,color='royalblue',linestyle='-.',linewidth=1.2)
z = np.loadtxt('A2.txt')
x = np.linspace(0,180,37)
y = np.linspace(2.7,7,44)
X,Y = np.meshgrid(x,y)
lines = [-150,-110,-80,-40,-15,-5,100,500,2000,7000,20000]
linesf = [z.min(),-150,-110,-80,-40,-15,-5,100,500,2000,7000,20000,z.max()]
lbs = [-150,-110,-80,-40,-15,-5,100,500,2000,7000,20000]
lb_pos = [(85,3.75),(85,4),(90,4.3),(90,4.9),(90,6), \
(90,6.8),(90,2.9),(140,4),(40,3.8),(140,3.5),(40,3.3), \
(150,3.3),(30,3.2),(20,3),(160,3.1),(10,4.5),(170,4.7)]
C = plt.contour(X, Y, z, levels=lines, linewidths=0.9, colors='k')
plt.clabel(C, lbs, inline=True, fontsize=12, inline_spacing=3.,manual=lb_pos,fmt='%d')
plt.contourf(X, Y, z, levels=linesf,norm=colors.PowerNorm(gamma=0.16),cmap='coolwarm')
fig.set_size_inches(7, 5.5)
figure_fig = plt.gcf()
figure_fig.savefig('A2.pdf', format='pdf', dpi=1000)

#plt.clabel(C, lbs, inline=True, fontsize=12, inline_spacing=3.,manual=lb_pos,fmt='%d')

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
fig, ax = plt.subplots(tight_layout=True,figsize=(5,3))
xtcs = [0,20,40,60,80,100,120,140,160,180]
plt.xticks(ticks=xtcs,fontsize=14)
ytcs = [-170,-150,-130,-110,-90]
plt.yticks(ticks=ytcs,fontsize=14)
ax.set_xlabel(r'$\theta/\mathrm{degrees}$', fontsize=16)
ax.set_ylabel(r'$E/\mathrm{cm}^{-1}$', fontsize=16)
e_pro = np.loadtxt('A2-energy_profile_curve.txt')
ept = e_pro[:,0]
epr = e_pro[:,1]
ax.plot(ept,epr,color='royalblue',linewidth=1.2)
ax.set_xlim(0, 180)
ax.set_ylim(-175, -85)
ax.grid(linestyle=':')
ax.set_facecolor('whitesmoke')
fig.set_size_inches(7, 4)
figure_fig = plt.gcf()
figure_fig.savefig('A2-energy_profile.pdf', format='pdf', dpi=1000)