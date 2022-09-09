#!/data/SW/anaconda3/envs/myenv/bin/python

from pylab import *
from matplotlib.patches import Patch
import seaborn as sns; sns.set(color_codes=True)
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable
from statannotations.Annotator import Annotator

df=pd.read_csv("log_odds_pair_pos_new.tsv", comment="#", sep="\t")
sns.set(style="ticks",rc={'axes.facecolor':'white','figure.figsize': (9, 6)})
fig, (ax, ax1) = plt.subplots(nrows=2, sharex=True)
plt.subplots_adjust(wspace=0, hspace=0)
subplot(2,1,1)
ax=sns.kdeplot(data=df, x='Gio')
y= ax.lines[-1].get_ydata()
x=ax.lines[-1].get_xdata()
med=df['Gio'].median()
plt.text(0.46, 0.5, 'median', horizontalalignment='center',
     verticalalignment='center', transform=ax.transAxes,rotation=90,fontsize=15)
ax.axes.yaxis.set_ticklabels([])
ax.tick_params(left=False)
ax.axes.yaxis.set_visible(False)
ax.set_xticks(np.arange(-7,7,1))
plt.xticks(rotation=90)
plt.fill_between(x,0,y, where = x<=-1, color='#a50026')
plt.fill_between(x,0,y, where = x>=1, color='#006837')
plt.fill_between(x,0,y, where = (x>=0) & (x<=1), color='#b7e075')
plt.fill_between(x,y, where = (x>=-1) & (x<=0), color='#fdbf6f')
plt.text(0.1, 0.5, 'Gio', horizontalalignment='center',
     verticalalignment='center', transform=ax.transAxes,rotation=90,fontsize=25)
plt.gca().invert_xaxis()
ax.vlines(med, 0,y.max(), linestyle='--', alpha=1,color='black')

subplot(2,1,2)
ax1=sns.kdeplot(data=df, x='Gs')
ax1.axes.yaxis.set_ticklabels([])
ax1.axes.yaxis.set_visible(False)
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(left=False)
ax1.set_xticks(np.arange(-7,7,0.5))
plt.xticks(rotation=90)
med=df['Gs'].median()
y1= ax1.lines[-1].get_ydata()
x1=ax1.lines[-1].get_xdata()
plt.fill_between(x1,0,y1, where = x1<-1, color='#a50026')
plt.fill_between(x1,0,y1, where = x1>1, color='#006837')
plt.fill_between(x1,0,y1, where = (x1>=0) & (x1<1), color='#b7e075')
plt.fill_between(x1,0,y1, where = (x1>-1) & (x1<=0), color='#fdbf6f')
plt.gca().invert_yaxis()
plt.xlabel("logodds",fontsize=25,rotation=180)
plt.text(0.1, 0.5, 'Gs', horizontalalignment='center',
     verticalalignment='center', transform=ax1.transAxes,rotation=90,fontsize=25)
plt.text(0.39, 0.5, 'median', horizontalalignment='center',
     verticalalignment='center', transform=ax1.transAxes,rotation=90,fontsize=15)
ax1.vlines(med, 0,y1.max(), linestyle='--', alpha=1,color='black')


pairs=[('Gs', 'Gio')]
stat,pval=stats.ranksums(df["Gs"], df["Gio"],nan_policy='omit')
plt.text(6.42, 0.25, 'Ranksum  pval='+"{:.2e}".format(pval), fontsize = 22,rotation=90)
fig.savefig('log_odds_pair.svg',dpi=600,bbox_inches='tight',format='svg')
