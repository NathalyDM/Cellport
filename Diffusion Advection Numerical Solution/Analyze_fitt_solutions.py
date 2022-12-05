from variables import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from bioinfokit.analys import stat
warnings.filterwarnings('ignore')



frames_golgi=[df_CD59FL_golgi,df_GPIFRA_golgi,df_GPICD59_golgi,df_PCX_golgi,df_Ecad_golgi]
results_golgi = pd.concat(frames_golgi)

frames_ER=[df_CD59FL_ER,df_GPIFRA_ER,df_GPICD59_ER,df_PCX_ER,df_Ecad_ER]
results_ER = pd.concat(frames_ER)

results_ER_df2 = results_ER.rename({'D(min-1)': 'D', 'V(min-1)': 'V',
                                    'kL(min-1)':'kL','kR(min-1)':'kR'}, axis='columns')


results_golgi_df2 = results_golgi.rename({'D(min-1)': 'D', 'V(min-1)': 'V',
                                    'kL(min-1)':'kL','kR(min-1)':'kR'}, axis='columns')


# Post hoc comparison
print('Diffusion')
res = stat()
res.tukey_hsd(df=results_ER_df2, res_var='D', xfac_var='System', anova_model='D ~ C(System)')
print(res.tukey_summary)

# Post hoc comparison
print('Velocity')
res = stat()
res.tukey_hsd(df=results_ER_df2, res_var='V', xfac_var='System', anova_model='V ~ C(System)')
print(res.tukey_summary)

# Post hoc comparison
print('Retrograde Transport')
res = stat()
res.tukey_hsd(df=results_ER_df2, res_var='kL', xfac_var='System', anova_model='kL ~ C(System)')
print(res.tukey_summary)

# Post hoc comparison
print('Anterograde Transport')
res = stat()
res.tukey_hsd(df=results_ER_df2, res_var='kR', xfac_var='System', anova_model='kR ~ C(System)')
print(res.tukey_summary)


print('\n','GOLGI')
# Post hoc comparison
print('Diffusion')
res = stat()
res.tukey_hsd(df=results_golgi_df2, res_var='D', xfac_var='System', anova_model='D ~ C(System)')
print(res.tukey_summary)

# Post hoc comparison
print('Velocity')
res = stat()
res.tukey_hsd(df=results_golgi_df2, res_var='V', xfac_var='System', anova_model='V ~ C(System)')
print(res.tukey_summary)

# Post hoc comparison
print('Retrograde Transport')
res = stat()
res.tukey_hsd(df=results_golgi_df2, res_var='kL', xfac_var='System', anova_model='kL ~ C(System)')
print(res.tukey_summary)

# Post hoc comparison
print('Anterograde Transport')
res = stat()
res.tukey_hsd(df=results_golgi_df2, res_var='kR', xfac_var='System', anova_model='kR ~ C(System)')
print(res.tukey_summary)

"""
sns.set_style("whitegrid")
# Difussion
plt.figure()
sns.FacetGrid(results_ER, hue = 'System', aspect=2).map(sns.distplot, 'D(min-1)')
plt.legend(loc='upper left')
plt.show()

plt.figure()
sns.FacetGrid(results_golgi, hue = 'System', aspect=2).map(sns.distplot, 'D(min-1)')
plt.legend(loc='upper left')
plt.show()

# Velocity
plt.figure()
g=sns.FacetGrid(results_ER, hue = 'System', aspect=2).map(sns.distplot, 'V(min-1)')
g.set(xlim=(0, 0.05))
plt.legend(loc='upper left')
plt.show()

plt.figure()
g=sns.FacetGrid(results_golgi, hue = 'System', aspect=2).map(sns.distplot, 'V(min-1)')
g.set(xlim=(0, 0.05))
plt.legend(loc='upper left')
plt.show()

# kL(min-1)
plt.figure()
sns.FacetGrid(results_ER, hue = 'System', aspect=2).map(sns.distplot, 'kL(min-1)')
plt.legend(loc='upper left')
plt.show()

plt.figure()
sns.FacetGrid(results_golgi, hue = 'System', aspect=2).map(sns.distplot, 'kL(min-1)')
plt.legend(loc='upper left')
plt.show()

# kR(min-1)
plt.figure()
sns.FacetGrid(results_ER, hue = 'System', aspect=2).map(sns.distplot, 'kR(min-1)')
plt.legend(loc='upper left')
plt.show()

plt.figure()
sns.FacetGrid(results_golgi, hue = 'System', aspect=2).map(sns.distplot, 'kR(min-1)')
plt.legend(loc='upper left')
plt.show()


"""



"""

results_golgi = results_ER
sns.set_style("whitegrid")
colors_list = ['blue', 'black',  'red','green']
plt.figure()
ax=sns.violinplot(x="System", y="D(min-1)", data=results_golgi,palette=colors_list, alpha=0.1)
plt.setp(ax.collections, alpha=.3)
#plt.set(ylim=(0, 0.4))
plt.show()

plt.figure()
ax=sns.boxplot(x="System", y="V(min-1)", data=results_golgi,palette=colors_list)
plt.setp(ax.collections, alpha=.3)
plt.show()

plt.figure()
ax = sns.violinplot(x="System", y="kL(min-1)", data=results_golgi,palette=colors_list, alpha=0.1)
plt.setp(ax.collections, alpha=.3)
#ax.set(ylim=(0, 0.1))
plt.show()

plt.figure()
ax = sns.violinplot(x="System", y="kR(min-1)", data=results_golgi,palette=colors_list, alpha=0.1)
plt.setp(ax.collections, alpha=.3)
#ax.set(ylim=(0, 0.1))
plt.show()

plt.figure()
ax = sns.barplot(x="System", y="Dt(min)", data=results_golgi,palette=colors_list, alpha=0.1)
plt.show()

"""

