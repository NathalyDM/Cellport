from variables import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

frames_golgi=[df_CD59FL_golgi,df_GPIFRA_golgi,df_GPICD59_golgi,df_PCX_golgi,df_Ecad_golgi]
results_golgi = pd.concat(frames_golgi)

frames_ER=[df_CD59FL_ER,df_GPIFRA_ER,df_GPICD59_ER,df_PCX_ER,df_Ecad_ER]
results_ER = pd.concat(frames_ER)

results_ER_df2 = results_ER.rename({'D(min-1)': 'D', 'V(min-1)': 'V',
                                    'kL(min-1)':'kL','kR(min-1)':'kR'}, axis='columns')


results_golgi_df2 = results_golgi.rename({'D(min-1)': 'D', 'V(min-1)': 'V',
                                    'kL(min-1)':'kL','kR(min-1)':'kR'}, axis='columns')

plt.figure()
# make boxplot with Seaborn
bplot = sns.barplot(y='kR', x='System',
                    data=results_golgi_df2,
                    width=0.5,
                    palette="pastel")

# add stripplot to boxplot with Seaborn
bplot = sns.stripplot(y='kR', x='System',
                      data=results_golgi_df2,
                      jitter=True,
                      marker='o',
                      alpha=0.5,
                      color='black')
bplot
plt.show()