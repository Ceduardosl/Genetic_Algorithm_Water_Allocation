#%%
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib import rcParams

def operation(v_max, vo, va, ev, ev_alfa, ev_beta, reg):
    
    if (vo + va) < 0.001:
        ev1 = 0
    else:
        ev1 = (ev/2)*ev_alfa*np.power(vo + va, ev_beta)
    
    if ((vo + va) - ev1) > reg/2:
        v_reg1 = reg/2
    else:
        v_reg1 = max([0, ((vo + va) - ev1)])
    
    v1 = (vo + va)-(ev1 + v_reg1)

    if v1 < 0.001:
        ev2 = 0
    else:
        ev2 = (ev/2)*ev_alfa*np.power(v1, ev_beta)

    if (v1 - ev2) > reg/2:
        v_reg2 = reg/2 
    else:
        v_reg2 = max([0, (v1 - ev2)])
    
    v2 = v1 - ev2 - v_reg2
    v_reg = v_reg1 + v_reg2
    evap = ev1 + ev2
    v_vert = max([0, v2 - v_max])
    v3 = v2 - v_vert

    if (v_reg < reg):
        falha = 1
    else:
        falha = 0
    return ([v3, v_reg, evap, v_vert, falha])

#%%
input_data = pd.read_excel("{}/Data_Oros.xlsx".format(os.getcwd()), sheet_name = "input_data")
vo = 970
v_max = 1940
ev_alfa = 0.337665240404597
ev_beta = 0.842665135414265
reg_90 = 49.0338371337891
reg = np.arange(49, 67, 2)

#%%
for i in reg:
    vol_list = []
    reg_list = []
    evap_list = []
    vert_list = []
    falha_list = []
    result_df = input_data.copy()

    for j in range(0, len(input_data)):
        vol_list.append(vo)
        reg_list.append(operation(v_max, vo, input_data.iloc[j,1], input_data.iloc[j,2], ev_alfa, ev_beta, i)[1])
        evap_list.append(operation(v_max, vo, input_data.iloc[j,1], input_data.iloc[j,2], ev_alfa, ev_beta, i)[2])
        vert_list.append(operation(v_max, vo, input_data.iloc[j,1], input_data.iloc[j,2], ev_alfa, ev_beta, i)[3])
        falha_list.append(operation(v_max, vo, input_data.iloc[j,1], input_data.iloc[j,2], ev_alfa, ev_beta, i)[4])
        vo = operation(v_max, vo, input_data.iloc[j,1], input_data.iloc[j,2], ev_alfa, ev_beta, i)[0]

    result_df.insert(len(result_df.columns),"storage", vol_list)
    result_df.insert(len(result_df.columns),"vol_reg", reg_list)
    result_df.insert(len(result_df.columns),"evap", evap_list)
    result_df.insert(len(result_df.columns),"vol_vert", vert_list)
    result_df.insert(len(result_df.columns),"falha", falha_list)
    df_anual = result_df.groupby(result_df["Data"].dt.year).sum()

    garantia = 1 - sum(result_df["falha"])/len(result_df)

    result_df.to_csv("{}/Outputs/Monthly_reg_{:.0f}_garan_{:.2f}.csv".format(os.getcwd(), i, garantia))
    df_anual.to_csv("{}/Outputs/Annual_reg_{:.0f}_garan_{:.2f}.csv".format(os.getcwd(), i, garantia))
    # fig, (ax1, ax2) = plt.subplots(2,1,dpi = 600)
    # ax1.set_title("Regularização = {:.2f} hm³/mês | Garantia = {:.2%}".format(i, garantia))
    # ax1.plot(result_df["Data"], result_df["storage"], ls = "-", c = "blue", zorder = 3, label = "1-Volume / 2-Regularização")
    # ax1.set_ylabel("Armazenado (hm³)")
    # ax2.plot(result_df["Data"], result_df["vol_reg"], ls = "-", c = "blue", zorder = 3, label = "Vazão Regularizada")
    # ax2.set_ylabel("Regularizado (hm³/mês)")
    # fig.savefig("{}/Figuras/Figura_reg_{:.2f}.png".format(os.getcwd(), i), dpi = 600, bbox_inches = "tight", facecolor = "w")
#%%












# %%
reg = np.arange(0, 1005, 5)
garant_list = []
for j in reg:
    falha_list = []
    for i in range(0, len(input_data)):
        vol_list.append(vo)
        falha_list.append(operation(v_max, vo, input_data.iloc[i,1], input_data.iloc[i,2], ev_alfa, ev_beta, j)[4])
        vo = operation(v_max, vo, input_data.iloc[i,1], input_data.iloc[i,2], ev_alfa, ev_beta, j)[0]

    garant_list.append(1 - sum(falha_list)/(len(falha_list)))

result_df = pd.DataFrame({"Regularizacao": reg, "Garantia": garant_list})
# %%
fig, ax = plt.subplots(dpi = 600)
ax.plot(result_df["Regularizacao"], result_df["Garantia"], c = "black", zorder = 1)
ax.plot(result_df["Regularizacao"].loc[result_df["Garantia"] > 0.8], result_df["Garantia"].loc[result_df["Garantia"] > 0.8], c = "blue", zorder = 3)
ax.scatter(reg_90, 0.9, c = "red", zorder = 3)
ax.annotate("({:.2f}, {:.2%})".format(reg_90, 0.9), (reg_90+25, 0.9))
ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
ax.set_xlabel("hm³/mês")
ax.set_ylabel("Garantia")
ax.set_title("Curva de Regularização", loc = "center")
fig.savefig("{}/Curva_Regularizacao.png".format(os.getcwd()), dpi = 600, bbox_inches = "tight", facecolor = "w")

# %%
