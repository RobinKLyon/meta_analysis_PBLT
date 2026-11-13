import pandas as pd
results_stats_all_patients = pd.read_csv(r'C:\Users\Robin\Documents\THESE_CODE\PROJET_META_ANALYSE\Data_manadgement\outputs\analyse_pblt\resultats_statistiques_tous_patients.csv')
print(results_stats_all_patients[results_stats_all_patients["bands"]=="HFA"].sum())