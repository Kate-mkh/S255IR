import pandas as pd
from tabulate import tabulate


df = pd.read_excel('/home/mikh_kate/kalenskii/CASA/lines.ods', usecols = [4, 3, 12, 11])
df2 = df[["Молекула", "Каталог", "Переход", "Комментарии"]]
df_sorted = df2.sort_values(by = ["Молекула", "Каталог"])


df_sorted = df_sorted.drop_duplicates()


df_sorted = df_sorted.fillna('')2

#pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

#print(df_sorted)



table_data = []
for _, row in df_sorted.iterrows():
    table_data.append([row["Молекула"], row["Каталог"], row["Переход"], row["Комментарии"]])


headers = ["Молекула", "Частота", "Переход", "Комментарии"]
latex_table = tabulate(table_data, headers=headers, tablefmt="latex")
print(latex_table)
