import pandas as pd
from tabulate import tabulate

df = pd.read_excel('/home/mikh_kate/kalenskii/CASA/lines_new.ods')


mol_column = df['Молекула']
mol_count = mol_column.value_counts()


unique_molecules = mol_count.index
counts = mol_count.values


table_data = []
for idx, mol in enumerate(unique_molecules, start=1):
    table_data.append([idx, mol, counts[idx - 1]])


headers = ["Номер молекулы", "Название молекулы", "Количество встречаний"]
latex_table = tabulate(table_data, headers=headers, tablefmt="latex")

print(table_data)
df = pd.DataFrame(table_data)

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

# Print the entire DataFrame
#print(df)



print(latex_table)
