import pandas as pd

#Sorting data

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

df1 = pd.read_csv('get_html.csv')
df = df1.sort_values(['Frequency (MHz)'])
print(df[['Frequency (MHz)', 'Quantum Numbers']])
#print(df)

df.to_csv('get_html_sorted.csv')
