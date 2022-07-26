import pandas as pd

pd.set_option('display.max_columns', None)

planet = 'Sarr'
df = pd.read_csv('outputs/2022-07-24-12-32/csvs/' + planet + '/ref' + planet + '.csv')

print(df)
print(df.loc[df['similarityIndex_' + planet].idxmax()])
