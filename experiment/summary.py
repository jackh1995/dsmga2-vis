from os import listdir
from os.path import isfile, join
import pandas as pd

csv_list = [f for f in listdir('./results') if isfile(join('./results', f))]

df_list = [pd.read_csv(f'./results/{csv}', names=['pop', 'gen', 'nfe']) for csv in csv_list]

for csv, df in zip(csv_list, df_list):
    print(f'[{csv}]')
    print(df.describe())
    print()