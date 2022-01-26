import csv 
import pandas as pd   
from collections import defaultdict
import sklearn as sk
import nltk
import string
from nltk import word_tokenize
from nltk.corpus import stopwords
from nltk.stem.porter import PorterStemmer
import re

# open the files into dataframes
fp = open("abt.csv")
reader = csv.reader(fp)
df1 = pd.DataFrame(reader, columns=next(reader))
df1 = df1.set_index('idABT')


fp = open("buy.csv")
reader = csv.reader(fp)
df2 = pd.DataFrame(reader, columns=next(reader))
df2 = df2.set_index('idBuy')


# A function that preprocesses names
def preprocess(name):
    name_p = "".join([char for char in name if char not in string.punctuation])
    words = word_tokenize(name_p)
    stop_words = stopwords.words('english')
    filtered_words = [word for word in words if word not in stop_words]
    porter = PorterStemmer()
    stemmed = [porter.stem(word) for word in filtered_words]
    final = " ".join(stemmed)
    return final

pattern = r'(?<= - )[A-Z0-9]{3,20}'
pattern_abt = defaultdict(int)
f = 0
for i in df1['name']:
    if len(re.findall(pattern, i)) == 1:
        key = re.findall(pattern, i)[0]
        index = df1[df1['name'] == i].index.values[0]
        pattern_abt[key] = int(index)
        
pattern = r'(?<= - )[A-Z0-9]{3,20}'
p = r'(?<= )[A-Z0-9]+-[A-Z0-9]+'
p2 = r'(?<= )[A-Z0-9]{4,20}'
p3 = r'[A-Z0-9]{4,20}$'
p4 = r'[A-Z0-9]{4,20} -'
found = 0
pattern_buy = defaultdict(int)
for i in df2['name']:
    index = int(df2[df2['name'] == i].index.values[0])
    if len(re.findall(p, i)) == 1:
        key = re.findall(p, i)[0]
        key = key.replace("-", "")
        pattern_buy[key] = int(index)
    elif len(re.findall(pattern, i)) == 1:
        key = re.findall(pattern, i)[0]
        pattern_buy[key] = int(index)
    elif len(re.findall(p2, i)) == 1:
        key = re.findall(p2, i)[0]
        pattern_buy[key] = int(index)
    elif len(re.findall(p3, i)) == 1:
        key = re.findall(p3, i)[0]
        pattern_buy[key] = int(index)
    elif len(re.findall(p4, i)) == 1:
        key = re.findall(p4, i)[0]
        key = key.replace(" -", "")
        pattern_buy[key] = int(index)
        
def_m = []
blocks_abt = defaultdict(list)
blocks_buy = defaultdict(list)
for i in pattern_abt.keys():
    if i in pattern_buy.keys():
        blocks_abt[i].append(pattern_abt[i])
        blocks_buy[i].append(pattern_buy[i])
        def_m.append((pattern_abt[i], pattern_buy[i]))
        
# these are the ids left - create blocks for these as well 
abt_ids = [int(i) for i in df1.index.values]
pos_abt = [i for i in abt_ids if i not in [x for x,y in def_m]]
buy_ids = [int(i) for i in df2.index.values]
pos_buy = [i for i in buy_ids if i not in [y for x,y in def_m]]

for i1 in pos_abt:
    name = df1.loc[str(i1)]['name']
    n1 = preprocess(name.split()[0])
    if i1 not in blocks_abt[n1]:
        blocks_abt[n1].append(i1)
        
for i2 in pos_buy:
    name = df2.loc[str(i2)]['name']
    n2 = preprocess(name.split()[0])
    if i2 not in blocks_buy[n2]:
        blocks_buy[n2].append(i2)
        
        
def generate_record_pairs(groups_abt, groups_buy):
    record_pairs = []
    missed = 0
    for block in groups_buy.keys():
        if block not in groups_abt.keys():
            missed += 1
        group1 = groups_abt[block]
        group2 = groups_buy[block]
        for id1 in group1:
            for id2 in group2:
                record_pairs.append((id1, id2))
    
    return record_pairs

record_pairs = generate_record_pairs(blocks_abt, blocks_buy)

# write output csv files 
def generate_key_id(dic):
    output = defaultdict(str)
    for key, ids in dic.items():
        for i in ids:
            output[i] = key
    return output

key_id_abt = generate_key_id(blocks_abt)
key_id_buy = generate_key_id(blocks_buy)

# write the csv for abt blocks
headings = ["block_key", "product_id"]
abt_csv = open("abt_blocks.csv", "w")
writer = csv.writer(abt_csv)
rows = [[y, x] for x, y in key_id_abt.items()]
writer.writerow(headings)
writer.writerows(rows)


# write the csv for buy blocks
headings = ["block_key", "product_id"]
buy_csv = open("buy_blocks.csv", "w")
writer = csv.writer(buy_csv)
rows = [[y, x] for x, y in key_id_buy.items()]
writer.writerow(headings)
writer.writerows(rows)
