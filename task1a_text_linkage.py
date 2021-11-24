import csv
import textdistance as td   
import pandas as pd   
from collections import defaultdict
import re
import string
from nltk import word_tokenize
from nltk.corpus import stopwords
from nltk.stem.porter import PorterStemmer

# open the files into dataframes
fp = open("abt_small.csv")
reader = csv.reader(fp)
df1 = pd.DataFrame(reader, columns=next(reader))
df1 = df1.set_index('idABT')
# display(df1)

fp = open("buy_small.csv")
reader = csv.reader(fp)
df2 = pd.DataFrame(reader, columns=next(reader))
df2 = df2.set_index('idBuy')
# display(df2)

# A function that preprocesses names
def prep(name):
    # punctuation removal
    name_p = "".join([char for char in name.lower() if char not in string.punctuation])
    # tokenization
    words = word_tokenize(name_p)
    # stop word removal 
    stop_words = stopwords.words('english')
    filtered_words = [word for word in words if word not in stop_words]
    porter = PorterStemmer()
    stemmed = [porter.stem(word) for word in filtered_words]
    final = " ".join(stemmed)
    return final

# preprocess names       
pre_name1 = []
for name in df1['name']:
    pre_name1.append(prep(name))
df1['prename'] = pre_name1

pre_name2 = []
for name in df2['name']:
    pre_name2.append(prep(name))
df2['prename'] = pre_name2

# sort based on their name
df1 = df1.sort_values(by=['name'])
df2 = df2.sort_values(by=['name'])

# generate all possible pairs
all_pairs = []
for i1 in df1.index.values:
    for i2 in df2.index.values:
        all_pairs.append((int(i1), int(i2)))
# find the model code in Abt names
pattern = r'(?<= - )[A-Z0-9]{3,20}'
pattern_abt = defaultdict(int)
f = 0
for i in df1['name']:
    if len(re.findall(pattern, i)) == 1:
        key = re.findall(pattern, i)[0]
        index = df1[df1['name'] == i].index.values[0]
        pattern_abt[key] = int(index)
    else:
        print(re.findall(pattern, i))

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

# filter out definite matches first     
def_m = []
for i in pattern_abt.keys():
    if i in pattern_buy.keys():
        def_m.append((pattern_abt[i], pattern_buy[i]))

abt_ids = [int(i) for i in df1.index.values]
pos_abt = [i for i in abt_ids if i not in [x for x,y in def_m]]
buy_ids = [int(i) for i in df2.index.values]
pos_buy = [i for i in buy_ids if i not in [y for x,y in def_m]]

# test for possible matches 
matches = []
fil = defaultdict(list)
for i1 in pos_abt:
    aname = df1.loc[str(i1)]['prename']
    for i2 in pos_buy:
        bname = df2.loc[str(i2)]['prename']
        if aname.split()[0] in bname.split():
            if td.cosine(aname,bname) > 0.5:
                fil[i1].append((td.cosine(aname,bname), i2))
for i1, i2s in fil.items():
    fil[i1].sort(reverse=True)
    fil[i1] = [y for x,y in fil[i1][:3]]   

pattern = r'[A-Z0-9]+[0-9]+'
final_filter = defaultdict(int)
for i1, i2s in fil.items():
    s1 = df1.loc[str(i1)]['name']
    s1 = s1.translate(str.maketrans('', '', string.punctuation))
    for i2 in i2s:
        s2 = df2.loc[str(i2)]['name']
        s2 = s2.translate(str.maketrans('', '', string.punctuation))
        score = len(set(re.findall(pattern, s1)).intersection(set(re.findall(pattern, s2))))
        if score >= 1:
            final_filter[i1] = i2
            
record_pairs = def_m + list(final_filter.items())  
# write the output file
headings = ["idAbt", "idBuy"]
task1 = open("task1a.csv", "w")
writer = csv.writer(task1)
rows = [[x, y] for x, y in record_pairs]
writer.writerow(headings)
writer.writerows(rows)

