# Put task2a.py code here
import csv
import pandas as pd
from sklearn import neighbors
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn import preprocessing
from collections import defaultdict
from sklearn.tree import DecisionTreeClassifier
import statistics
import scipy.stats as stats

# open files
fp = open("world.csv")
reader = csv.reader(fp)     
world = pd.read_csv('world.csv', index_col='Country Name')
life = pd.read_csv('life.csv', index_col='Country') 

# discard country that does not present in life.csv and store it in new_world
discard_country = [country for country in world.index if country not in life.index]
new_world = world.drop(discard_country, axis=0)
new_world = new_world.sort_values(by=['Country Code'])

discard_countries = [country for country in life.index if country not in world.index]
new_life = life.drop(discard_countries, axis=0)
new_life = new_life.sort_values(by=['Country Code'])

assert(list(new_world.index.values) == list(new_life.index.values))

org_features = new_world.columns[2:]

# merge the two into a new df
df = new_world.merge(new_life).set_index(new_world.index.values)

# calculate median of each column
stats = defaultdict(list)

for heading in org_features:
    data_list = []
    for d in new_world[heading]:
        if d != "..":
            data_list.append(float(d))
        else:
            data_list.append(0)
    med = statistics.median(data_list) 
    stats[heading].append(round(med, 3))

# implement median to missing values
new_df = defaultdict(list)
for country in new_world.index:
    new_df['Country Name'].append(country)
for head in stats.keys():
    data_list = []
    for d in new_world[head]:
        if d != "..":
            data_list.append(float(d))
        else:
            data_list.append(stats[head][0])
    new_df[head].extend(data_list)

# update median implementation in the orginial dataframe
new = pd.DataFrame(new_df)
new_world.update(new.set_index('Country Name'))

 
# get just the features
data = new_world[org_features]

# get just the class labels
classlabel = df['Life expectancy at birth (years)']

# randomly select 70% of the instances to be training and the rest to be testing
X_train, X_test, y_train, y_test = train_test_split(data, classlabel, train_size=0.70, test_size=0.30, random_state=200)

# train the decision tree model 
dt = DecisionTreeClassifier(criterion="entropy",random_state=200, max_depth=2)
dt.fit(X_train, y_train)

y_pred = dt.predict(X_test)
print(f'Accuracy of decision tree: {accuracy_score(y_test, y_pred):.3f}')

# normalise the data to have no mean and unit variance using the library functions.
scaler = preprocessing.StandardScaler(with_mean=True, with_std=True).fit(X_train)
X_train=scaler.transform(X_train)
X_test=scaler.transform(X_test)

# train the k-nearest neighbours classification model k = 3
knn = neighbors.KNeighborsClassifier(n_neighbors=3)
knn.fit(X_train, y_train)
y_pred=knn.predict(X_test)
print(f'Accuracy of k-nn (k=3): {accuracy_score(y_test, y_pred):.3f}')

# train the k-nearest neighbours classification model k = 7
knn = neighbors.KNeighborsClassifier(n_neighbors=7)
knn.fit(X_train, y_train)
y_pred=knn.predict(X_test)
print(f'Accuracy of k-nn (k=7): {accuracy_score(y_test, y_pred):.3f}')

# update the mean and variances that we used for classfication to the dictionary 
for i in range(len(scaler.mean_)):
    stats[org_features[i]].append(round(list(scaler.mean_)[i], 3))
for i in range(len(scaler.var_)):
    stats[org_features[i]].append(round(list(scaler.var_)[i], 3))

# write to an output csv file containing headings feature, median, mean, variance.
headings = ["feature", "median", "mean", "variance"]
Task2 = open("task2a.csv", "w")
writer = csv.writer(Task2)
rows = [[x, y[0],y[1],y[2]] for x, y in stats.items()]
writer.writerow(headings)
writer.writerows(rows)
writer.writerow(["\n"])
Task2.close()

fp = open("task2a.csv", "r")
# print(fp.read())

