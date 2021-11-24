## final

# for task 2b
# 10.3 version
import csv
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import neighbors
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn import preprocessing
from collections import defaultdict
import statistics
import sklearn.preprocessing as preproc
from sklearn.cluster import KMeans
import numpy as np
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.preprocessing import StandardScaler

fp = open("world.csv")
reader = csv.reader(fp)  
world = pd.read_csv('world.csv', index_col='Country Name')
life = pd.read_csv('life.csv', index_col='Country') 


# discard country that does not present in life.csv and store it in new_world
discard_country = [country for country in world.index if country not in life.index]
new_world = world.drop(discard_country, axis=0)
discard_countries = [country for country in life.index if country not in world.index]
new_life = life.drop(discard_countries, axis=0)
new_life = new_life.sort_values(by=['Country Code'])
org_features = new_world.columns[2:]


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

# update median implementation in the orginial dataframe, and sort based on Country Code
new = pd.DataFrame(new_df)
new_world.update(new.set_index('Country Name'))
new_world = new_world.sort_values(by=['Country Code'])

# merge the two into a new df
df = new_world.merge(new_life).set_index(new_world.index.values) 

# first, split the training and test data and standarize them 

# X is a df containing original 20 features
X = new_world.drop(['Time', 'Country Code'], axis=1).astype(float)

# the class labels
target = df['Life expectancy at birth (years)']

X_train, X_test, y_train, y_test = train_test_split(X, target, train_size=0.70, test_size=0.30, random_state=200)

# generate interaction term pairs for both training and testing set
X_train_int = preproc.PolynomialFeatures(interaction_only = True, include_bias=False).fit_transform(X_train)
X_test_int = preproc.PolynomialFeatures(interaction_only = True, include_bias=False).fit_transform(X_test)

# create separate dfs for training and testing datasets
train_df = pd.DataFrame(X_train_int[:, :])
test_df = pd.DataFrame(X_test_int[:, :])

# a dictionary of country names
country_names = defaultdict(list)
country_names['Country Name'].extend(new_world.index)

# name the columns based on feature names + create dfs 
train_df.columns = [list(org_features) + ['Feature ' + str(x) for x in range(1, 191)]]
test_df.columns = [list(org_features) + ['Feature ' + str(x) for x in range(1, 191)]]

# Evidence  of  correct  implementation  of  interaction  term  pair
print(train_df[['Feature ' + str(x) for x in range(1, 191)]].sample(10))
print(test_df[['Feature ' + str(x) for x in range(1, 191)]].sample(10))
 
x = np.array(train_df.values)

# use elbow method to select the optimal number of clusters
wcss = []
for i in range(1,11): 
    kmeans = KMeans(n_clusters=i, init ='k-means++', max_iter=300, n_init=10,random_state=200)
    kmeans.fit(x)
    wcss.append(kmeans.inertia_)

plt.plot(range(1,11),wcss)
plt.title('The Elbow Method Graph')
plt.xlabel('Number of clusters')
plt.ylabel('WCSS')
plt.savefig("task2bgraph1.png")

# the optimal cluster number is 3

# we need to do cluster labels separately for training set and testing set 
x_t = np.array(train_df.values)
scaler = StandardScaler()
x_t = scaler.fit_transform(x_t)

# create one feature using clustering 
kmeans = KMeans(n_clusters=3, init ='k-means++', max_iter=300, n_init=10,random_state=200 )
y_kmeans = kmeans.fit_predict(x_t)

# we add cluster labels to a separate dictionary
clusters = defaultdict(list)

for i in range(len(x_t)):
    predict_me = np.array(x_t[i].astype(float))
    predict_me = predict_me.reshape(-1, len(predict_me))
    prediction = kmeans.predict(predict_me)
    clusters['c'].append(prediction[0])

train_df.insert(0, "cluster label", list(clusters['c']), True)

# Evidence  of  correct  implementation  of  feature  engineered  from  clustering
print(train_df[["cluster label"]].sample(10))

# for the test set
x_test = np.array(test_df.values)
scaler = StandardScaler()
x_test = scaler.fit_transform(x_test)

# create one feature using clustering 
kmeans = KMeans(n_clusters=3, init ='k-means++', max_iter=300, n_init=10,random_state=200 )
y_kmeans = kmeans.fit_predict(x_test)

# we add this feature
test_clusters = defaultdict(list)

for i in range(len(x_test)):
    predict_me = np.array(x_test[i].astype(float))
    predict_me = predict_me.reshape(-1, len(predict_me))
    prediction = kmeans.predict(predict_me)
    test_clusters['c'].append(prediction[0])

test_df.insert(0, "cluster label", list(test_clusters['c']), True)
print(test_df[["cluster label"]].sample(10))

# feature selection using random forests
X_train_rf = train_df.values
y_train_rf = y_train.values

# feature selection using random forest
sel = SelectFromModel(RandomForestClassifier(n_estimators = 100))
sel.fit(X_train_rf, y_train_rf)
selected_feat= train_df.columns[(sel.get_support())]

# evidence of selected 4 features
print(selected_feat[:4])

x_best = train_df[selected_feat[:4]]


# plot feature importances for four selected features
features = ['F1', 'F2', 'F3', 'F4']
importances = sel.estimator_.feature_importances_
indices = np.argsort(importances)[::-1][:4]

plt.title('Feature Importances of Selected four most important features')
plt.barh(range(len(indices)), importances[indices], color='b', align='center')
plt.xlabel('Relative Importance')
plt.ylabel("Feature indexes")
plt.savefig("task2bgraph2.png")

importances = sel.estimator_.feature_importances_
indices = np.argsort(importances)[::-1]

# X is the train data used to fit the model 
plt.figure()
plt.title("Feature importances of all features")
plt.bar(range(X_train_rf.shape[1]), importances[indices],color="r", align="center")
plt.xticks(range(X_train_rf.shape[1]), indices)
plt.xlim([-1, X_train_rf.shape[1]])
plt.xlabel("Feature indexes")
plt.ylabel("Relative importance")
plt.savefig("task2bgraph3.png")


# visualize the pair plots of four selected features
x = np.array(x_best)

rd_df = pd.DataFrame(x, columns=['F1', 'F2', 'F3', 'F4'])

labels = []
for i in list(y_train.values):
    if i == 'Low':
        labels.append(0)
    if i == 'Medium':
        labels.append(1)
    if i == 'High':
        labels.append(2)

grr = pd.plotting.scatter_matrix(rd_df, c=labels, figsize=(10,10), marker='o', hist_kwds={'bins':20}, s=10, alpha=1.5, cmap='brg')

handles = [plt.plot([],[],color=plt.cm.brg(i/2.), ls="", marker=".", \
                    markersize=np.sqrt(10))[0] for i in range(3)]
classes=["Low", "Medium", "High"]
plt.legend(handles, classes, loc=(1.02,0))
plt.suptitle( "Pair Plots for 4 Features selected using Random Forest")
plt.savefig("task2bgraph6.png")

# standarize before applying k-nearest neighbours
scaler = preprocessing.StandardScaler(with_mean=True, with_std=True).fit(x_best)
X_tr=scaler.transform(x_best)
X_test=scaler.transform(test_df[selected_feat[:4]])

# train the k-nearest neighbours classification model k = 3 using 4 selected features by RFs
knn = neighbors.KNeighborsClassifier(n_neighbors=3)
knn.fit(X_tr, y_train_rf)
y_pred=knn.predict(X_test)

print(f'Accuracy of feature engineering: {accuracy_score(y_test, y_pred):.3f}')

# then, we do PCA analysis of original 20 features

# first, split the training and test data and standarize them 
X = new_world.drop(['Time', 'Country Code'], axis=1).astype(float)

target = df['Life expectancy at birth (years)']

X_train, X_test, y_train, y_test = train_test_split(X, target, train_size=0.70, test_size=0.30, random_state=200)

## Performing standardization before applying PCA
scaler = preprocessing.StandardScaler(with_mean=True, with_std=True).fit(X_train)
X_train=scaler.transform(X_train)
X_test=scaler.transform(X_test)

pca = PCA(n_components = 4)
fit = pca.fit(X_train)
pca.fit(X_test)

# print the explained variance ratios
print(pca.explained_variance_ratio_)

PCA(copy=True, iterated_power='auto', n_components=23, random_state=None,svd_solver='auto', tol=0.0, whiten=False)

X_train_pca = pca.transform(X_train)
X_test_pca = pca.transform(X_test)

# a dataframe containing the 4 PCs, evidence for the PCs generated
pcDF = pd.DataFrame(data=X_train_pca, columns = ['pc1', 'pc2', 'pc3', 'pc4'])
print(pcDF.sample(10))

# train the k-nearest neighbours classification model k = 3 using 4 pca features 
knn = neighbors.KNeighborsClassifier(n_neighbors=3)
knn.fit(X_train_pca, y_train)
y_pred=knn.predict(X_test_pca)

print(f'Accuracy of PCA: {accuracy_score(y_test, y_pred):.3f}')

# we visualize the interactions between 4PCs using the pandas library
x = np.array(pcDF.values)

pcdf = pd.DataFrame(x, columns=['PC1', 'PC2', 'PC3', 'PC4'])

labels = []
for i in list(y_train.values):
    if i == 'Low':
        labels.append(0)
    if i == 'Medium':
        labels.append(1)
    if i == 'High':
        labels.append(2)


grr = pd.plotting.scatter_matrix(pcdf, c=labels, figsize=(10,10), marker='o', hist_kwds={'bins':20}, s=10, alpha=1.5, cmap='brg')

handles = [plt.plot([],[],color=plt.cm.brg(i/2.), ls="", marker=".", \
                    markersize=np.sqrt(10))[0] for i in range(3)]
classes=["Low", "Medium", "High"]
plt.legend(handles, classes, loc=(1.02,0))
plt.suptitle( "Pair Plots for 4PCs selected using PCA")
plt.savefig("task2bgraph4.png")

# Next, take the first four features from the original dataset as a sample of the original 20 features.  
# Perform 3-NN classification

x_data = new_world[org_features[:4]]
y = target

# randomly select 70% of the instances to be training and the rest to be testing
X_train, X_test, y_train, y_test = train_test_split(x_data, y, train_size=0.70, test_size=0.30, random_state=200)

# normalise the data to have no mean and unit variance using the library functions.
scaler = preprocessing.StandardScaler(with_mean=True, with_std=True).fit(X_train)
X_train=scaler.transform(X_train)
X_test=scaler.transform(X_test)

# train the k-nearest neighbours classification model k = 3 using the first four features
knn = neighbors.KNeighborsClassifier(n_neighbors=3)
knn.fit(X_train, y_train)
y_pred=knn.predict(X_test)

print(f'Accuracy of first four features: {accuracy_score(y_test, y_pred):.3f}')

# visualize the pair plots of first four features
d = pd.DataFrame(X_train)
x = np.array(d[list(d.columns)[:4]].values)

# we name the frist 4 features as 'A', 'B', 'C', 'D'
feature_df = pd.DataFrame(x, columns=['A', 'B', 'C', 'D'])

# we colour code our points on the graph
labels = []
for i in list(y_train.values):
    if i == 'Low':
        labels.append(0)
    if i == 'Medium':
        labels.append(1)
    if i == 'High':
        labels.append(2)

# we plot the pair plots for the four features
grr = pd.plotting.scatter_matrix(feature_df, c=labels, figsize=(10,10), marker='o', hist_kwds={'bins':20}, s=10, alpha=1.5, cmap='brg')

handles = [plt.plot([],[],color=plt.cm.brg(i/2.), ls="", marker=".", \
                    markersize=np.sqrt(10))[0] for i in range(3)]
classes=["Low", "Medium", "High"]
plt.legend(handles, classes, loc=(1.02,0))
plt.suptitle( "Pair Plots for 4 original features")
plt.savefig("task2bgraph5.png")
