#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 06:06:08 2020

@author: lynnezhao
"""
### !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  5 06:13:55 2020

@author: lynnezhao
"""

import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin
import re
import matplotlib.pyplot as plt
from collections import defaultdict
import csv
import json

page_limit = 200

#Specify the initial page to crawl
base_url = 'http://comp20008-jh.eng.unimelb.edu.au:9889/main/'
seed_item = 'Hodg001.html'

# create a seed url and parse the page
seed_url = base_url + seed_item
page = requests.get(seed_url)
soup = BeautifulSoup(page.text, 'html.parser')

# a dictionary of visited urls
visited = {}; 
visited[seed_url] = True
pages_visited = 1

# Remove index.html
links = soup.findAll('a')
seed_link = soup.findAll('a', href=re.compile("^index.html"))
to_visit_relative = [l for l in links if l not in seed_link]


# Resolve to absolute urls
to_visit = []
for link in to_visit_relative:
    to_visit.append(urljoin(seed_url, link['href']))

   
# Find all outbound links on succsesor pages and explore each one 
while (to_visit):
    
    # Impose a limit to avoid breaking the site 
    if pages_visited == page_limit :
        break
        
    # consume the list of urls
    link = to_visit.pop(0)

    # need to concat with base_url
    page = requests.get(link)
    
    # scarping the page 
    soup = BeautifulSoup(page.text, 'html.parser')
    
    # add the item to visited list and remove from to_visit
    visited[link] = True
    to_visit
    
    # find links on existing website
    new_links = soup.findAll('a')
    for new_link in new_links:
        new_item = new_link['href']
        new_url = urljoin(link, new_item)
        if new_url not in visited and new_url not in to_visit:
            to_visit.append(new_url)
        
    pages_visited = pages_visited + 1


# find and extract the headline of the article 
for url in visited.keys():
    page = requests.get(url)
    soup = BeautifulSoup(page.text, 'html.parser')
    headline = soup.find(id='headline').get_text()
    visited[url] = headline.strip('"')


# write the url and headline into an output csv file 
headings = ['url', 'headline']
rows = [[x, y] for x, y in visited.items()]
task1 = open("Task1.csv", 'w')
writer = csv.writer(task1)
writer.writerow(headings)
writer.writerows(rows)
task1.close()

# f = open("Task1.csv", "r")
# print(f.read())
# f.close()


# task 2 question 1 - find the name of the first team mentioned in each article

teams = []
with open ('rugby.json') as f:
    data = json.load(f)

   # find a list of team names in the json file
    for p in data['teams']:
        team_name = list(p.items())[0][1]
        teams.append(team_name)
   

# extract the name of the first team mentioned in the article [this means both the headline and the article!]

# a dictionary of url vs (headline and first team mentioned)
articles = defaultdict(list)

for url in visited.keys():
    page = requests.get(url)
    soup = BeautifulSoup(page.content, 'html.parser')
    for txts in soup.find_all('p'):
        txt = txts.get_text()
        headline = soup.find(id='headline').get_text()
        
        # if team name found in headline, then exit
        for word in list(headline.split()):
            
            # RE used to remove any non-word characrters
            word_list = re.split(r'\W+', word)
            for w in word_list:
                
                # for each token, if they are team names, then add them
                if w in teams:
                    articles[url].append(visited[url])
                    articles[url].append(w)
                    break
                
        # if not, continue iterate through the main body part 
        for word in list(txt.split()):
            word_list = re.split(r'\W+', word)
            for w in word_list: 
                if w in teams:
                    articles[url].append(visited[url])
                    articles[url].append(w)
                    break
                
             
                

# task 2 question 2 - find the largest match score for each team

# a dictionary of team name vs all match scores
lscore = defaultdict(list)

for url in articles.keys():
    page = requests.get(url)
    soup = BeautifulSoup(page.content, 'html.parser')
    scores = []
    
    # iterate through the article body 
    for txts in soup.find_all('p'):
        txt = txts.get_text()
        
        # RE used to search for valid scores within 3 digits 
        pattern = r'^\d{1,3}-\d{1,3}$'
        for word in re.findall(r"\w+(?:[-']\w+)*|'|[-.(]+|\S\w*", txt):
            if re.search(pattern,word):
                scores.append(word)
                
    scores.sort()
    if scores:
        
        match_scores = []
        for score in scores:
            
            # find the total score of each match score
            scs = re.split("-", score)
            first = int(scs[0])
            second = int(scs[1])
            match_score = first + second
            match_scores.append([match_score, score])
            
        # sort based on the total score
        match_scores.sort()
        
        # add the article headline, name of first mentioned team and largest score respectively
        lscore[url].append(articles[url][0])
        lscore[url].append(articles[url][1])
        
        max_total_score = match_scores[-1][0]
        score_list = []
        for x,y in match_scores:
            if x == max_total_score:
                scs = re.split("-", y)
                first = int(scs[0])
                second = int(scs[1])
                game_differences = abs(first - second)
                score_list.append([game_differences, y])
        score_list.sort()  
        lscore[url].append(score_list[-1][-1])
        

# write those results into an output csv file
task2 = open("Task2.csv", "w") 
headings = ['url', 'headline', 'team', 'score']
rows = [[x, y[0], y[1], y[2]] for x, y in lscore.items()]
task2 = open("Task2.csv", 'w')
writer = csv.writer(task2)
writer.writerow(headings)
writer.writerows(rows)
task2.close()

# read the file
fi = open("Task2.csv", "r") 
# print(fi.read())
fi.close()   

# task 3 - find the average match score difference for each team

game_diff = defaultdict(list)


# iterate through match scores for each team, and find their differences
for url in lscore.keys():
    score = lscore[url][2]
    team = lscore[url][1]
    scs = re.split("-", score)
    first = int(scs[0])
    second = int(scs[1])
    score_diff = abs(first-second) 
    game_diff[team].append(score_diff)

# a dictionary of team vs average score difference
av_game_diff = defaultdict(int)


# find the average match score difference for each team 
for team in game_diff.keys():
    avgamediff = sum(game_diff[team])/len(game_diff[team])
    av_game_diff[team] = avgamediff

    
# write the average match score differences into an output csv
task3 = open("Task3.csv", "w") 
headings = ['url', 'headline', 'team', 'score']
rows = [[x, y] for x, y in av_game_diff.items()]
task2 = open("Task3.csv", 'w')
writer = csv.writer(task3)
writer.writerow(headings)
writer.writerows(rows)
task3.close()

# read the file
fi = open("Task3.csv", "r") 
# print(fi.read())
fi.close()

# task 4: graph the five teams that articles are most frequently written about and the number of times an article is written about that team.

team_freq = defaultdict(int) # team vs freqeuency of mentioned 
 
teams_ment = defaultdict(list)  # url vs teams mentioned

# extract team names from articles using RE
for url in articles.keys():
    page = requests.get(url)
    soup = BeautifulSoup(page.content, 'html.parser')
    for txts in soup.find_all('p'):
        txt = txts.get_text()
        for word in list(txt.split()):
            word_list = re.split(r'\W+', word)
            for w in word_list: 
                if w in teams:
                    teams_ment[url].append(w)
                    
# count the number of times that each team is mentioned                    
for url in teams_ment.keys():
    teams = teams_ment[url]
    for te in teams:
        team_freq[te] += 1
    

# find the five teams that articles are most frequently written about 

# sort the dictionary based on the frequency of teams mentioned 
new_team_freq = [[y,x] for x,y in team_freq.items()]
new_team_freq.sort()


# find the five teams with greatest frequency 
teams = [x[1] for x in new_team_freq[-5:]]

# find the number of times that they are mentioned
freq = [x[0] for x in new_team_freq[-5:]]

# plot a bar graph
fig = plt.figure(figsize=(10,5))

plt.bar(teams, freq, width = 0.4)
plt.xlabel("Teams Mentioned")
plt.ylabel("Number of Times Mentioned")
plt.title("Top Five Frequently Mentioned Articles")
plt.ylim(0 * 25, 5 * 80)
plt.show()

# task 5 

# plot a scatter plot showing relationship between average game difference and n. of articles mentioned about it

# a dictionary of team vs n. of articles written about it
team_art = defaultdict(int)

for url in lscore.keys():
    team  = lscore[url][1]
    team_art[team] += 1

teams = list(team_art.keys())
y = team_art.values()
x = [av_game_diff[team] for team in teams]

scale_factor = 5

xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()

plt.xlim(0 * scale_factor, 5 * scale_factor)
plt.ylim(0 * scale_factor, 5 * scale_factor)

plt.scatter(x, y, marker='o')

plt.xlabel("Average Game Difference")
plt.ylabel("Number of Articles Written about the team")
plt.grid(True)
plt.title("Number of Articles Written about Each Team vs The Average Game Difference")
plt.show()


