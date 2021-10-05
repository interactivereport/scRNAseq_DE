#!/usr/bin/env python

#Activate VIP env, such as: conda activate cellxgeneVIP

import sys
import sqlite3
import pandas as pd

db1 = sys.argv[1]
db2 = sys.argv[2]
tag1 = sys.argv[3]
tag2 = sys.argv[4]
prefix = sys.argv[5]

## db1
conn = sqlite3.connect(db1)
D1 = pd.read_sql_query("select * from DEG", conn)
conn.close()

## db2
conn = sqlite3.connect(db2)
D2 = pd.read_sql_query("select * from DEG", conn)
conn.close()

## add additional tags
if len(tag1)>0:
    D1['tags'] = D1['tags']+';'+tag1
if len(tag2)>0:
    D2['tags'] = D2['tags']+';'+tag2

df = pd.concat([D1,D2],ignore_index=True)
D = df[["log2fc","pval","qval"]]
D.index = pd.MultiIndex.from_frame(df[["gene","contrast","tags"]])

conn = sqlite3.connect(prefix+'.db')
D.to_sql("DEG",conn,if_exists="replace")
conn.close()

print("sqlite3 db (%s) with %d records was created successfully!"%(prefix+'.db',D.shape[0]))

#conn = sqlite3.connect("test.db")
#DB = pd.read_sql_query("select * from DEG;", conn)
#conn.close()


