import math
import sqlite3
from random import random
import sys

# Create or connect to SQLite database
conn = sqlite3.connect('points.db')
cursor = conn.cursor()

# Create a table for points on the unit circle
cursor.execute('''
CREATE TABLE IF NOT EXISTS unit_circle_points (
    id INTEGER PRIMARY KEY,
    x REAL,
    y REAL
)
''')
cursor.execute('''
DELETE FROM unit_circle_points
''')
conn.commit()

num_points = 10

if (len(sys.argv) == 2 and sys.argv[1]):
    num_points = int(sys.argv[1])

# Insert points into the table
for i in range(num_points):
    x = random() * 0.5
    y = random() * 0.5
    norm = math.sqrt(x**2 + y**2)
    if (norm == 0 or norm > 1):
        conn.close()
        print("Point not on the unit circle")
        break
    cursor.execute('''
    INSERT INTO unit_circle_points (x, y)
    VALUES (?, ?)
    ''', (x, y))
    conn.commit()

# Commit and close the connection
conn.close()

# print("SQLite database with points on the unit circle created.")
