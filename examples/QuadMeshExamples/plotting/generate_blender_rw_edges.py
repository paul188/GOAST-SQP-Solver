#!/usr/bin/env python3

import json
from seaborn import kdeplot
import matplotlib.pyplot as plt
import numpy as np

obj = {
    "min_len" : 0.02,
    "avg_len" : 0.04,
    "max_len" : 0.06,
    "colorscheme" : [
        {
            "value_to_map_to_this_color" : 0.75,
            "color_space" : "sRGB",
            "format" : "hex",
            "color" : "#DC143C",
            "kwargs" : {}
        },
        {
            "value_to_map_to_this_color" : 0.25 + 0.5 * 1.0 / 3.0,
            "color_space" : "sRGB",
            "format" : "hex",
            "color" : "#1474DC",
            "kwargs" : {}
        },
        {
            "value_to_map_to_this_color" : 0.25 + 0.5 * 2.0 / 3.0,
            "color_space" : "sRGB",
            "format" : "num",
            "color" : [ 1.0, 0.0, 1.0 ],
            "kwargs" : {}
        },
        {
            "value_to_map_to_this_color" : 0.25,
            "color_space" : "sRGB",
            "format" : "hex",
            "color" : "#80FFC8",
            "kwargs" : {}
        }
    ],
    "vectors": []
}

positions = []
rw_edges = []

n = None

# Open the file and read it line by line
with open('rw_edges.txt', 'r') as file:
    n = sum(1 for _ in file)

    # Reset file pointer to the beginning
    file.seek(0)

    # Read the first n lines as positions
    for _ in range(n // 2):  # You should know n, the number of position vectors
        position_line = file.readline().strip()
        #numbers = [float(x) for x in position_line.split()]
        positions.append([float(x) for x in position_line.split()])

    # Read the second n lines as normal vectors
    for _ in range(n // 2):  # The same n for the normals
        rw_edges_line = file.readline().strip()
        rw_edges.append([float(x) for x in rw_edges_line.split()])

for i in range(n // 2):
    obj["vectors"].append({
        "pos" : positions[i],
        "vec" : rw_edges[i],
    })

with open("rw_edges-field.json", "w") as writer:
    json.dump(obj, writer)