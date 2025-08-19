#!/usr/bin/env python3

import json
from seaborn import kdeplot
import matplotlib.pyplot as plt
import numpy as np

obj_1 = {
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

obj_2 = {
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
rw_rulings_1 = []
rw_rulings_2 = []

n = None

# Open the file and read it line by line
with open('/home/s24pjoha_hpc/goast_old_old/goast/examples/QuadMeshExamples/plotting/rw_rulings.txt', 'r') as file:
    n = sum(1 for _ in file)

    # Reset file pointer to the beginning
    file.seek(0)

    # Read the first n lines as positions
    for _ in range(n // 3):  # You should know n, the number of position vectors
        position_line = file.readline().strip()
        #numbers = [float(x) for x in position_line.split()]
        positions.append([float(x) for x in position_line.split()])

    # Read the second n lines as normal vectors
    for _ in range(2*(n // 3) ):  
        rw_ruling = file.readline().strip()
        if _ % 2 == 0:
            rw_rulings_1.append([float(x) for x in rw_ruling.split()])
        else:
            rw_rulings_2.append([float(x) for x in rw_ruling.split()])


    for i in range((n // 3)):
        obj_1["vectors"].append({
            "pos" : positions[i],
            "vec" : rw_rulings_1[i],
        })

    for i in range((n // 3)):
        obj_2["vectors"].append({
            "pos" : positions[i],
            "vec" : rw_rulings_2[i],
        })

with open("rw-rulings-field_1.json", "w") as writer:
    json.dump(obj_1, writer)



with open("rw-rulings-field_2.json", "w") as writer:
    json.dump(obj_2, writer)