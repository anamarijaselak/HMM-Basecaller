import numpy as np
import os

def get6mers():
    with open(os.getcwd()+"/r9_250bps.nucleotide.6mer.template.model", "r") as ins:
        array = []
        for line in ins:
            array.append(line.split()[0])
