import subprocess
import re
import matplotlib.pyplot as plt
import numpy as np

def interface(file):
    #it is better, given the enormous stream this creates, to write to a file and use gnuplot.
    x_vals_mat = []
    psi_vals_mat = []
    x_vals = []
    psi_vals = []
    try:
        out = subprocess.run([file], capture_output=True, text=True, check=True) #this runs the file, the rest is a parser
        
        output_lines = out.stdout.strip().split('\n') # remove whitespace and chop by line
        pattern = re.compile(r"x = (.*?),\s*psi = (.+)$") # god-awful regex that strips the x and psi values out.
        for line in output_lines:
            line = line.strip()
            if line == "-": # HACK for multiple states, very janky and will break if you blink near it.
                x_vals_mat.append(np.array(x_vals))
                psi_vals_mat.append(np.array(psi_vals))
                x_vals = []
                psi_vals = []
            match = pattern.match(line)
            if match:
                x=float(match.group(1))
                psi=float(match.group(2))
                x_vals.append(x)
                psi_vals.append(psi)
    except Exception as e:
        print(f"bad stuff happened, {e}")
    return x_vals_mat, psi_vals_mat

x_list, psi_list = interface("PATH") # replace with a path to your binary, for most reliable results, use absolute path.
x = np.array(x_list)
psi = np.array(psi_list)

# plot results: 
n_rows = x.shape[0]
n_cols = x.shape[1]

fig, ax= plt.subplots(n_rows, 1, figsize=(15, 4*n_rows))
if n_rows == 1: # HACK for if only one row
    ax = [ax]
for i in range(n_rows):
    ax[n_rows-i-1].plot(x[i, :], psi[i, :], label=f'state {i}') # TODO change label for e.g. 0 = ground state.
    #plt.show()
