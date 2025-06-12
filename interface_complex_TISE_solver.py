import subprocess
import re
import ast
import matplotlib.pyplot as plt
import numpy as np

def interface(file):
    x_vals_mat = []
    psi_vals_mat = []
    x_vals = []
    psi_vals = []
    eigvals = []
    try:
        out = subprocess.run([file], capture_output=True, text=True, check=True) #this runs the file, the rest is a parser
        
        output_lines = out.stdout.strip().split('\n') # remove whitespace and chop by line
        pattern = re.compile(r"x = (.*?),\s*psi = (.+)$") # god-awful regex that strips the x and psi values out.
        eigval_regex = re.compile(r"^" + re.escape("eigval: ") + r"(.*)")
        for line in output_lines:
            line = line.strip()
            #print(line)
            #print("AA")
            if line == "-": # HACK for multiple states, very janky and will break if you blink near it.
                x_vals_mat.append(np.array(x_vals))
                psi_vals_mat.append(np.array(psi_vals))
                x_vals = []
                psi_vals = []
            match = pattern.match(line)
            #print(line)
            eigval = eigval_regex.match(line)
            #print("MATCH: ", match)
            if match:
                #print(match, match.group(0), match.group(1))
                if match.group(1)[0] == "(":
                    x=float(ast.literal_eval(match.group(1))[0])
                else:
                   x=float(match.group(1))
                if match.group(2)[0] == "(":
                    psi=float(ast.literal_eval(match.group(2))[0])
                else:
                    psi=float(match.group(2))
                x_vals.append(x)
                psi_vals.append(psi)
            if eigval:
                #print("VAL: ", ast.literal_eval(eigval.group(1))[0])
                
                eigvals.append(ast.literal_eval(eigval.group(1)))
    except Exception as e:
        print(f"bad stuff happened, {e}")
    #print(eigvals)
    return x_vals_mat, psi_vals_mat, eigvals

x_list, psi_list, eigvals_list = interface("/home/cj/git_uni/M2R/Quantum/code") # replace with a path to your binary, for most reliable results, use absolute path. Must be an exe for Windows.
x = np.array(x_list)
psi = np.array(psi_list)

# plot results: 
n_rows = x.shape[0]
n_cols = x.shape[1]

fig, ax= plt.subplots(n_rows, 1, figsize=(15, 4*n_rows))
if n_rows == 1: # HACK for if only one row
    ax = [ax]
for i in range(n_rows):
    ax[n_rows-i-1].plot(x[i, :], psi[i, :]) # TODO change label for e.g. 0 = ground state.
    #(eigvals_list[i], type(eigvals_list[i]))
    ax[n_rows-i-1].legend([str(eigvals_list[i][0])], loc = "upper left")
    #plt.show()
l1 = []
l2 = []
for i in eigvals_list:
    if i[0]<3000:
        l1.append(i)
    else:
        l2.append(i)

#plt.plot([i[0] for i in eigvals_list], [i[1] for i in eigvals_list])
#plt.plot([i[0] for i in l1], [i[1] for i in l1])
#plt.plot([i[0] for i in l2], [i[1] for i in l2])
#for i in l1:
#    print(i[1]/i[0])