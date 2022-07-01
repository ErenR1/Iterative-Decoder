# Belief Propagation for LDPC under BSC
import csv
import math
from random import choices
import numpy as np

with open('H1.txt', 'r') as H:
    sparse = csv.reader(H)
    sparse_data = list(sparse)

n = int(sparse_data[-1][1])
x = len(sparse_data)  # also, the number of edges
CN = {}
VN = {}
for i in range(x):
    if sparse_data[i][0] in CN:
        CN[sparse_data[i][0]].append(sparse_data[i][1])
    else:
        CN[sparse_data[i][0]] = [sparse_data[i][1]]
    if sparse_data[i][1] in VN:
        VN[sparse_data[i][1]].append(sparse_data[i][0])
    else:
        VN[sparse_data[i][1]] = [sparse_data[i][0]]
CN = {int(k): [int(i) for i in v] for k, v in CN.items()}
VN = {int(k): [int(i) for i in v] for k, v in VN.items()}


H_matrix = np.zeros((len(CN), n))
for i in range(len(sparse_data)):
    H_matrix[int(sparse_data[i][0]) - 1][int(sparse_data[i][1]) - 1] = 1

e = float(input("Enter Bit Flip Probability: "))
bi = [0, 1]
we = [e, 1 - e]
y = choices(bi, we, k=n)
Ly0 = math.log((1 - e) / e)
Ly1 = -Ly0
iter_n = int(input("Enter the number of iterations: "))

# Step 1 Initialization:

Vnval = np.zeros(len(VN))  # Info in ith VN
for key in VN:
    if y[key-1] == 0:
        Vnval[key-1] = Ly0
    else:
        Vnval[key-1] = Ly1
Cnval = np.zeros(len(CN))  # Info in ith CN

# Steps 2&3 Info Exchange:

for _ in range(iter_n):

    for u in CN:
        prod = 1
        for u1 in range(len(CN[u])):
            prod *= np.tanh(0.5 * Vnval[CN[u][u1]-1])
        Cnval[u-1] = prod

    for v in VN:
        if y[v-1] == 0:
            temp = Ly0
        else:
            temp = Ly1
        for v1 in range(len(VN[v])):
            t = np.tanh(0.5 * Vnval[v - 1])
            if t == 1 or t == -1:
                t = 0.98
            temp += 2 * np.arctanh(Cnval[VN[v][v1] - 1] / t)
        Vnval[v-1] = temp

# Step 4 Termination:
#print(Vnval)
#print(Cnval)
c_ha = np.zeros(n)
for i in range(len(Vnval)):
    if Vnval[i] < 0:
        c_ha[i] = 1
c_hat = c_ha.astype(int)
# print(y)
print(c_hat)  # Decoded value
print(np.matmul(c_hat, np.transpose(H_matrix)))  # Syndrome
