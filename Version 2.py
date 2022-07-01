# Belief Propagation for LDPC codes under BSC channel with/without noise
import csv
import math
from random import choices
import numpy as np

 # IF needed AWGN Channel noise = np.random.normal(0,1,10) 

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
VN = {int(k): [int(i) for i in v] for k, v in VN.items()} # Creating Ajacency list type struct for Tanner Graph

CN_end = np.zeros(x, dtype='int16')
VN_end = np.zeros(x, dtype='int16')

index = 0
for j in CN:
    for k in range(len(CN[j])):
        index += 1
        VN_end[index-1] = CN[j][k]
        CN_end[index-1] = j

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

VN_message = np.zeros(x)
CN_message = np.zeros(x)

for _ in range(iter_n):

    for ed in range(x):
        i = VN_end[ed]
        dv = len(VN[i])
        if y[i-1] == 1:
            temp = Ly1
        else:
            temp = Ly0

        for idx in range(dv):
            edge_idx = VN[i][idx]
            if edge_idx != ed:
                temp += CN_message[edge_idx]
        VN_message[ed] = temp

    for ed in range(x):
        j = CN_end[ed]
        dc = len(CN[j])
        prod = 1
        for idx in range(dc):
            edge_indx = CN[j][idx]
            if edge_indx != ed:
                prod *= 2 * np.arctanh(VN_message[edge_indx])
        CN_message[ed] = prod

c_hat = np.zeros(n)
for i in range(n):
    if y[i] == 0:
        si = Ly0
        for q in range(len(VN[i])):
            si += CN_message[VN[i][q]]
    else:
        si = Ly1
        for q in range(len(VN[i])):
            si += CN_message[VN[i][q]]
    if si < 0:
        c_hat[i] = 1
    else:
        c_hat[i] = 0

print(np.matmul(c_hat, np.transpose(H_matrix)))  # Checking the syndrome

