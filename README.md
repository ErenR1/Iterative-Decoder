# Iterative-Decoder

Given two files are the working prototypes for the LDPC(Low Density Parity Check) code's decoder.
One example is uploaded in a txt file ("H1.txt") which is nothing but the indices where the sparse matrix H[i][j] = 1;

Example.txt file contains the input code and corresponding output


Both use similar techinique but a small varitaion in space complexity.


As in both cases the preferred channel was BEC (Binary Errasure Channel) with a bit-flip probability of "e". So that the "Log Likelihood Ratio" (LLR) is user defined.

Using binary addition (XOR) and Tanh rule we would change the information stored in the graph and update it for a number of iterations (usually between 30-90 times). 
Based on the final numbers and their sign we can conclude the decoded values of the codes at the receiver end.

Please Note that these are only prototypes the decoder for Berman Codes and their performance analysis are under work.
