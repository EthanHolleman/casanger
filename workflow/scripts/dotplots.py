def plot_dot_plot(axis, seq):
    
    def delta(x,y):
        return 0 if x == y else 1
    
    def M(seq1,seq2,i,j,k):
        return sum(delta(x,y) for x,y in zip(seq1[i:i+k],seq2[j:j+k]))


    def makeMatrix(seq,k):
        n = len(seq)
        m = len(seq)
        return [[M(seq,seq,i,j,k) for j in range(m-k+1)] for i in range(n-k+1)]

    return np.array(makeMatrix(seq,1))
   # axis.xticks(np.arange(len(list(seq))),list(seq))
   # axis.yticks(np.arange(len(list(seq))),list(seq))
    
    return dotplot