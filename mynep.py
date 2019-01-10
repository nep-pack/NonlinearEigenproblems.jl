import numpy as np;
import cmath as m;

def compute_M(s):
    """Compute the matrix M(s) for a given eigenvalue approximation"""
    A=np.matrix('1 2; 3 4');  B=np.matrix('0 0; 0 1');   C=np.matrix('1 1; 1 1');
    M=A+s*B+m.exp(s)*C
    return M

def compute_Mlincomb(s,X):
    """Compute the linear combination of derivatives for value s"""
    A=np.matrix('1 2; 3 4');  B=np.matrix('0 0; 0 1');   C=np.matrix('1 1; 1 1');

    X=np.matrix(X) # Explicitly convert to matrix
    z=np.zeros((2,1));
    # Zeroth derivative
    z=z+A*X[:,0]
    z=z+B*(s*X[:,0])
    z=z+C*(m.exp(s)*X[:,0])

    # First derivative
    if (np.size(X,1)>1):
        z=z+B*(X[:,1])+C*(m.exp(s)*X[:,1])
    # Higher derivatives
    if (np.size(X,1)>1):
        for k in range(2,np.size(X,1)):
            z=z+C*(m.exp(s)*X[:,k])
    return z
