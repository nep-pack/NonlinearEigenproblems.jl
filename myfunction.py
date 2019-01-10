

import numpy as np;
import cmath as m;
def compute_M(s):
    """Compute the matrix M(s) for a given eigenvalue approximation"""
    A=np.matrix('1 2; 3 4')
    B=np.matrix('0 0; 0 1');
    C=np.matrix('1 1; 1 1');
    M=A+s*B+m.exp(s)*C
    return M

def compute_Mlincomb(s,X):
    """Compute the linear combination of derivatives for value s"""
    A=np.matrix('1 2; 3 4')
    B=np.matrix('0 0; 0 1');
    C=np.matrix('1 1; 1 1');
    X=np.matrix(X)
    print(X)
    print("size(X,1)=",np.size(X,1))

    z=np.zeros((2,1));
    print("size(z,0)=",np.size(z,0))
    print("size(z,1)=",np.size(z,1))
    z=z+A*X[:,0]
    z=z+B*(s*X[:,0])
    z=z+C*(m.exp(s)*X[:,0])

    if (np.size(X,1)>1):
        z=z+B*(X[:,1])+C*(m.exp(s)*X[:,1])

    if (np.size(X,1)>1):
        for k in range(2,np.size(X,1)):
            print("k=",k)
            z=z+C*(m.exp(s)*X[:,k])
    return z




print(compute_M(3+3j))
X=np.matrix('1 2 -1 0 ; 3 4 -2 1')
print(compute_Mlincomb(3,X))
X2=np.matrix(' 1.+0.j  1.+0.j  1.+0.j  1.+0.j  1.+0.j; 1.+0.j  1.+0.j  1.+0.j  1.+0.j  1.+0.j')
z=compute_Mlincomb(3,X2);
print(z)
