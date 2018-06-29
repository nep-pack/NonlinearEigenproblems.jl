import numpy as np

#n=4;   # degree
# coefficients in monomial basis
#a=np.random.random((n+1,1));
a=np.array([4,2,-1,5,-2]);
n=a.size-1; # degree

# coefficients in Chebishev basis
b=np.zeros(n+3);

rho=2.0;      # shift
gamma=1.0;    # rescale

# constants in the three term recurrence
alpha=1/(2*rho);    beta=-gamma/rho;

bb=np.zeros(n+3);   # bb are the coefficients in the next iteration

# remember that python do not take the right-most element in a loop
for j in range(n,-1,-1):
    # perform one iteration
    bb[0]=alpha*b[1]+beta*b[0]+a[j];
    bb[1]=beta*b[1]+alpha*b[2]+2*alpha*b[0];
    # remember that python do not take the right-most element in a loop
    for k in range(2,n-j-1):
        #print("Iteration:",j,"Inner loop:",k)
        bb[k]=alpha*b[k-1]+beta*b[k]+alpha*b[k+1];
    if n-j-1>1:
        bb[n-j-1]=alpha*b[n-j-2]+beta*b[n-j-1];
    if n-j>1:
        bb[n-j]=alpha*b[n-j-1];
    print("filling till:",n-j)
    print(bb)
    b=bb;
    bb=np.zeros(n+3);

print("Coefficients in monomial basis\n")
print(a)
print("Coefficients in Chebishev basis\n")
b=b[:5]
print(b)
