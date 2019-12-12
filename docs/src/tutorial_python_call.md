# Tutorial: Calling NEP-PACK from python

## PyJulia

The previous tutorial illustrated how a NEP defined
in python code can be solved with NEP-PACK.
Python is a more mature (read: older) language than Julia,
and there are considerable packages and features
in python not present in Julia. If you need these features,
it may be more convenient to call NEP-PACK
from  python, rather than calling python from julia.

The python package [PyJulia](https://github.com/JuliaPy/pyjulia)
gives us that possibility.

Installation on Ubuntu linux:
```
$ python3 -m pip install julia # Only necessary first time you run it
...
$ python3
>>> import julia
>>> jl = Julia(compiled_modules=False) # compilation flag necessary on ubuntu
>>> julia.install()               # Only necessary first time you run it
>>> from julia import Base
>>> Base.MathConstants.golden  # Julia's definition of golden ratio
1.618033988749895
```


## Using PyJulia and NEP-PACK

import numpy as np;
import cmath as m;

def my_compute_M(s,der):
    """Compute the matrix M(s) for a given eigenvalue approximation"""
    A=np.matrix('1 2; 3 4');  B=np.matrix('2 0; 0 1');   C=np.matrix('1 1; 1 1');
    M=m.exp(s)*C
    if (der==0):
        M=M+A+s*B;
    elif (der==1):
        M=M+B;
    return M


![To the top](http://jarlebring.se/onepixel.png?NEPPACKDOC_PYTHON2)
