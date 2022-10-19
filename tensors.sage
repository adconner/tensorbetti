def tensors22c(c,tensor=True):
    Ts = []
    Ts.append( { (0,0,0) : 1 } )
    Ts.append( { (0,0,0) : 1, (1,1,0) : 1} )
    Ts.append( { (0,0,0) : 1, (0,1,1) : 1} )
    Ts.append( { (0,0,0) : 1, (1,0,1) : 1} )
    Ts.append( { (0,0,0) : 1, (0,1,1) : 1, (1,0,1) : 1} )
    Ts.append( { (0,0,0) : 1, (1,1,1) : 1} )
    Ts.append( { (0,0,0) : 1, (0,1,2) : 1, (1,0,1) : 1} )
    Ts.append( { (0,0,0) : 1, (0,1,1) : 1, (1,0,1) : 1, (1,1,2) : 1 })
    if c > 3:
        Ts.append( { (0,0,0) : 1, (0,1,2) : 1, (1,0,1) : 1, (1,1,3) : 1 })
    if tensor:
        return [ dict_to_tensor(2,2,c,T) for T in Ts ]
    else:
        return Ts
        
def tensors23c(c,tensor=True):
    Ts = tensors22c(c,False)
    Ts.append( { (0,0,0) : 1, (0,1,1) : 1, (0,2,2) : 1} )
    Ts.append( { (0,0,0) : 1, (0,2,1) : 1, (1,1,0) : 1} )
    Ts.append( { (0,0,0) : 1, (0,1,1) : 1, (1,1,0) : 1, (1,2,1) : 1} )
    Ts.append( { (0,0,0) : 1, (0,1,2) : 1, (1,0,1) : 1, (1,2,2) : 1} )
    Ts.append( { (0,0,0) : 1, (0,1,1) : 1, (1,2,2) : 1} )
    Ts.append( { (0,0,0) : 1, (0,1,1) : 1, (0,2,2) : 1, (1,0,1) : 1} )
    Ts.append( { (0,0,0) : 1, (0,1,1) : 1, (0,2,2) : 1, (1,0,1) : 1, (1,1,2) : 1} )
    Ts.append( { (0,0,0) : 1, (0,1,1) : 1, (1,0,1) : 1, (1,2,2) : 1} )
    Ts.append( { (0,0,0) : 1, (0,1,1) : 1, (1,1,1) : 1, (1,2,2) : 1} )
    if tensor:
        return [ dict_to_tensor(2,3,c,T) for T in Ts ]
    else:
        return Ts

def dict_to_tensor(a,b,c,Tdict):
    T = [matrix(ZZ,b,c) for i in range(a)]
    for (i,j,k),e in Tdict.items():
        T[i][j,k] += e
    return T
