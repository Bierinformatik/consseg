## SLOW R IMPLEMENTATION
## TODO: compare incremental and slow R implementations
## TODO: relative error in test function, and count deviations
## potential


## print breakpoints
cons_r$breakpoints


### TEST MAIN WRAPPER
consensus(b,n=50,w=w)

## test compiling potential function
e <- "exp(L/2)-1)"
e <- "L*L*L/3"
consensus(b,n=50,w=w,e=e)

## test pre-compiled potential function
ec <- compileEquation(e)
consensus(b,n=50,w=w,e=ec)

## check with direct calculation above
cons_c$breakpoints
cons_r$breakpoints
