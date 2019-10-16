#  This file shows some example tests of the library
#
push!(LOAD_PATH, "D:/ettore/Git/dual/julia")
# Start the test
using DualNumbers, LinearAlgebra, SparseArrays, SuiteSparse
using dual_LinearAlgebra


####################  dual_vec ####################
p=[0,0,0]
v=[1,0,0]
dv=dual_vec(p,v)
println("Dual vector")
test=display(dv)
####################  dual_pinv ####################
    n = 3
    A = randn(n,n)     # real-valued │
    A0 = randn(n,n)     # real-valued ├─ vector
    dA = A + ε * A0    # dual-valued │
    println("Dual Matrix dA=")
    display(dA)
    println("Real part of dA=")
    display(realpart.(dA))
    println("Dual part of dA=")
    display(dualpart.(dA))
    println("Moore-Penrose pseudoinverse of dA=")
    dG=dual_pinv(dA)
    println("Test correctness")
    test=display(dG*dA)
###### Test1
###### Test2
<<<<<<< HEAD

=======
###### Branch Ettore test1
>>>>>>> 9bb18e289f1ee0aa7e70a85e28a487023b4c59f8
