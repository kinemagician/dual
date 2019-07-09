module dual_LinearAlgebra
using DualNumbers, LinearAlgebra, SparseArrays, SuiteSparse
export dual_pinv
export dual_vec
export dual_decodevec
#
function dual_vec(p,v)
#  DUAL_VEC
#  Form a dual vector.
#  Dual vector expresses a vector with a given line of action
#
#    INPUT
# p: Array with the coordinates of a point on the line of action
#    OUTPUT
# v: Array with vector components
    cc=cross(p,v)
    dual_v=v+ε*cc
    return dual_v
end # of function dual_vec

function dual_decodevec(dv)
#  DUAL_DECODEVEC
#  Decode a dual vector.
#  Dual vector expresses a vector with a given line of action
#
#    INPUT
#  dv
#    OUTPUT
# p: point on line vector closest to origin
# v: vector components
    v=realpart.(dv)
    u=dualpart.(dv)
    A=[0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
    p=-pinv(A)*u
    return p,v
end # of function dual_decodevec


function dual_pinv(dA)
#  DUAL_PINV
#  Compute the Moore-Penrose pseudoinverse of a dual matrix
#  References: The Moore–Penrose Dual Generalized Inverse Matrix With Application
#  to Kinematic Synthesis of Spatial Linkages
# July 2018Journal of Mechanical Design 140(10)
# DOI: 10.1115/1.4040882
#    INPUT
# dA: A matrix of dual elements
#    OUTPUT
# dG: Dual Moore-Penrose dual matrix
# Separate real and dual parts
    A=realpart.(dA)
    A0=dualpart.(dA)
# Get the size of the matrix
    m,n=size(A)
    G=pinv(A)
    At=transpose(A)
    AP1=pinv(At*A)
    Id=one(rand(m,m))
    G0=AP1*(transpose(A0)*(Id-A*G)-At*A0*G)
    dG=G+ε* G0
return dG
end  # End of function dual_pinv

end # End of module
