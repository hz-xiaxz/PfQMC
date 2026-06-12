using XDiag
using LinearAlgebra

function check_capabilities()
    N = 4
    block = Spinhalf(N)
    ops = OpSum()
    ops += Op("Sz", 1)
    _, psi = eig0(ops, block)
    
    println("Checking capabilities for psi::", typeof(psi))
    
    # 1. Check iteration/conversion
    try
        v = Vector(psi)
        println("Vector(psi) works. Eltype: ", eltype(v))
    catch e
        println("Vector(psi) failed.")
    end
    
    # 2. Check indexing
    try
        val = psi[1]
        println("psi[1] works. Value: ", val)
    catch e
        println("psi[1] failed.")
    end
    
    # 3. Check in-place apply
    try
        op = Op("Sz", 1)
        out = similar(psi)
        # Try different signatures for in-place apply if they exist in common libraries
        # apply!(out, op, psi)
        # XDiag might not have it, but let's check if 'apply!' is defined
        if isdefined(XDiag, :apply!)
            println("XDiag.apply! exists!")
        else
            println("XDiag.apply! does NOT exist.")
        end
    catch e
        println("In-place check failed: ", e)
    end
end

check_capabilities()
