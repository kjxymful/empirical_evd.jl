using DynamicalSystems

# rhs for the bursting neuron model by Durstewitz 2009
@inline @inbounds function loop_burstn(u, p, t)
    I = p[1]
    Cₘ = p[2]
    gₗ = p[3]
    Eₗ = p[4]
    gₙₐ = p[5]
    Eₙₐ = p[6]
    Vₕₙₐ = p[7]
    kₙₐ = p[8]
    gₖ = p[9]
    Eₖ = p[10]
    Vₕₖ = p[11]
    kₖ = p[12]
    τₙ = p[13]
    gₘ = p[14]
    Vₕₘ = p[15]
    kₘ = p[16]
    τₕ = p[17]
    gₙₘ₀ₐ = p[18]
    Eₙₘ₀ₐ = 0 # as far as i could find

    V = u[1]
    n = u[2]
    h = u[3]
    s∞(V) = 1 / (1 + 0.33 * exp(-0.0625 * V))
    m∞(V, Vₕₙₐ, kₙₐ) = 1 / (1 + exp((Vₕₙₐ - V) / kₙₐ))
    n∞(V, Vₕₖ, kₖ) = 1 / (1 + exp((Vₕₖ - V) / kₖ))
    h∞(V, Vₕₘ, kₘ) = 1 / (1 + exp((Vₕₘ - V) / kₘ))

    du1 = (I - gₗ * (V - Eₗ) - gₙₐ * m∞(V, Vₕₙₐ, kₙₐ) * (V - Eₙₐ)
           -
           gₖ * n * (V - Eₖ) - gₘ * h * (V - Eₖ)
           -
           gₙₘ₀ₐ * s∞(V) * (V - Eₙₘ₀ₐ)) / Cₘ
    du2 = (n∞(V, Vₕₖ, kₖ) - n) / τₙ
    du3 = (h∞(V, Vₕₘ, kₘ) - h) / τₕ
    return SVector{3}(du1, du2, du3)
end


"""
implements the bursting neuron model by Durstewitz 2009

# Parameters:
-----------
- "gₙₘ₀ₐ::Float": the bifurcation parameter
- "u0::Vector{Flot}": initial condition
- every other parameter specified accessible as in Paper (capital subscript\n
                                         letters are small at variables)
...

# Returns:
--------
ds (ContinuousDynamicalSystem): dynamical system with the bursting neuron ODES
...
"""
function bursting_neuron(; u0=[-24.4694, 0.0386, 0.0231],
    I=0,
    Cₘ=6,
    gₗ=8,
    Eₗ=-80,
    gₙₐ=20,
    Eₙₐ=60,
    Vₕₙₐ=-20,
    kₙₐ=15,
    gₖ=10,
    Eₖ=-90,
    Vₕₖ=-25,
    kₖ=5,
    τₙ=1,
    gₘ=25,
    Vₕₘ=-15,
    kₘ=5,
    τₕ=200,
    gₙₘ₀ₐ=10.2)
    p = [I, Cₘ, gₗ, Eₗ, gₙₐ, Eₙₐ, Vₕₙₐ, kₙₐ, gₖ, Eₖ, Vₕₖ, kₖ, τₙ, gₘ, Vₕₘ, kₘ, τₕ, gₙₘ₀ₐ]
    ds = ContinuousDynamicalSystem(loop_burstn, u0, p)
    return ds
end