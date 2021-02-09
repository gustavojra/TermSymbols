module TermSymbols
using Combinatorics

export Determinant
export get_dets
export latex_that_shit
export latex

orb_per_shell = Dict(
        "s" => 1,
        "p" => 3,
        "d" => 5,
        "f" => 7,
        "h" => 9,
        "g" => 11
)

function num_of_orb(Orb::String)

    if length(Orb) == 1
        return orb_per_shell[Orb]
    elseif length(Orb) == 2
        return orb_per_shell[Orb[2:2]]
    else
        error("Invalid orbital label: $Orb")
    end
end

struct Determinant
    Ml::Int
    Ms::Rational{Int64}
    α::Dict{String,BitArray{1}}
    β::Dict{String,BitArray{1}}
end

function Determinant(αstr::Dict{String,BitArray{1}}, βstr::Dict{String,BitArray{1}})
    
    @assert keys(αstr) == keys(βstr)

    Ms = 0//1
    Ml = 0

    for k in keys(αstr)

        α = αstr[k]
        β = βstr[k]

        Ms += (sum(α) - sum(β))//2

        Norb = length(α)

        L = Int((Norb - 1)/2)

        Angular = [i-L for i = 0:(Norb-1)]

        Ml += sum(Angular[α]) + sum(Angular[β])

    end

    return Determinant(Ml, Ms, αstr, βstr)

end

# Return a Tuple with all BitArrays representing determinants that can be formed for a given shell
function get_dets_in_a_shell(Nα::Int, Nβ::Int, Orb::String)

    Norb = num_of_orb(Orb)

    α = vcat(repeat([1],Nα),repeat([0],Norb-Nα))
    β = vcat(repeat([1],Nβ),repeat([0],Norb-Nβ))

    αperms = multiset_permutations(α, Norb)
    βperms = multiset_permutations(β, Norb)


    out = Tuple{BitArray{1},BitArray{1}}[]

    for p1 in αperms
        for p2 in βperms
            push!(out, (BitArray(p1), BitArray(p2)))
        end
    end

    return out
end

function get_dets(config::Dict{String, Tuple{Int,Int}})

    out = Determinant[]

    shells = Dict()

    for k in keys(config)

        Orb = k
        Nα = config[k][1]
        Nβ = config[k][2]

        shells[k] = get_dets_in_a_shell(Nα, Nβ, Orb)
    end

    # Assuming just one shell
    
    if length(config) == 1

        k, = keys(config)

        N = length(shells[k])

        for i = 1:N
            push!(out, Determinant(Dict(k=>shells[k][i][1]), Dict(k=>shells[k][i][2])))
        end

    elseif length(config) == 2

        k1,k2 = keys(shells)

        s1,s2 = shells[k1], shells[k2]

        for i = eachindex(s1)
            for j = eachindex(s2)
                αstr = Dict(k1=>s1[i][1],
                            k2=>s2[j][1])

                βstr = Dict(k1=>s1[i][2],
                            k2=>s2[j][2])

                push!(out, Determinant(αstr, βstr))
            end
        end
    else
        error("Three shells not available now sorry :(")
    end

    return out
end

function latex(D::Determinant, core="")

    out = core

    for k in keys(D.α)

        Norb = num_of_orb(k)

        maxL = Int((Norb - 1)/2)

        Angular = [i-maxL for i = 0:(Norb-1)]

        α = D.α[k]
        β = D.β[k]

        for i in eachindex(Angular)

            ml = Angular[i]

            α[i] ? out = out*k*raw"$_{"*"$ml"*raw"}^\alpha$ " : nothing

            β[i] ? out = out*k*raw"$_{"*"$ml"*raw"}^\beta$ " : nothing
        end
    end

    return out
end

function latex_that_shit(dets::Array{Determinant,1}, Ml::Int)

    dets = filter(d->d.Ml == Ml, dets)

    base = raw"$[Ar]$ 4s$^2$ " 
    
    dorbs = [raw"3d$_{-2}",raw"3d$_{-1}",raw"3d$_{0}",raw"3d$_{+1}",raw"3d$_{+2}"]
    #dorbs = [raw"3p$_{-1}",raw"3p$_{0}",raw"3p$_{+1}"]

    L = length(dets[1].α)

    out = ""
    for D in dets
        line = base
        for i = 1:L
            if D.α[i]
                line = line*dorbs[i]*raw"\alpha$ "
            end
            if D.β[i]
                line = line*dorbs[i]*raw"\beta$ "
            end
        end

        println(line*raw"\\\\")
        println(raw"\noindent")

    end

end

end # module
