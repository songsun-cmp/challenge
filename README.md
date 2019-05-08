##Some explaination for the "challenge.jl"
#Basic Run procedure
julia> using LinearAlgebra
julia> include("ncon.jl")
julia> include("challenge.jl")
julia> fed(1.0) #This is the ans1
julia> for beta in 1:10 printl(degeneracy(beta)) end #the smallest one is the ans2
#Collaborator
张华琛 孙松 王昊昕
