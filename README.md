# Some explaination for the "challenge.jl"

## Basic Run procedure

###Question1

julia> include("challenge.jl")

julia> fed(1.0) #This is the ans1

###Question2
julia> x=zeros(100,1);

julia> y=zeros(100,1);

julia> for n=1:100
            y[n]=degeneracy(4+5/100*n,1e-6);
            x[n]=4+5/100*n;
       end

julia> using Plots

julia> plotly()

julia> plot(x,y,linewidth=2)

From the asymptoty of the line we get the degeneracy of the 

###Checks for our code 

julia> for n=1:100
                y[n]=fe(4+6.8/100*n);
                x[n]=4+6.8/100*n;
       end

julia> plot(x,y,linewidth=2)

From the plot we find that the Free energy is approching the accurate ground stat energy -66 when beta is increasing. The result make sence.

# Collaborator

张华琛 孙松 王昊昕
