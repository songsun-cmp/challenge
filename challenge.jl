include("ncon.jl");
#using LinearAlgerba
function tensor3(beta)
    tensor2(beta)=[exp(-beta) exp(beta); exp(beta) exp(-beta)]; #for AFM Ising model
    eig=eigen(tensor2(beta));
    half=diagm(0=>.√(complex(eig.values)))*eig.vectors;
    #transpose(half)*half-tensor2(beta)#test
    tmp=zeros(2,2,2);
    for ii=1:2
        for jj=1:2
	   for kk=1:2
	      tmp[ii,jj,kk]=half[ii,1]*half[jj,1]*half[kk,1]+half[ii,2]*half[jj,2]*half[kk,2];
	   end
	end
    end
    #tmp=ncon([half,half,half],[[-1,1],[-2,1],[-3,1]];cont_order=[1],check_network=true)
    return tmp
end

function tensor5(beta)
    T3=tensor3(beta);
    ncon(Any[T3,T3,T3,T3,T3],Any[[-1,1,5],[-2,2,1],[-3,3,2],[-4,4,3],[-5,5,4]])
end

function tensor9(beta)
    T5=tensor5(beta);
    ncon(Any[T5,T5,T5],Any[[-1,-2,2,1,-9],[-3,-4,-5,3,2],[-6,-7,-8,1,3]]);
end

function tensor12(beta)
    T9=tensor9(beta);
    ncon(Any[T9,T9],Any[[-1,1,2,3,-8,-9,-10,-11,-12],[-7,3,2,1,-2,-3,-4,-5,-6]]);
 end

function partition(beta)
    T12=tensor12(beta);
    ncon(Any[T12,T12],Any[[1,2,3,4,5,6,7,8,9,10,11,12],[3,2,1,12,11,10,9,8,7,6,5,4]]);
end

fed(beta)=-(1/(60*beta))*log(partition(beta));
fe(beta)=-(1/(beta))*log(partition(beta));
delta=1e-5;
entropy(beta)=((fe(beta+0.5*delta)-fe(beta-0.5*delta))/delta)*beta^2;
degeneracy(beta)=exp(entropy(beta));

#ncon([pauli[:,:,1],pauli[:,:,3]],[[-2,1],[1,-1]];cont_order=[1],check_network=true) #recording
#using PyCall
#@pyimport numpy as np


