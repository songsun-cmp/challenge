"""
    ncon(tensor_list, connect_list_in; cont_order=[], check_network=true)
------------------------
    by Glen Evenbly (c) for www.tensors.net, (v1.2) - last modified 19/1/2019
------------------------
Network CONtractor. Input is an array of tensors 'tensor_list' and an array of
vectors 'connect_list_in',  with each vector labelling the indices of the
corresponding tensor. Labels should be  positive integers for contracted indices
and negative integers for free indices. Optional input 'cont_order'  can be used
to specify order of index contractions (otherwise defaults to ascending order of
the positive indices). Checking of the consistancy of the input network can be
disabled for slightly faster operation.

Further information can be found at: https://arxiv.org/abs/1402.0939
"""
function ncon(tensor_list, connect_list_in; cont_order=[], check_network=true)

    # copy original list to avoid destroying
    connect_list = deepcopy(connect_list_in);

    # put inputs into an array if necessary
    if (tensor_list[1] isa Real) | (tensor_list[1] isa Complex)
        tensor_list = Any[tensor_list];
    end
    if !(connect_list[1] isa Array)
        connect_list = Any[connect_list];
    end

    # generate contraction order if necessary
    flat_connect = vcat(connect_list...);
    if isempty(cont_order)
       cont_order = sort(unique(flat_connect[flat_connect .> 0]));
    end

    # check inputs if enabled
    if check_network
       dims_list = Array{Any,1}(undef,length(tensor_list));
       for ik = 1:length(tensor_list)
           dims_list[ik] = [size(tensor_list[ik])...];
       end
       ncon_check_inputs(connect_list, flat_connect, dims_list, cont_order);
    end

    # do all partial traces
    for ip = 1:length(connect_list)
        num_cont = length(connect_list[ip]) - length(unique(connect_list[ip]));
        if num_cont > 0
            tensor_list[ip], connect_list[ip], cont_label = ncon_partial_trace(tensor_list[ip], connect_list[ip]);
            cont_order = setdiff(cont_order,cont_label);
        end
    end

    # do all binary contractions
    while !isempty(cont_order)
        # identify tensors to be contracted
        cont_ind = minimum(cont_order);
        locs = [ele for ele in collect(1:length(connect_list)) if sum(connect_list[ele] .== cont_ind) > 0];

        # do a binary contraction
        cont_many = intersect(connect_list[locs[1]],connect_list[locs[2]]);
        A_cont = [findfirst(connect_list[locs[1]] .== x) for x in cont_many];
        B_cont = [findfirst(connect_list[locs[2]] .== x) for x in cont_many];
        A_free = deleteat!(collect(1:length(connect_list[locs[1]])), sort(A_cont));
        B_free = deleteat!(collect(1:length(connect_list[locs[2]])), sort(B_cont));
        push!(tensor_list, ncon_tensordot(tensor_list[locs[1]],tensor_list[locs[2]],A_cont,B_cont));
        push!(connect_list, vcat(connect_list[locs[1]][A_free],connect_list[locs[2]][B_free]));

        # remove contracted tensors from list and update cont_order
        deleteat!(connect_list,locs);
        deleteat!(tensor_list,locs);
        cont_order = setdiff(cont_order,cont_many);
    end

    # do all outer products
    while length(tensor_list) > 1
        s1 = size(tensor_list[end-1]);
        s2 = size(tensor_list[end]);
        tensor_list[end-1] = reshape(reshape(tensor_list[end-1],prod(s1))*reshape(tensor_list[end],1,prod(s2)),(s1...,s2...));
        connect_list[end-1] = vcat(connect_list[end-1],connect_list[end]);
        deleteat!(connect_list,length(connect_list));
        deleteat!(tensor_list,length(tensor_list));
    end

    # do final permutation
    if length(connect_list[1]) > 0
        return permutedims(tensor_list[1],sortperm(connect_list[1],by=abs));
    else
        return tensor_list[1][1];
    end
end

"""
ncon_tensordot: contracts a pair of tensors via matrix multiplication, similar to the
Numpy function of the same name
"""
function ncon_tensordot(A,B,A_cont,B_cont)

    A_free = deleteat!(collect(1:ndims(A)), sort(A_cont));
    B_free = deleteat!(collect(1:ndims(B)), sort(B_cont));
    A_perm = vcat(A_free, A_cont);
    B_perm = vcat(B_cont, B_free);

    return reshape(reshape(permutedims(A,A_perm),prod(size(A)[A_free]),prod(size(A)[A_cont]))*
        reshape(permutedims(B,B_perm),prod(size(B)[B_cont]),prod(size(B)[B_free])),(size(A)[A_free]...,size(B)[B_free]...))
end

"""
ncon_partial_trace: partial trace on tensor A over repeated labels in A_label
"""
function ncon_partial_trace(A,A_label)

    num_cont = length(A_label) - length(unique(A_label));
    if num_cont > 0
        dup_list = []
        for ele in unique(A_label)
            if sum(A_label .== ele) > 1
                dup_list = vcat(dup_list, findall(A_label .== ele));
            end
        end

        cont_ind = reshape(permutedims(reshape(dup_list,2,num_cont),[2,1]),2*num_cont);
        free_ind = deleteat!(collect(1:length(A_label)),sort(dup_list));
        cont_dim = prod(size(A)[cont_ind[1:num_cont]]);
        free_dim = size(A)[free_ind];

        cont_label = unique(A_label[cont_ind]);
        B = zeros(prod(free_dim));
        perm_tot = [free_ind; cont_ind];

        A_dims = size(A);
        A = reshape(permutedims(reshape(A[:],A_dims),vcat(free_ind, cont_ind)),prod(free_dim),cont_dim,cont_dim);
        for ip = 1:cont_dim
            B = B + A[:,ip,ip];
        end

        return reshape(B,free_dim), deleteat!(A_label, sort(cont_ind)), cont_label
    else
        return A, A_label
    end
end

"""
ncon_check_inputs: check consistency of input tensor network
"""
function ncon_check_inputs(connect_list, flat_connect, dims_list, cont_order)

    pos_ind = flat_connect[flat_connect .> 0];
    neg_ind = flat_connect[flat_connect .< 0];

    # check that lengths of lists match
    if length(dims_list) != length(connect_list)
        error("NCON error: $(length(dims_list)) tensors given but $(length(connect_list)) index sublists given")
    end

    # check that tensors have the right number of indices
    for ik = 1:length(dims_list)
        if length(dims_list[ik]) != length(connect_list[ik])
            error("NCON error: number of indices does not match number of labels on tensor $(ik): $(length(dims_list[ik]))-indices versus $(length(connect_list[ik]))-labels")
        end
    end

    # check that contraction order is valid
    if !(sort(cont_order) == sort(unique(pos_ind)))
        error("NCON error: invalid contraction order")
    end

    # check that negative indices are valid
    for ind = -1:-1:-length(neg_ind)
        if sum(neg_ind .== ind) == 0
            error("NCON error: no index labelled $(ind)")
        elseif sum(neg_ind .== ind) > 1
            error("NCON error: more than one index labelled $(ind)")
        end
    end

    # check that positive indices are valid and contracted tensor dimensions match
    flat_dims = [];
    for ele in dims_list
        flat_dims = vcat(flat_dims,ele)
    end
    for ind in unique(pos_ind)
        if sum(pos_ind .== ind) == 1
            error("NCON error: only one index labelled $(ind)")
        elseif sum(pos_ind .== ind) > 2
            error("NCON error: more than two indices labelled $(ind)")
        end
        cont_dims = flat_dims[flat_connect .== ind];
        if cont_dims[1] != cont_dims[2]
            error("NCON error: tensor dimension mismatch on index labelled $(ind): dim-$(cont_dims[1]) versus dim-$(cont_dims[2])")
        end
    end

    return true
end