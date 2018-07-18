###print

using ArgParse
using HDF5

using PyCall
@pyimport h5py

# julia really doesn't want to read a 4-index tensor written to hdf5
# from C++ the right way. It's annoying. I should file a bug report.
# not type-stable!
#
function read_4index_cxx_hdf5(fn)
    f = h5py.File("$(fn)", "r+")
    μcc = get(f, "tensor")[:value]
    f[:close]()

    @show size(μcc)
    μcc = permutedims(μcc, [2,1,4,3])
    L = size(μcc,1)
    N = size(μcc,2)
    reshape(μcc, (N,L,N,L))
    if 5 == length(size(μcc))
        μcc  =  μcc[:,:,:,:,1] + im*μcc[:,:,:,:,2]
    elseif 4 == length(size(μcc))
    else
        error("wrong rank")
    end
    return μcc
end

#julia h5read is fundamentally broken here
function read_2index_cxx_hdf5(fn)
    #@show fn
    f = h5py.File("$(fn)", "r+")
    μcc = get(f, "tensor")[:value]
    f[:close]()
    if 3 == length(size(μcc))
        μcc = μcc[:,:,1] + im * μcc[:,:,2]
    elseif 2 != length(size(μcc))
        error("wrong rank")
    end
    return μcc
end

function pauli_matrices_sparse(L :: Int64)
    sigx = sparse([0 1; 1 0])
    sigy = sparse([0 -im; im 0])
    sigz = sparse([1 0; 0 -1])
    sigp = sparse([0 1; 0 0])
    sigm = sparse([0 0; 1 0])

    X = [reduce(kron, (speye(2^(j-1)), sigx, speye(2^(L - j)))) for j in 1:L]
    Y = [reduce(kron, (speye(2^(j-1)), sigy, speye(2^(L - j)))) for j in 1:L]
    Z = [reduce(kron, (speye(2^(j-1)), sigz, speye(2^(L - j)))) for j in 1:L]
    
    P = [reduce(kron, (speye(2^(j-1)), sigp, speye(2^(L - j)))) for j in 1:L]
    M = [reduce(kron, (speye(2^(j-1)), sigm, speye(2^(L - j)))) for j in 1:L]

    X = convert(Array{SparseMatrixCSC{Float64,Int64}, 1}, X)
    Y = convert(Array{SparseMatrixCSC{Complex{Float64},Int64}, 1}, Y)
    Z = convert(Array{SparseMatrixCSC{Float64,Int64}, 1}, Z)
    P = convert(Array{SparseMatrixCSC{Float64,Int64}, 1}, P)
    M = convert(Array{SparseMatrixCSC{Float64,Int64}, 1}, M)
    return (X,Y,Z,P,M)
end

function exact_μ{T}(H :: Array{T,2}, j :: Array{Float64,2}, N :: Int)
    μ = zeros(N,N)
    d,v = eig(H)
    H = zeros((2,2))
    gc()
    
    jv1 = copy(v)
    A_mul_B!(jv1,j,v)
    jv = Ac_mul_B(v,jv1)
    
    jv1 = zeros(2,2)
    v = zeros(2,2)
    gc()
    for n in 1:N
        TnH = cos.((n-1)*acos.(d))
        TnHjvT = transpose(TnH .* jv)
        for m in 1:N
            TmH = cos.((m-1)*acos.(d))
            μ[n,m] = μ[m,n]= vecdot(TnHjvT, (TmH .* jv))
        end
        TnHjvT = zeros(2,2)
        gc()
    end
    return μ
end

function correlation{T,S,R}(H :: Array{T,2}, A :: Array{S,2}, B :: Array{R,2}, N :: Int)
    μ = zeros(Complex{Float64},N,N)
    d,v = eig(H)
    H = zeros((2,2))
    gc()

    Av = v'*A*v
    Bv = v'*B*v
    
    for n in 1:N
        TnH = diagm(cos.((n-1)*acos.(d)))
        for m in n:N
            TmH = diagm(cos.((m-1)*acos.(d)))
            μ[n,m] = μ[m,n]= trace(TnH*Av*TmH*Bv)
        end
        gc()
    end
    return μ
end

function rfheis(hz)
    L = length(hz)
    X,Y,Z,P,M = pauli_matrices_sparse(L)
    H  = 0.25*sum(map((x1,x2) -> x1 * x2, X[1:end-1], X[2:end]))
    H += 0.25*sum(map((x1,x2) -> x1 * x2, Y[1:end-1], Y[2:end]))
    H += 0.25*sum(map((x1,x2) -> x1 * x2, Z[1:end-1], Z[2:end]))
    H += 0.5*sum(hz.*Z)
    H = H/(3*(L-1)*0.25 + sum(abs.(hz))*0.5)
    Hfull = H |> full |> real
    return Hfull
end

function rf_2NJW(w :: Array)
    
    #hardcode V
    
    # https://arxiv.org/abs/cond-mat/0610854
    # JW transform notes of 2017-12-21
    # factors of 2,4 are S/σ translation
    
    L = size(w,1)
    X,Y,Z,P,M = pauli_matrices_sparse(L)
    H  = sum(j -> (0.25+0im)*X[j]*X[j+1], 1:L-1)
    H += sum(j -> (0.25+0im)*Y[j]*Y[j+1], 1:L-1)
    H += sum(j -> (0.25+0im)*Z[j]*Z[j+1], 1:L-1)
    
    H += sum(j -> (0.25+0im)*X[j]*Z[j+1]*X[j+2], 1:L-2)
    H += sum(j -> (0.25+0im)*Y[j]*Z[j+1]*Y[j+2], 1:L-2)
    
    H += sum(j -> 0.5*w[j]*Z[j], 1:L )
    H = H/(3*0.25*(L-1) + 2*0.25*(L-2) + sum(abs.(w))*0.5)
    return H |> full
end

function check_dos_trace(H, ifn, ofn, L)
    tnf = open("$ifn.chtrre")
    trTncc = 2^(L/2)*[parse(Float64, s) for s in split(readline(tnf))]
    N = length(trTncc)
    close(tnf)
    d,v = H |> full |> eig
    trTned = zeros(N)
    for n = 1:N
        trTned[n] = sum(cos.((n-1)*acos.(d)))
    end
    #@show L
    #@show trTned[1:3]
    #@show trTncc[1:3]
    trTndiff = trTned - trTncc
    @show abs.(trTndiff) |> maximum
end


function check_conductivity(H, j, ifn,ofn,L)
    μcc = readdlm("$ifn.re")
    N = size(μcc,1)
    μed = -exact_μ(H,j,N) #minus because I took imag part of j

    cond_diff = 2.0^(-L)*μed - μcc
    #@show L
    #@show μed[1,1], μed[1,2]
    #@show 2.0^(L)*μcc[1,1], 2.0^(-L)*μcc[1,2]
    #@show L, (μed./μcc)[1,1], (μed./μcc)[1,2]
    @show abs.(cond_diff) |> maximum 
end

function check_Sz_correlation(H, ifn, L)
    X,Y,Z,P,M = pauli_matrices_sparse(L)
    Z = map(full, Z)

    #julia's h5read is broken. Need the reshape to fix.
    μcc = h5read("$(ifn)SzSz.h5", "/tensor")
    N = size(μcc, 2)
    μcc = reshape(μcc, (N,L,N,L))
    
    j1 = 1
    j2 = 1
    μed = correlation(H, 0.5*Z[j1], 0.5*Z[j2], N)
    Sz1Sz1_diff = 2.0^(-L)* μed - μcc[:,j1, :,j2]
    @show abs.(Sz1Sz1_diff) |> maximum

    μed = zeros(N, L, N, L)
    for j1 = 1:L
        for j2 = 1:L
            μed[:,j1,:,j2] = correlation(H, 0.5*Z[j1], 0.5*Z[j2], N)
        end
    end

    @show 
    SzSz_diff = 2.0^(-L)* μed - μcc
    @show abs.(SzSz_diff) |> maximum
end
    
function check_Sz_fourier_correlation(H, ifn, L)
    X,Y,Z,P,M = pauli_matrices_sparse(L)
    Z = map(full, Z)

    #julia's h5read is broken. Need the reshape to fix.
    for q = 0:Int(L/2)-1
        μcc  =    readdlm("$(ifn).SzqSzq.$(q).re")
        μcc += im*readdlm("$(ifn).SzqSzq.$(q).im")[:,1:end-1]
        N = size(μcc, 2)

        Zq = zeros(2^L, 2^L)
        for j = 1:L
            Zq += exp(im*j*(2*pi*q/L)) * Z[j]
        end
        μed = correlation(H, Zq, Zq', N)/4
        SzqSzq_diff = 2.0^(-L)* μed - μcc
        @show q, abs.(SzqSzq_diff) |> maximum
    end
end

s = ArgParseSettings()
@add_arg_table s begin
    "-i"
        help     = "input filename base"
        arg_type = String 
    "-m"
        help     = "model"
        arg_type = String
    "-o"
        help     = "output filename base"
        arg_type = String
end

function main(args)
    O = parse_args(args, s) |> Dict
    ifn = O["i"]
    ofn = O["o"]
    model = O["m"]

    hzf = open("$ifn.dis")
    hz = [parse(Float64, s) for s in split(readline(hzf))]
    L = length(hz)
    close(hzf)
    if L >= 14
        println("L = $L >= 14: skipping ED check")
    else
        
        X,Y,Z,P,M = pauli_matrices_sparse(L)
        if "rfheis" == model
            H = rfheis(hz)
            j = imag(sum([0.5*im*(P[l]*M[l+1] - M[l]*P[l+1]) for l in 1:L-1]) |> full)
        elseif "2NJW" == model
            H = rf_2NJW(hz)
            j  = imag(sum([im*(P[l]*M[l+1] - M[l]*P[l+1]) for l in 1:L-1]) |> full)
            j += imag(sum([im*(P[l]*Z[l+1]*M[l+2] - M[l]*Z[l+1]*P[l+2]) for l in 1:L-2]) |> full)
        else
            error("Bad model $model")
        end
        #@show trace(H^2)
        check_dos_trace(H, ifn, ofn, L)
        check_conductivity(H, j, ifn, ofn, L)
        check_Sz_fourier_correlation(H,ifn,L)
        check_Sz_correlation(H, ifn, L)
    end
end
main(ARGS)
