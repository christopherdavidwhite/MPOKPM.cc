###print

using ArgParse
using PyPlot
using PyCall
@pyimport seaborn as sns
sns.set_style("whitegrid")


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
    @time d,v = eig(H)
    H = zeros((2,2))
    gc()
    
    jv1 = copy(v)
    A_mul_B!(jv1,j,v)
    jv = Ac_mul_B(v,jv1)
    
    jv1 = zeros(2,2)
    v = zeros(2,2)
    gc()
    println("start loop")
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
    @show trTned
    @show trTncc
    trTndiff = trTned - trTncc
    @show abs.(trTndiff) |> maximum
    semilogy(abs.(trTned - trTncc), ".")
    xlabel("n", size=20)
    ylabel("err in trace Tn",size=20)
    savefig("$ofn-plt-trTnerr.pdf", bbox_inches="tight")
    cla()
    clf()
end


function check_conductivity(H, j, ifn,ofn,L)
    μcc = readdlm("$ifn.re")
    N = size(μcc,1)
    μed = -exact_μ(H,j,N) #minus because I took imag part of j

    diff = 2.0^(-2L)*μed - μcc
    vmax = abs.(diff) |> maximum
    @show μcc
    @show 2.0^(-2L)*μed
    pcolor(diff, cmap="seismic", vmin=-vmax, vmax=vmax)
    savefig("$ofn-plt-trTnjTmjerr.pdf", bbox_inches="tight")
    cla()
    clf()
        
    #@show μed
    #@show μcc
    @show abs.(diff) |> maximum 
end 

PyDict(pyimport("matplotlib")["rcParams"])["xtick.labelsize"] = 20
PyDict(pyimport("matplotlib")["rcParams"])["ytick.labelsize"] = 20
PyPlot.rc("text", usetex=true)

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

    #plot bond dim. as fn of chebyshev order
    nM = readdlm("$ifn.chM")
    semilogy(nM[:,1], nM[:,2], ".")
    xlabel("Chebyshev order \$n\$", size=20)
    ylabel("bond dimension M", size=20)
    savefig("$ofn-plt-chM.pdf", bbox_inches="tight")

    #plot sing. vals. of every fourth Chebyshev


    # f = open("$ifn.chs")
    # for (n, ln) in enumerate(eachline(f))
    #     if n % 4 == 1
    #         semilogy([parse(Float64, s) for s in split(ln)], ".-", label="n=$n")
    #     end
    # end
    # close(f)
    # legend()
    # ylabel("Singular value \$s_j\$", size=20)
    # xlabel("Index \$j\$", size=20)
    # xlim(0,520)
    
    # savefig("$ofn-plt-chs.pdf", bbox_inches="tight")
    # cla()
    # clf()

    hzf = open("$ifn.dis")
    hz = [parse(Float64, s) for s in split(readline(hzf))]
    L = length(hz)
    H = rfheis(hz)
    close(hzf)
    L = length(hz)
    if L >= 14
        println("L = $L >= 14: skipping ED check")
    else
        #@show trace(H^2)

        X,Y,Z,P,M = pauli_matrices_sparse(L)
        if "rfheis" == model
            H = rfheis(hz)
            j = imag(sum([2*(Y[l]*X[l+1] - X[l]*Y[l+1]) for l in 1:L-1]) |> full)
        elseif "2NJW" == model
            H = rf_2NJW(hz)
            j =  imag(sum([2*(Y[l]*X[l+1] - X[l]*Y[l+1]) for l in 1:L-1]) |> full) 
            j += imag(sum([1*(Y[l]*Z[l+1]X[l+2] - X[l]*Z[l+1]*Y[l+2]) for l in 1:L-2]) |> full) 
        else
            error("Bad model $model")
        end
        check_dos_trace(H, ifn, ofn, L)
        check_conductivity(H, j, ifn, ofn, L)
        #check conductivity coeffs
    end
end
main(ARGS)
