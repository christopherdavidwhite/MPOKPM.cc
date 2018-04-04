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

PyDict(pyimport("matplotlib")["rcParams"])["xtick.labelsize"] = 20
PyDict(pyimport("matplotlib")["rcParams"])["ytick.labelsize"] = 20
PyPlot.rc("text", usetex=true)

s = ArgParseSettings()
@add_arg_table s begin
    "-i"
        help     = "input filename base"
        arg_type = String 
    "-o"
        help     = "output filename base"
        arg_type = String
end

function main(args)
    O = parse_args(args, s) |> Dict
    ifn = O["i"]
    ofn = O["o"]


    #plot bond dim. as fn of chebyshev order
    nM = readdlm("$ifn.chM")
    @show nM
    semilogy(nM[:,1], nM[:,2], ".")
    xlabel("Chebyshev order \$n\$", size=20)
    ylabel("bond dimension M", size=20)
    savefig("$ofn-plt-chM.pdf", bbox_inches="tight")

    #plot sing. vals. of every fourth Chebyshev

    cla()
    clf()
    
    f = open("$ifn.chs")
    for (n, ln) in enumerate(eachline(f))
        if n % 4 == 1
            semilogy([parse(Float64, s) for s in split(ln)], ".-", label="n=$n")
        end
    end
    close(f)
    legend()
    ylabel("Singular value \$s_j\$", size=20)
    xlabel("Index \$j\$", size=20)
    xlim(0,520)
    
    savefig("$ofn-plt-chs.pdf", bbox_inches="tight")
    cla()
    clf()

    hzf = open("$ifn.dis")
    hz = [parse(Float64, s) for s in split(readline(hzf))]
    close(hzf)
    L = length(hz)
    if L > 14
        println("L = $L > 14: skipping ED check")
    else
        tnf = open("$ifn.re")
        trTncc = 2^(L/2)*[parse(Float64, s) for s in split(readline(tnf))]
        N = length(trTncc)
        close(tnf)

        L = length(hz)
        X,Y,Z,P,M = pauli_matrices_sparse(L)
        H  = 0.25*sum(map((x1,x2) -> x1 * x2, X[1:end-1], X[2:end]))
        H += 0.25*sum(map((x1,x2) -> x1 * x2, Y[1:end-1], Y[2:end]))
        H += 0.25*sum(map((x1,x2) -> x1 * x2, Z[1:end-1], Z[2:end]))
        H += 0.5*sum(hz.*Z)
        H = H/(3*(L-1)*0.25 + sum(abs.(hz))*0.5)
        #@show trace(H^2)
        d,v = H |> full |> eig

        trTned = zeros(N)
        for n = 0:N-1
            trTned[n+1] = sum(cos.(n*acos.(d)))
        end
        semilogy(abs.(trTned - trTncc), ".")
        xlabel("n", size=20)
        ylabel("err in trace Tn",size=20)
        savefig("$ofn-plt-trTnerr.pdf", bbox_inches="tight")
        cla()
        clf()
    end
end
main(ARGS)
