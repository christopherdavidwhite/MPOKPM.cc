###print

using ArgParse
using PyPlot
using PyCall
@pyimport seaborn as sns
sns.set_style("whitegrid")

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
    plot(nM[:,1], nM[:,2], ".")
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
end

main(ARGS)
