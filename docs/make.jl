using Documenter
using BioCCP
using STMOZOO.LocalSearch

makedocs(sitename="BioCCP",
    format = Documenter.HTML(),
    modules=[Example,
            SingleCellNMF,], # add your module
    pages=Any[
        "Example"=> "man/example.md",  # add the page to your documentation
        "GenProgAlign" => "man/gen_prog_align.md",
	       
    ])

#=
deploydocs(
            repo = "github.com/kirstvh/BioCCP.jl.git",
        )
=#
