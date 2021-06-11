using Documenter
using BioCCP
 

makedocs(sitename="BioCCP",
    format = Documenter.HTML(),
    modules=[BioCCP],  
    pages=Any["BioCCP"=> "man/BioCCP.md"])

#=
deploydocs(
            repo = "github.com/kirstvh/BioCCP.jl.git",
        )

=#

 