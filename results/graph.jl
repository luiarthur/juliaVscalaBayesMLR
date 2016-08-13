using PlotlyJS, Rsvg, DataFrames

dat = readtable("results.dat")


lang = dat[:,1]
cols = dat[:,2]
secs = dat[:,3]

langColor = [l == "sc" ? :mediumseagreen : :crimson for l in lang]
scalaInds = lang .== "sc"
juliaInds = lang .== "jl"

sS = scatter(;x=cols[scalaInds],y=secs[scalaInds],name="Scala",
             marker=attr(color=:mediumseagreen,size=20))
sJ = scatter(;x=cols[juliaInds],y=secs[juliaInds],name="julia",
             marker=attr(color=:crimson,size=20))
l = Layout(width=670,height=470,
           margin=attr(r=20,l=50,b=50,t=20),
           xaxis=attr(title="Columns in X matrix"),
           yaxis=attr(title="Execution Time (seconds)",zeroline=false))
p = plot([sS,sJ],l)

savefig(p, "results.html")
run(`sed -i.bak '8,13!d' results.html`)
run(`sed -i.bak 's/newPlot('\''[a-z|0-9|-]*'\''/newPlot('\''jlvsc'\''/' results.html`)
run(`rm -f results.html.bak`)
#plot([Plot(b1), Plot(b2)])
#=
include("graph.jl")
=#
