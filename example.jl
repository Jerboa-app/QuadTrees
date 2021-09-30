include("QuadTrees/QuadTrees.jl")

using Main.QuadTrees, Plots

# construct an empty tree from at the region (-0.5,0), (0.5,0), (0,0.5), (0,-0.5)
Q = QuadTrees.QuadTree(.5,.5,0.5)
# random points in a square: (0,1), (1,1), (1,0), (0,0)
p = rand(2^9,2)

for i in 1:size(p,1)
    insert!(Point(p[i,:]...),Q)
end

scatter(p[:,1],p[:,2],label="",dpi=300,
    aspect_ratio=:equal,xaxis=nothing,yaxis=nothing,axis=:off
)

Draw!(Q)
plot!()

savefig("quadtree.png")
