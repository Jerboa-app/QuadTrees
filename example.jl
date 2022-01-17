include("QuadTrees.jl")

using Main.QuadTrees, GLMakie, ProgressMeter

function Draw!(a::AxisAlignedBox,ax)
    linesegments!(ax,[
        Point2(a.x[1]-a.h,a.x[2]-a.h),
        Point2(a.x[1]-a.h,a.x[2]+a.h),
        Point2(a.x[1]-a.h,a.x[2]+a.h),
        Point2(a.x[1]+a.h,a.x[2]+a.h),
        Point2(a.x[1]+a.h,a.x[2]+a.h),
        Point2(a.x[1]+a.h,a.x[2]-a.h),
        Point2(a.x[1]+a.h,a.x[2]-a.h),
        Point2(a.x[1]-a.h,a.x[2]-a.h)
    ])
end

function Draw!(Q::QuadTree,ax)
    Draw!(Q.boundary,ax)
    for q in [Q.nw,Q.ne,Q.sw,Q.se]
        q != nothing ? Draw!(q,ax) : nothing
    end
end

#construct an empty tree from at the region (-0.5,0), (0.5,0), (0,0.5), (0,-0.5)
N = 256
Q = QuadTrees.QuadTree(.5,.5,.5)
Q.minWidth=0.025
Q.maxObjects=10
# random points in a square: (0,1), (1,1), (1,0), (0,0)
p = rand(N,2)

for i in 1:size(p,1)
    insert!(CircularParticle(QuadTrees.Point(p[i,:]...),0.025,i),Q)
end

fig = Figure()
ax = Axis(fig[1,1])
ax.aspect=AxisAspect(1)
scatter!(ax,p[:,1],p[:,2],markersize=(0.025,0.025),markerspace=SceneSpace)

Draw!(Q,ax)

save("quadtree.png",fig)

Q2 = QuadTrees.QuadTree(.5,.5,0.5)
Q2.minWidth=0.01
Q2.maxObjects=100
p = [CircularParticle(QuadTrees.Point(rand(2)...),0.001,i) for i in 1:2^12]
@info 1.0/(@timed [insert!(n,Q2) for n in p])[2]

p = [CircularParticle(QuadTrees.Point(2.0.+(18.0-2.0) .*rand(2)...),1.0,i) for i in 1:64]

pos = Vector{Point2}()
angles = rand(length(p)).*2.0*π

for i in 1:size(p,1)
    push!(pos,Point2(p[i].pos...))
end

pos = Node(pos)

fig = Figure()
ax = Axis(fig[1,1])

ax.aspect=AxisAspect(1)
limits!(ax,0.,20.,0.,20.)

scatter!(ax,pos,markersize=(2,2),markerspace=SceneSpace)

T = 10000
dt = 0.001
prog = Progress(T)

record(fig, "example.mp4",collect(1:T),framerate=60) do i

    Q = QuadTrees.QuadTree(20/2.,20/2.,20.)
    Q.minWidth=1.0
    Q.maxObjects=100
    @assert sum([insert!(n,Q) for n in p])==64

    C = collisions(Q)

    F = zeros(length(p),2)

    for c in C
        r = (p[c.j].r+p[c.i].r)
        rx = p[c.j].pos[1] - p[c.i].pos[1]
        ry = p[c.j].pos[2] - p[c.i].pos[2]
        d = sqrt(c.d)
        k = -100. * (r-d)
        F[c.i,:] .+= [k*rx/d,k*ry/d]
        F[c.j,:] .-= [k*rx/d,k*ry/d]
    end

    newpos = []
    for n in 1:length(p)
        angles[n] = angles[n]+sqrt(2.0*0.005*dt)*randn()
        x = p[n].pos[1] + dt*(1.0*cos(angles[n])+F[n,1])
        y = p[n].pos[2] + dt*(1.0*sin(angles[n])+F[n,2])

        r = p[n].r

        vx = x-p[n].pos[1]
        vy = y-p[n].pos[2]
        ux = 0.
        uy = 0.
        flag = false
        ang = 0.

        if (x < r) || (x > 20.0-r)
            ux = -1.0*vx
            flag = true
        end

        if (y<r)||(y > 20.0-r)
            uy = -1.0*vy
            flag = true
        end

        if (flag)
            x += ux
            y += uy
        end

        push!(newpos,Point2(x,y))

        p[n] = CircularParticle(QuadTrees.Point(x,y),p[n].r,p[n].id)
    end

    pos[] = newpos

end
