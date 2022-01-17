module QuadTrees

    import Base.∈, Base.∉, Plots.Shape, Base.insert!, Base.-, Base.+

    export QuadTree, Point, insert!, ∈, ∉, Draw, Draw!, Query,
        CircularParticle, AxisAlignedBox, CollisionEvent, collisions

    abstract type Object end

    const Point = Tuple{Float64,Float64} # cartesian point in 2D

    Point(x::Float64,y::Float64) = Point((x,y))

    +(p::Point,q::Point) = Point(p[1]+q[1],p[2]+q[2])
    -(p::Point,q::Point) = Point(p[1]-q[1],p[2]-q[2])

    struct CircularParticle <: Object
        pos::Point
        r::Float64
        id::UInt128
    end

    struct AxisAlignedBox
        # box at x, with half-size-length h
        x::Point
        h::Float64
    end

    width(a::AxisAlignedBox)::Float64 = a.h*2.0

    mutable struct QuadTree
        """
            A QuadTree defined by a node structure
        """
        root::Bool
        boundary::AxisAlignedBox
        objects::Vector{Object}
        minWidth::Float64
        maxObjects::Int64
        nw::Union{Nothing,QuadTree}
        ne::Union{Nothing,QuadTree}
        sw::Union{Nothing,QuadTree}
        se::Union{Nothing,QuadTree}
    end

    width(q::QuadTree)::Float64 = width(q.boundary)

    function ∈(p::Point,b::AxisAlignedBox)::Bool
        """
            Test if a point is inside the given box (inc perimeter)
        """
        if (p[1] >= b.x[1]-b.h && p[1] <= b.x[1]+b.h
            && p[2] >= b.x[2]-b.h && p[2] <= b.x[2]+b.h)
            return true
        else
            return false
        end
    end

    function ∈(o::CircularParticle,b::AxisAlignedBox)::Bool
        """
            Test if CircularParticle is inside the given box,
                accounts for radius. Must be strictly inside.
        """
        r = o.r
        p = o.pos
        if (p[1]-r > b.x[1]-b.h && p[1]+r < b.x[1]+b.h
            && p[2]-r > b.x[2]-b.h && p[2]+r < b.x[2]+b.h)
            return true
        else
            return false
        end
    end


    ∉(p::Point,b::AxisAlignedBox)::Bool = !∈(p,b)

    ∈(o::Object,b::AxisAlignedBox)::Bool = ∈(o,b)
    ∉(o::Object,b::AxisAlignedBox)::Bool = !∈(o,b)

    function points(a::AxisAlignedBox)::Vector{Point}
        """
            Short hand for the corners of an AxisAlignedBox
        """
        return  Point.([(a.x[1]-a.h,a.x[2]-a.h),
                (a.x[1]-a.h,a.x[2]+a.h),
                (a.x[1]+a.h,a.x[2]-a.h),
                (a.x[1]+a.h,a.x[2]+a.h)
            ])
    end

    function intersects(a::AxisAlignedBox,b::AxisAlignedBox)::Bool
        """
            Tests if any of the corners of box a are in box b
        """
        for p in points(a)
            if p ∈ b
                return true
            end
        end
        return false
    end

    function subdivide(b::AxisAlignedBox)::Vector{AxisAlignedBox}
        """
            Breaks a box into 4 equal size boxes and returns these as
            new instances
        """
        h = b.h/2.0
        return  [
                    AxisAlignedBox(Point(b.x[1]-h,b.x[2]+h),h),
                    AxisAlignedBox(Point(b.x[1]+h,b.x[2]+h),h),
                    AxisAlignedBox(Point(b.x[1]-h,b.x[2]-h),h),
                    AxisAlignedBox(Point(b.x[1]+h,b.x[2]-h),h),
                ]
    end

    ∈(p::Point,q::QuadTree)::Bool = ∈(p,q.boundary)
    ∉(p::Point,q::QuadTree)::Bool = ∉(p,q.boundary)
    ∈(o::Object,q::QuadTree)::Bool = ∈(o,q.boundary)
    ∉(o::Object,q::QuadTree)::Bool = !∈(o,q.boundary)
    isempty(q::QuadTree)::Bool = q.objects|>length == 0
    full(q::QuadTree)::Bool = length(q.objects) >= q.maxObjects

    # useful constructors
    QuadTree(a::AxisAlignedBox)::QuadTree = QuadTree(false,a,Vector{CircularParticle}(),0.1,100,nothing,nothing,nothing,nothing)
    QuadTree(r::Bool,a::AxisAlignedBox)::QuadTree = QuadTree(r,a,Vector{CircularParticle}(),0.1,100,nothing,nothing,nothing,nothing)
    QuadTree(
        x::Float64,
        y::Float64,
        h::Float64)::QuadTree = QuadTree(AxisAlignedBox(Point(x,y),h))
    QuadTree(
        r::Bool,
        x::Float64,
        y::Float64,
        h::Float64)::QuadTree = QuadTree(r,AxisAlignedBox(Point(x,y),h))


    function subdivide!(q::QuadTree)::Nothing
        """
            Calls subdivide on this nodes boundary box, and
            replaces the four segments with new boxes
        """
        boxes = subdivide(q.boundary)
        q.nw = QuadTree(boxes[1])
        q.ne = QuadTree(boxes[2])
        q.sw = QuadTree(boxes[3])
        q.se = QuadTree(boxes[4])

        for c in [q.nw,q.ne,q.sw,q.se]
            c.minWidth = q.minWidth
            c.maxObjects = q.maxObjects
        end
        nothing
    end

    function insert!(o::Object,q::QuadTree)::Bool
        """
            Attempt to insert the point p into the quadtree
        """
        if (o ∉ q)
            return false # outside bounds
        end

        if width(q) > q.minWidth && q.nw == nothing
            # subdivide
            subdivide!(q)
        end

        # attempt to go lower
        for c in [q.nw,q.ne,q.sw,q.se]
            if (c != nothing && insert!(o,c))
                return true
            end
        end

        # this is the lowest we can go, place here
        if !full(q)
            push!(q.objects,o)
            return true
        end

        return false # could not place (should not happen)
    end

    function query(a::AxisAlignedBox, Q::QuadTree)
        """
            Return a list of points (in the QuadTree Q) which are also
            inside the box a
        """
        if ~intersects(a,Q.boundary)
            return []
        end

        if Q.point == nothing
            return []
        end

        result = []

        for o in Q.objects
            if (o ∈ a)
                push!(result,o)
            end
        end

        for q in [Q.nw,Q.ne,Q.sw,Q.se]
            if q != nothing
                for p in query(a,q)
                    push!(result,p)
                end
            end
        end

        return result
    end

    function size(Q::QuadTree)
        s = 1
        for q in [Q.nw,Q.ne,Q.sw,Q.se]
            if q != nothing
                s += size(q)
            end
        end
        return s
    end

    # a pair of indices, i, and j, and a distance d
    struct CollisionEvent
        i::UInt128
        j::UInt128
        d::Float64
    end

    """
        For a given object find all collisions with it
    """
    function collisions(q::QuadTree,o::Object)
        col = Vector{CollisionEvent}()
        for j in 1:length(q.objects)
            o.id == q.objects[j].id ? continue : nothing
            rc = o.r+q.objects[j].r
            r = q.objects[j].pos-o.pos
            d2 = r[1]*r[1]+r[2]*r[2]
            if (d2 < rc*rc)
                push!(col,CollisionEvent(
                    o.id,
                    q.objects[j].id,
                    d2
                ))
            end
        end
        # down the children
        for c in [q.nw,q.ne,q.sw,q.se]
            if (c != nothing)
                for ccol in collisions(c,o)
                    push!(col,ccol)
                end
            end
        end

        return col
    end

    """
        Collect all collision events in the tree
    """
    function collisions(q::QuadTree)
        col = Vector{CollisionEvent}()
        for i in 1:length(q.objects)
            for j in i+1:length(q.objects)
                rc = q.objects[i].r+q.objects[j].r
                r = q.objects[j].pos-q.objects[i].pos
                d2 = r[1]*r[1]+r[2]*r[2]
                if (d2 < rc*rc)
                    push!(col,CollisionEvent(
                        q.objects[i].id,
                        q.objects[j].id,
                        d2
                    ))
                end
            end
        end

        # go down the children
        for c in [q.nw,q.ne,q.sw,q.se]
            if c == nothing
                continue
            end
            for i in 1:length(q.objects)
                for ccol in collisions(c,q.objects[i])
                    push!(col,ccol)
                end
            end
            for ccol in collisions(c)
                push!(col,ccol)
            end
        end

        return col
    end
end
