module QuadTrees

    using Plots

    import Base.∈, Base.∉, Plots.Shape

    export QuadTree, Point, insert!, ∈, ∉, Draw, Draw!, Query

    const Point = Tuple{Float64,Float64} # cartesian point in 2D

    Point(x::Float64,y::Float64) = Point((x,y))

    struct AxisAlignedBox
        # box at x, with half-size-length h
        x::Point
        h::Float64
    end

    mutable struct QuadTree
        """
            A QuadTree defined by a node structure
        """
        root::Bool
        boundary::AxisAlignedBox
        point::Union{Nothing,Point}
        nw::Union{Nothing,QuadTree}
        ne::Union{Nothing,QuadTree}
        sw::Union{Nothing,QuadTree}
        se::Union{Nothing,QuadTree}
    end

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

    ∉(p::Point,b::AxisAlignedBox)::Bool = !∈(p,b)

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
    isempty(q::QuadTree)::Bool = q.point == nothing

    # useful constructors
    QuadTree(a::AxisAlignedBox)::QuadTree = QuadTree(false,a,nothing,nothing,nothing,nothing,nothing)
    QuadTree(r::Bool,a::AxisAlignedBox)::QuadTree = QuadTree(r,a,nothing,nothing,nothing,nothing,nothing)
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
        nothing
    end

    function insert!(p::Point,q::QuadTree)::Bool
        """
            Attempt to insert the point p into the quadtree
        """
        if (p ∉ q)
            return false # outside bounds
        end

        if isempty(q) && q.nw == nothing
            # subdivide and place point in whichever child
            subdivide!(q)
            if p ∈ q.nw
                q.nw.point = p
            end

            if p ∈ q.ne
                q.ne.point = p
            end
            if p ∈ q.sw
                q.sw.point = p
            end
            if p ∈ q.se
                q.se.point = p
            end
            return true
        end

        if (q.nw == nothing)
            subdivide!(q)
        end

        if insert!(p,q.nw)
            return true
        end
        if insert!(p,q.ne)
            return true
        end
        if insert!(p,q.sw)
            return true
        end
        if insert!(p,q.se)
            return true
        end

        return false # could not place (should not happen)
    end

    Shape(a::AxisAlignedBox) = Shape(
        a.x[1] .+ [-a.h,a.h,a.h,-a.h],
        a.x[2] .+ [a.h,a.h,-a.h,-a.h]
    )

    function Draw!(a::AxisAlignedBox)
        plot!(Shape(a),label="",fillalpha=0.0)
    end

    function Draw(Q::QuadTree,points=true)
        p = plot(aspect_ratio=:equal,label="")
        Draw!(Q.boundary)
        for q in [Q.nw,Q.ne,Q.sw,Q.se]
            q != nothing ? Draw!(q,points) : nothing
        end
        return p
    end

    function Draw!(Q::QuadTree,points=true)
        Draw!(Q.boundary)
        for q in [Q.nw,Q.ne,Q.sw,Q.se]
            q != nothing ? Draw!(q) : nothing
        end
    end

    function query(a::AxisAlignedBox, Q::QuadTree)
        """
            Return a list of points (in the QuadTree Q) which are also
            inside the box a
        """
        if ~intersects(a,Q.boundary)
            return Vector{Point}([])
        end

        if Q.point == nothing
            return Vector{Point}([])
        end

        result = Vector{Point}()

        if (Q.point ∈ a)
            push!(result,Q.point)
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
end
