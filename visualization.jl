using Plots
using DataFrames
using CSV

frameRate = 1/24
global timeInSec=0.0-frameRate
global secsToPlot = 3

function getDataFromFile()
    global timeInSec += frameRate
    df = CSV.read("data/"*string(timeInSec)*".csv")
    stats = convert(Matrix{Float64}, df)
    numOfBodies = size(stats)[1]
    xs = [stats[1,2]]
    ys = [stats[1,3]]
    zs = [stats[1,4]]
    ms = [stats[1,1]]
    for i in 2:numOfBodies
        push!(xs,stats[i,2])
        push!(ys,stats[i,3])
        push!(zs,stats[i,4])
        push!(ms,stats[i,1])
    end
    rm("data/"*string(timeInSec)*".csv")
    return xs, ys, zs, ms
end

function getMarkerSize(ms::Array{Float64,1})
    size = []
    for i in 1:length(ms)
        dot = log(ms[i]/M)
        push!(size, dot)
    end
    return size
end
function plotting()
    # initialize a 3D plot with 1 empty series
    global timeInSec = nothing
    global timeInSec = 0.0-frameRate
    xs, ys, zs, ms = getDataFromFile()
    len = length(xs)
    minX = minimum(xs)
    minY = minimum(ys)
    minZ = minimum(zs)
    maxX = maximum(xs)
    maxY = maximum(ys)
    maxZ = maximum(zs)
    plt = Plots.plot([xs[1]],[ys[1]], [zs[1]],
        seriestype=:scatter, marker= (:yellow, stroke(3, 0.2, :yellow)), markersize = 0, camera = (30, 30),background_color = :black,
        xlim = (minX - 100, maxX + 100),
        ylim = (minY - 100, maxY + 100),
        zlim = (minZ - 100, maxZ + 100)
    )
    size = getMarkerSize(ms[2:len])
    scatter!(xs[2:len],ys[2:len], zs[2:len],marker= (:white, stroke(3, 0.2, :white)), markersize = size, camera = (30,30),background_color = :black)

    display(plt)


    # build an animated gif by pushing new points to the plot, saving every 10th frame
    @gif for i=2:secsToPlot*24+1
        xs, ys, zs, ms = getDataFromFile()
        len = length(xs)

        plot([xs[1]],[ys[1]], [zs[1]],
            seriestype=:scatter, marker= (:yellow, stroke(3, 0.2, :yellow)), markersize = 0, camera = (30,30),background_color = :black,
            xlim = (minX - 1e2, maxX + 1e2),
            ylim = (minY - 1e2, maxY + 1e2),
            zlim = (minZ - 1e2, maxZ + 1e2)
        )

        size = getMarkerSize(ms[2:len])
        scatter!(xs[2:len],ys[2:len], zs[2:len],marker= (:white, stroke(3, 0.2, :white)), markersize  = size, camera = (30,30),background_color = :black)
    end

end
