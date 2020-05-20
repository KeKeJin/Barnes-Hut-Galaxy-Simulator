global timeInSec=0.0
frameRate = 1/60

function getDataFromFile()
    global timeInSec+=frameRate
    df = CSV.read("data/"*string(timeInSec)*".csv")
    stats = convert(Matrix{Float64}, df)
    numOfBodies = size(stats)[1]
    xs = [stats[1,2]]
    ys = [stats[1,3]]
    zs = [stats[1,4]]
    for i in 2:numOfBodies
        push!(xs,stats[i,2])
        push!(ys,stats[i,3])
        push!(zs,stats[i,4])
    end

    return xs, ys, zs
end

function plotting()
    # initialize a 3D plot with 1 empty series
    xs, ys, zs = getDataFromFile()
    plt = Plots.plot(xs,ys, zs,
         seriestype=:scatter, markersize = 2,
         xlim = (-50, 1050),
    ylim = (-50, 1050),
    zlim = (-50, 1050))
    display(plt)


    # # build an animated gif by pushing new points to the plot, saving every 10th frame
    @gif for i=1:1800-50
        xs, ys, zs = getDataFromFile()

        plot(xs,ys, zs,
             seriestype=:scatter, markersize = 2,xlim = (-50, 1050),
        ylim = (-50, 1050),
        zlim = (-50, 1050))
    end

end
