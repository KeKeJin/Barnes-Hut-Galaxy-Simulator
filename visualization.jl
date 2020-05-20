global timeInSec=75.13333333333036
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
         seriestype=:scatter, markersize = 2)
    display(plt)


    # # build an animated gif by pushing new points to the plot, saving every 10th frame
    @gif for i=1:850
        xs, ys, zs = getDataFromFile()

        plot(xs,ys, zs,
             seriestype=:scatter, markersize = 2)
    end every 10

end

# function makie()
#     N = 10
#     r = [(rand(7, 2) .- 0.5) .* 25 for i = 1:N]
#     xs, ys, zs = getDataFromFile()
#     scene = scatter(xs, ys, zs, markersize = 1)
#     record(scene, "output.mp4", r) do m
#         for i = 1:800
#             func!(scene)     # animate scene
#             recordframe!(io) # record a new frame
#         end
#     end
# end
