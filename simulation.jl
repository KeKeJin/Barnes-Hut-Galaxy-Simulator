includet("barnesHut.jl")
includet("visualization.jl")

# this is a helper function called by other specific simulator functions
function basicSimulation(C::Node, num::Int64)
    global timesInSec = 0.00


    creep(C)
    # keep a local dataframe to update
    df = writeStatsToFrame(C)
    # write this stats to a csv file
    CSV.write("data/"*string(timesInSec)*".csv", df)


    secs = 30
    frames = secs*60
    # simulate new positions for "frames" number of frames
    for i in 1:frames
        C = makeTreeFromFrame(df)
        creep(C)
        df = writeStatsToFrame(C)
        global timesInSec += frameRate
        CSV.write("data/"*string(timesInSec)*".csv", df)
    end
    plotting()
end

# this function is a helper function that returns minBounds and maxBounds based on the type
function getBoundsBasedOnTypes(type::String, vector::Array{Float64,1}, minX::Float64, maxX::Float64)
    minBounds = SVector{3, Float64}(zeros(Float64,3))
    maxBounds = SVector{3, Float64}(zeros(Float64,3))
    if type == "random"
        minBounds = SVector{3, Float64}(minX,minX,minX)
        maxBounds = SVector{3, Float64}(maxX,maxX,maxX)
    elseif type == "disk"
        normal = vector
        # calcuate the possible max and min value of z for initializing a tree
        if normal[3] == 0 && normal[1] != 0 && normal[2] != 0
            minBounds = SVector{3, Float64}(-maxX,-maxX,minX)
            maxBounds = SVector{3, Float64}(maxX,maxX,maxX)
        else
            println(-minX*normal[1]-minX*normal[2],-maxX*normal[1]-maxX*normal[2])
            minZ = min(-minX*normal[1]-minX*normal[2],-maxX*normal[1]-maxX*normal[2])/normal[3]
            maxZ = max(-minX*normal[1]-minX*normal[2],-maxX*normal[1]-maxX*normal[2])/normal[3]
            if minZ == maxZ
                minZ = -maxX
                maxZ = maxX
            end
            minBounds = SVector{3, Float64}(minX,minX,minZ)
            maxBounds = SVector{3, Float64}(maxX,maxX,maxZ)
        end
    elseif type == "line"
        # calcuate the possible max and min value of y and z for initializing a tree
        dir = vector
        minY = min(minX/dir[1] , maxX/dir[1])*dir[2]
        maxY = max(minX/dir[1] , maxX/dir[1])*dir[2]
        minZ = min(minX/dir[1] , maxX/dir[1])*dir[3]
        maxZ = max(minX/dir[1] , maxX/dir[1])*dir[3]
        minBounds = SVector{3, Float64}(minX,minY,minZ)
        maxBounds = SVector{3, Float64}(maxX,maxY,maxZ)
    end
    return minBounds, maxBounds
end

# this function is a helper function that adds bodies based on different system
function addBodies(C::Node, type::String, num::Int64, vector::Array{Float64,1}, minn::Int64,maxx::Int64)
    if type == "random"
        # add random bodies
        for i in 1:num
            mass = exp(randn())*M
            position = [rand(minn:maxx)+rand(Float64), rand(minn:maxx)+rand(Float64) ,rand(minn:maxx)+rand(Float64)]
            newBody = Body(i, mass, position, SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
            insertBody(C,newBody)
        end
    elseif type == "disk"
        normal = vector
        # cases for different normal vector
        if normal[3]==0 && normal[1]!= 0 && normal[2]!= 0
            # add bodies in a disk shaped space
            for i in 1:num
                mass = exp(randn())*M
                x = rand(minn:maxx)+rand(Float64)
                y = (0-normal[1]*x)/normal[2]
                z = rand(minn:maxx)+rand(Float64)
                newBody = Body(i, mass, [x,y,z], SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
                insertBody(C,newBody)
            end
        elseif normal[3]==0 && normal[1]== 0 && normal[2]!= 0
            # add bodies in a disk shaped space
            for i in 1:num
                mass = exp(randn())*M
                x = rand(minn:maxx)+rand(Float64)
                y = 0
                z = rand(minn:maxx)+rand(Float64)
                newBody = Body(i, mass, [x,y,z], SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
                insertBody(C,newBody)
            end
        elseif normal[2]==0 && normal[1]!= 0 && normal[3]!= 0
            # add bodies in a disk shaped space
            for i in 1:num
                mass = exp(randn())*M
                x = rand(minn:maxx)+rand(Float64)
                y = rand(minn:maxx)+rand(Float64)
                z = (0-normal[1]*x)/normal[3]
                newBody = Body(i, mass, [x,y,z], SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
                insertBody(C,newBody)
            end
        elseif normal[2]==0 && normal[1]== 0 && normal[3]!= 0
            # add bodies in a disk shaped space
            for i in 1:num
                mass = exp(randn())*M
                x = rand(minn:maxx)+rand(Float64)
                y = rand(minn:maxx)+rand(Float64)
                z = 0
                newBody = Body(i, mass, [x,y,z], SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
                insertBody(C,newBody)
            end
        elseif normal[1]==0 && normal[2]!= 0 && normal[3]!= 0
            # add bodies in a disk shaped space
            for i in 1:num
                mass = exp(randn())*M
                x = rand(minn:maxx)+rand(Float64)
                y = rand(minn:maxx)+rand(Float64)
                z = (0-normal[2]*y)/normal[3]
                newBody = Body(i, mass, [x,y,z], SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
                insertBody(C,newBody)
            end
        else
            for i in 1:num
                mass = exp(randn())*M
                x = rand(minn:maxx)+rand(Float64)
                y = rand(minn:maxx)+rand(Float64)
                z = (0-x*normal[1]-y*normal[2])/normal[3]
                newBody = Body(i, mass, [x,y,z], SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
                insertBody(C,newBody)
            end
        end
    elseif type == "line"
        dir = vector
        if dir[1]!=0
            for i in 1:num
                mass = exp(randn())*M
                x = rand(minn:maxx)+rand(Float64)
                y = dir[2]*x/dir[1]
                z = dir[3]*x/dir[1]
                newBody = Body(i, mass, [x,y,z], SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
                insertBody(C,newBody)
            end
        elseif dir[2]!=0
            for i in 1:num
                mass = exp(randn())*M
                y = rand(minn:maxx)+rand(Float64)
                x = dir[1]*y/dir[2]
                z = dir[3]*y/dir[2]
                newBody = Body(i, mass, [x,y,z], SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
                insertBody(C,newBody)
            end
        end
    end
end

# this is a random galaxy simulator similating a space of 1000 pc * 1000 pc *1000 pc
# optional argument: num -> number of bodies in the system
function randomGalaxy(num::Int64=100)
    minBounds, maxBounds = getBoundsBasedOnTypes("random",[0.0,0.0,0.0], -1.0, 1001.0)
    # making an empty tree
    C= Node(
           false, Array{Node,1}(undef, 8),0,
           Body(), 0.0, SVector{3, Float64}(zeros(3)),
           minBounds,maxBounds
           )
    # add random bodies
    addBodies(C,"random",num, [0.0,0.0,0.0],1,1000)
    # assign initial velocity for all the body
    caluculateInitialSpeed(C)

    # tuning for different cases
    if num== 2
        global timeScalar = nothing
        global timeScalar = 60*60*60*24*365*7e7
        global strengthOfInteraction = nothing
        global strengthOfInteraction = 3
    else
        global timeScalar = nothing
        global timeScalar = 60*60*60*24*365*5e6
        global strengthOfInteraction = nothing
        global strengthOfInteraction = 1.25
    end
    basicSimulation(C,num)

end

# this is a disk-shaped galaxy simulator similating a space of 1000 pc * 1000 pc *1000 pc
# optional arguments: normal-> the normal vector of the plane;
#       num -> number of bodies in the system
function diskGalaxy(normal::Array{Float64,1} = [1.0,1.0,0.0], num::Int64=100)
    if length(normal)!=3
        error("Normal vector must have 3 components")
    end

    # invert the sign to positive
    if normal[1] < 0 && normal[2] <0 && normal[3] <0
        normal = -normal
    end


    minBounds, maxBounds = getBoundsBasedOnTypes("disk",normal, -1.0, 1001.0)
    # making an empty tree
    C= Node(
           false, Array{Node,1}(undef, 8),0,
           Body(), 0.0, SVector{3, Float64}(zeros(3)),
           minBounds,maxBounds
           )
    # add bodies in a plane
    addBodies(C,"disk",num,normal,1,1000)
    # assign initial velocity for all the body
    caluculateInitialSpeed(C)
    if num== 2
        global timeScalar = nothing
        global timeScalar = 60*60*60*24*365*7e7
        global strengthOfInteraction = nothing
        global strengthOfInteraction = 3
    else
        global timeScalar = nothing
        global timeScalar = 60*60*60*24*365*5e6
        global strengthOfInteraction = nothing
        global strengthOfInteraction = 1
    end
    basicSimulation(C,num)
end

# this is a line-shaped galaxy simulator similating a space of 1000 pc * 1000 pc *1000 pc
# optional arguments: dir-> the direction of the line;
#                       num -> number of bodies on the line
function lineGalaxy(dir::Array{Float64,1}=[1.0,1.0,0.0], num::Int64=100)
    if length(dir)!=3
        error("Direction vector must have 3 components")
    end

    # invert the sign to positive
    if dir[1] < 0 && dir[2] <0 && dir[3] <0
        dir = -dir
    end

    minBounds, maxBounds = getBoundsBasedOnTypes("line",dir, -1.0, 1001.0)
    # making an empty tree
    C= Node(
           false, Array{Node,1}(undef, 8),0,
           Body(), 0.0, SVector{3, Float64}(zeros(3)),
           minBounds,maxBounds
           )
    # add bodies in a line
    addBodies(C,"line",num,dir,1,1000)
    # assign initial velocity for all the body
    caluculateInitialSpeed(C)
    if num== 2
        global timeScalar = nothing
        global timeScalar = 60*60*60*24*365*7e7
        global strengthOfInteraction = nothing
        global strengthOfInteraction = 3
    else
        global timeScalar = nothing
        global timeScalar = 60*60*60*24*365*5e5
        global strengthOfInteraction = nothing
        global strengthOfInteraction = 0.1
    end
    basicSimulation(C,num)
end

# this is a system of two galaxies simulator similating a space of 1000 pc * 1000 pc *1000 pc
# optional arguments: num1 -> number of bodies in the first galaxy,
#                     num2 -> number of bodies in the second galaxy,
#                     velocity1 -> the velocity of the first galaxy,
#                     velocity2 -> the velocity of the second galaxy,
#                     type1 -> the type of the first galaxy,
#                     type2 -> the type of the second galaxy,
#                     vector1 -> the normal vector or direction vector of the first galaxy,
#                     vector2 -> the normal vector or direction vector of the second galaxy]
function twoCollapseGalaxy(num1::Int64=100,num2::Int64=100, velocity1::SVector{3,Float64}=SVector{3,Float64}(50.0,50.0,50.0),
    velocity2::SVector{3,Float64}=SVector{3,Float64}(-100.0,-100.0,-100.0), type1::String="disk",type2::String="disk",
    vector1::Array{Float64,1}=[1.0,1.0,0.0], vector2::Array{Float64,1}=[1.0,1.0,0.0])
    if length(velocity2)!=3 || length(velocity1)!=3
        error("Velocity vector(s) must have 3 components")
    elseif length(vector1) != 3 || length(vector2) != 3
        error("direcion/normal vector(s) must have 3 components")
    elseif type1 != "disk" && type1 != "random" && type1 != "line"
        error("Type1 not supported")
    elseif type2 != "disk" && type2 != "random" && type2 != "line"
        error("Type2 not supported")
    end


    #initialize a new empty tree
    minBounds1, maxBounds1 = getBoundsBasedOnTypes(type1,vector1, -1.0, 500.0)
    minBounds2, maxBounds2 = getBoundsBasedOnTypes(type2,vector2, 500.0, 1001.0)

    netMinBounds = SVector{3, Float64}(min(minBounds1[1],minBounds2[1]),
                                       min(minBounds1[2],minBounds2[2]),
                                       min(minBounds1[3],minBounds2[3]))

    netMaxBounds = SVector{3, Float64}(max(maxBounds1[1],maxBounds2[1]),
                                       max(maxBounds1[2],maxBounds2[2]),
                                       max(maxBounds1[3],maxBounds2[3]))
    C= Node(
           false, Array{Node,1}(undef, 8),0,
           Body(), 0.0, SVector{3, Float64}(zeros(3)),
           netMinBounds, netMaxBounds
           )


    if type1 == "disk"
        addBodies(C,"disk",num1, vector1,1,500)
    elseif type1 == "random"
        addBodies(C,"random",num1, [0.0,0.0,0.0],1,500)
    elseif type1 == "line"
        addBodies(C,"line",num1, vector1,1,500)
    end
    println(C.numOfChild)
    caluculateInitialSpeedPlus(C,velocity1)

    if type2 == "disk"
        println("here")
        addBodies(C,"disk",num2,vector2,500,1000)
    elseif type2 == "random"
        addBodies(C,"random",num2, [0.0,0.0,0.0],500,1000)
    elseif type2 == "line"
        addBodies(C,"line",num1, vector2,500,1000)
    end

    caluculateInitialSpeedPlus(C,velocity2)

    println(C.numOfChild)
    global timeScalar = nothing
    global timeScalar = 60*60*60*24*365*5e6
    global strengthOfInteraction = nothing
    global strengthOfInteraction = 1
    basicSimulation(C,num1+num2)
end

# this function changes theta (s/d)
function changeTheta(theta::Float64)
    global THETA = nothing
    global THETA = theta
end
