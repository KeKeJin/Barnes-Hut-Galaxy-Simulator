includet("barnesHut.jl")
includet("visualization.jl")
global secs = 3

# this is a helper function called by other specific simulator functions
function basicSimulation(C::Node)
    global timesInSec = 0.00

    # keep a local dataframe to update
    df = writeStatsToFrame(C)
    # write this stats to a csv file
    CSV.write("data/"*string(timesInSec)*".csv", df)

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

# this function gets random position based on the type of the galaxy
function getRandomPosition(type::String, num, vector::Array{Float64,1}, minn::Int64,maxx::Int64)
    if type == "random"
        position = [rand(minn:maxx)+rand(Float64), rand(minn:maxx)+rand(Float64) ,rand(minn:maxx)+rand(Float64)]
    elseif type == "disk"
        normal = vector
        x = 0.0
        y = 0.0
        z = 0.0
        if normal[3]==0 && normal[1]!= 0 && normal[2]!= 0
            x = rand(minn:maxx)+rand(Float64)
            y = (0-normal[1]*x)/normal[2]
            z = rand(minn:maxx)+rand(Float64)
        elseif normal[3]==0 && normal[1]== 0 && normal[2]!= 0
            x = rand(minn:maxx)+rand(Float64)
            y = (0-normal[1]*x)/normal[2]
            z = rand(minn:maxx)+rand(Float64)
        elseif normal[2]==0 && normal[1]!= 0 && normal[3]!= 0
            x = rand(minn:maxx)+rand(Float64)
            y = rand(minn:maxx)+rand(Float64)
            z = (0-normal[1]*x)/normal[3]
        elseif normal[1]==0 && normal[2]!= 0 && normal[3]!= 0
            x = rand(minn:maxx)+rand(Float64)
            y = rand(minn:maxx)+rand(Float64)
            z = (0-normal[2]*y)/normal[3]
        else
            x = rand(minn:maxx)+rand(Float64)
            y = rand(minn:maxx)+rand(Float64)
            z = (0-x*normal[1]-y*normal[2])/normal[3]
        end
        position = [x,y,z]
    elseif type == "line"
        dir = vector
        x = 0.0
        y = 0.0
        z = 0.0
        if dir[1]!=0
            x = rand(minn:maxx)+rand(Float64)
            y = dir[2]*x/dir[1]
            z = dir[3]*x/dir[1]
        elseif dir[2]!=0
            y = rand(minn:maxx)+rand(Float64)
            x = dir[1]*y/dir[2]
            z = dir[3]*y/dir[2]
        end
        position = [x,y,z]
    end
    return position
end
# this function is a helper function that adds bodies based on different system
function addBodies(C::Node, type::String, num::Int64, vector::Array{Float64,1}, minn::Int64,maxx::Int64)
    mid = (maxx - minn)/5
    if type == "random"
        # add random bodies
        for i in 1:Int(floor(0.9num))
            mass = exp(randn())*M
            position = getRandomPosition("random", num, [0.0,0.0,0.0], minn,maxx)
            newBody = Body(i, mass, position, SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
            insertBody(C,newBody)
        end
        for i in 1:(num - Int(floor(0.9num)))-1
            mass = (rand(1:12)+rand(Float64))*10M
            position = getRandomPosition("random", num, [0.0,0.0,0.0], Int64(ceil(minn+mid)),Int64(floor(maxx-mid)))
            newBody = Body(i, mass, position, SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
            insertBody(C,newBody)
        end
        mass = (rand(1:10)+rand(Float64))*50000M
        position = rand()/100 .+C.centerOfMass
        newBody = Body(num, mass, position, SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
        insertBody(C,newBody)

    elseif type == "disk"
        for i in 1:Int(floor(0.9num))
            mass = exp(randn())*M
            position = getRandomPosition("disk", num, vector, minn,maxx)
            newBody = Body(i, mass, position, SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
            insertBody(C,newBody)
        end
        for i in 1:(num - Int(floor(0.9num)))-1
            mass = (rand(1:10)+rand(Float64))*10M
            position = getRandomPosition("disk", num, vector, Int64(ceil(minn+mid)),Int64(floor(maxx-mid)))
            newBody = Body(i, mass, position, SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
            insertBody(C,newBody)
        end

        mass = (rand(1:10)+rand(Float64))*50000M
        println("CENTER OF MASS BEFORE ", C.centerOfMass)

        println("total mass before, ", C.mass)
        position = rand()/100 .+C.centerOfMass
        println("position is, ", position)
        println("big mass is, ", mass)
        newBody = Body(-1, mass, position, SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
        println("Force is ", G*C.mass*mass./(dis(C.centerOfMass,position)*distanceRate)^2)
        insertBody(C,newBody)
        println("CENTER OF MASS AFTER ", C.centerOfMass)
        println("total mass after, ", C.mass)


    elseif type == "line"
        for i in 1:Int(floor(0.9num))
            mass = exp(randn())*M
            position = getRandomPosition("line", num, vector, minn,maxx)
            newBody = Body(i, mass, position, SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
            insertBody(C,newBody)
        end
        for i in 1:(num - Int(floor(0.9num)))-1
            mass = (rand(1:10)+rand(Float64))*10M
            position = getRandomPosition("line", num, vector, Int64(ceil(minn+mid)),Int64(floor(maxx-mid)))
            newBody = Body(i, mass, position, SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
            insertBody(C,newBody)
        end
        mass = (rand(1:10)+rand(Float64))*100M
        position = getRandomPosition("line", num, vector, Int64(ceil(minn+2mid)),Int64(floor(maxx-2mid)))
        newBody = Body(num, mass, position, SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
        insertBody(C,newBody)
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
        global timeScalar = 60*60*60*24*365*7e14
        global strengthOfInteraction = nothing
        global strengthOfInteraction = 2e33
    end
    basicSimulation(C)

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
    changeRotatingPlane(normal)
    #caluculateInitialSpeed(C)
    caluculateInitialSpeed(C)
    if num== 2
        global timeScalar = nothing
        global timeScalar = 60*60*60*24*365*7e7
        global strengthOfInteraction = nothing
        global strengthOfInteraction = 3
    else
        global timeScalar = nothing
        global timeScalar = 60*60*60*24*365*7e14
        global strengthOfInteraction = nothing
        global strengthOfInteraction = 2e33
        # or 4e33
    end
    # addEnvironment(C,"disk",10*num,normal,1,1000)
    basicSimulation(C)
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
    basicSimulation(C)
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
    global timeScalar = 60*60*60*24*365*1e7
    global strengthOfInteraction = nothing
    global strengthOfInteraction = 2e50
    basicSimulation(C)
end

function jupiterSring()
    minBounds, maxBounds = getBoundsBasedOnTypes("random",[0.0,0.0,0.0], 499.999999995, 500.000000005)
    # making an empty tree
    C= Node(
           false, Array{Node,1}(undef, 8),0,
           Body(), 0.0, SVector{3, Float64}(zeros(3)),
           minBounds,maxBounds
           )
    # add 100 small random bodies
    for i in 1:100
        mass = rand()*M*5e-17
        position = [500+rand()/250000000,500+rand()/250000000,500+rand()/250000000 ]
        newBody = Body(i, mass, position, SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
        insertBody(C,newBody)
    end

    # add jupitar
    newBody = Body(0, 0.001M, [500.0,500.0,500.0], SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
    insertBody(C,newBody)
    # assign initial velocity for all the body
    caluculateInitialSpeed(C)

    # tuning for different cases

    global timeScalar = nothing
    global timeScalar = 60*60*60
    global strengthOfInteraction = nothing
    global strengthOfInteraction = 1

    basicSimulation(C)

end

# this function add 10X random stars
function addEnvironment(C::Node, type::String, num::Int64, vector::Array{Float64,1}, minn::Int64,maxx::Int64)
    for i in 1:num
        mass = rand()*M
        position = getRandomPosition("disk", num, vector, minn,maxx)
        newBody = Body(i, mass, position, SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))
        insertBody(C,newBody)
    end
    printSpeed(C)
    caluculateInitialSpeedForEnvironment(C, vector)

end
# this function changes theta (s/d)
function changeTheta(theta::Float64)
    global THETA = nothing
    global THETA = theta
end

# this function changes how long the simulation is
function changeDuration(factor::Float64)
    oldSecs =  secs
    global secs = nothing
    global secs = oldSecs*factor
    global secsToPlot = nothing
    global secsToPlot = oldSecs*factor

end

# this function changes how long the simulation is
function changeDuration(factor::Int64)
    oldSecs =  secs
    global secs = nothing
    global secs = oldSecs*factor
    global secsToPlot = nothing
    global secsToPlot = oldSecs*factor

end

function getDuration()
    return secsToPlot
end

function changeSpeed(factor::Float64)
    oldSpeed = timeScalar
    global timeScalar = nothing
    global timeScalar = oldSpeed*factor
end

function changeRotatingPlane(vector::Array{Float64,1})
    if length(vector) != 3
        error("vector must be 3-dimensional")
    end
    global rotatingPlane = nothing
    global rotatingPlane = vector
end

function printSpeed(C::Node)
    if C.hasChildren && C.numOfChild != 1
        for i in 1:8
            printSpeed(C.children[i])
        end
    elseif C.hasChildren && C.numOfChild == 1
        println(C.body.velocity)
    end
end
