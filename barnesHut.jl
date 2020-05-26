using StaticArrays
using Revise
using LinearAlgebra
using CSV
using DataFrames

global G = 6.67e-11 # m^3/(kg*s^2)

global THETA = 0.7

global timesInSec = 0.00

global frameRate = 1/60

global M = 2e30 # 1 = 1 solar mass

global timeScalar = 60*60*60*24*365*5e6 #1000000 year in frames

global distanceRate = 3.3e16 # 1 = 1 parsec = 3.3 light year = 3.3e16 m

global strengthOfInteraction = 1

global rotatingPlane = [1.0,-1.0,-1.0]
mutable struct Body
    id::Int64
    mass::Float64
    position::SVector{3, Float64}
    force::SVector{3, Float64}
    velocity::SVector{3, Float64}
end
Body() = Body(0, 0.0, SVector{3, Float64}(zeros(Float64,3)),
SVector{3, Float64}(zeros(Float64,3)), SVector{3, Float64}(zeros(Float64,3)))


mutable struct Node
    hasChildren::Bool
    children::Array{Node,1}
    numOfChild::Int64
    body::Body
    mass::Float64
    centerOfMass::SVector{3, Float64}
    minBounds::SVector{3, Float64}
    maxBounds::SVector{3, Float64}
end

Node() = Node(
false, Array{Node,1}(undef, 8),0,
Body(), 0.0, SVector{3, Float64}(zeros(3)),
SVector{3, Float64}(zeros(3)),SVector{3, Float64}(zeros(3))
)



#this function returns which octants will the position belong to in the area defined by minBounds and maxBounds
# this function will output an integer: 1--> top, left, back ; 1--> top, left, back; 2--> top, right, back
#   3--> top, left, front ; 4--> top, right, front; 5--> down, left, back; 6--> down, right, back ;
# 7--> down, left, front; 1--> down, right, front
function determineWhichoctants(position::SArray{Tuple{3},Float64,1,3},
                                minBounds::SArray{Tuple{3},Float64,1,3},
                                maxBounds::SArray{Tuple{3},Float64,1,3})
    x = position[1]
    y = position[2]
    z = position[3]
    midX = (maxBounds[1] - minBounds[1])/2 + minBounds[1]
    midY = (maxBounds[2] - minBounds[2])/2 + minBounds[2]
    midZ = (maxBounds[3] - minBounds[3])/2 + minBounds[3]
    maxX = maxBounds[1]
    maxY = maxBounds[2]
    maxZ = maxBounds[3]
    minX = minBounds[1]
    minY = minBounds[2]
    minZ = minBounds[3]

    if minX <= x < midX && minY <= y < midY && minZ <= z < midZ
        return 5
    elseif minX <= x < midX && minY <= y < midY && maxZ >= z >= midZ
        return 7
    elseif minX <= x < midX && maxY >= y >= midY && minZ <= z < midZ
        return 1
    elseif minX <= x < midX &&  maxY >= y >= midY && maxZ>= z >= midZ
        return 3
    elseif maxX >= x >= midX && minY <= y < midY && minZ <= z < midZ
        return 6
    elseif maxX >= x >= midX && minY <= y < midY && maxZ >= z >= midZ
        return 8
    elseif maxX >= x >= midX && maxY >= y >= midY && minZ <= z < midZ
        return 2
    elseif maxX >= x >= midX && maxY >= y >= midY && maxZ >= z >= midZ
        return 4
    end
end

# this function returns two lists: a list of minBounds for 8 octantss, and a list of maxBounds for 8 octantss
function createSubBounds(minBounds::SVector{3, Float64},maxBounds::SVector{3, Float64})
    # ::Array{SArray{Tuple{3},Float64,1,3},1},::Array{SArray{Tuple{3},Float64,1,3},1}
    minX = minBounds[1]
    maxX = maxBounds[1]
    midX = minX + (maxX-minX)/2
    minY = minBounds[2]
    maxY = maxBounds[2]
    midY = minY + (maxY-minY)/2
    minZ = minBounds[3]
    maxZ = maxBounds[3]
    midZ = minZ + (maxZ-minZ)/2
    quadMinBounds = [[minX,midY,minZ],
    [midX,midY,minZ],
    [minX,midY,midZ],
    [midX,midY,midZ],
    [minX,minY,minZ],
    [midX,minY,minZ],
    [minX,minY,midZ],
    [midX,minY,midZ]]
    quadMaxBounds = [[midX,maxY,midZ],
    [maxX,maxY,midZ],
    [midX,maxY,maxZ],
    [maxX,maxY,maxZ],
    [midX,midY,midZ],
    [maxX,midY,midZ],
    [midX,midY,maxZ],
    [maxX,midY,maxZ]]
    return quadMinBounds, quadMaxBounds
end

# this function divides a space into 8 octantss and makes a new node for each octants
function divideToOctants(node::Node)
    quadMinBounds, quadMaxBound = createSubBounds(node.minBounds, node.maxBounds)
    for i in 1:8
        newNode = Node(
        false, Array{Node,1}(undef, 8),0,
        Body(), 0.0, SVector{3, Float64}(zeros(Float64,3)),
        quadMinBounds[i],quadMaxBound[i]
        )
        node.children[i] = newNode
    end
end

# this functions simulates adding a body into a node, subjecting to the quadtree structure
function insertBody(node::Node, body::Body)

    # which octants does the original body locates
    octantsOriginal = determineWhichoctants(node.centerOfMass, node.minBounds, node.maxBounds)
    # which octants does the new body locates
    octantsNew = determineWhichoctants(body.position, node.minBounds, node.maxBounds)
    # get the bounds for different Octants
    quadMinBounds, quadMaxBound = createSubBounds(node.minBounds, node.maxBounds)

    # update mass, centerOfMass
    totalMass = node.mass + body.mass
    node.centerOfMass = node.mass/totalMass.*node.centerOfMass + body.mass/totalMass.*body.position
    node.mass = totalMass
    # where to insert the body?
    if node.hasChildren && node.numOfChild == 1 # we need to make new subnodes

        node.numOfChild += 1
        divideToOctants(node)
        insertBody(node.children[octantsOriginal], node.body)
        insertBody(node.children[octantsNew],body)
    elseif node.numOfChild > 1 # we already have subnodes
        insertBody(node.children[octantsNew], body)
        node.numOfChild += 1
    else # this is a previously empty node
        node.numOfChild += 1
        node.body = body
        node.centerOfMass = body.position
        node.mass = body.mass
        node.hasChildren = true
    end
end

# this function returns the norm given two coordinates
function dis(x::SVector{3, Float64}, y::SVector{3, Float64})
    return distanceRate*sqrt((x[1]-y[1])^2 + (x[2]-y[2])^2 + (x[3]-y[3])^2)
end

# this function returns the direction subject to the plane of rotation (x-y)
function getDirectionOfVelocity(x::SVector{3, Float64}, y::SVector{3, Float64})
    displacement = distanceRate.*[y[1]-x[1], y[2]-x[2],y[3]-x[3]]
    velocityDir = cross(rotatingPlane, displacement)
    return velocityDir/dis(x,y)
end

# this function calculate the force of a body acted by a node
function calculateForce(node::Node, body::Body)
    currentForce =  calculateForceHelper(node, body)
    body.force = strengthOfInteraction*currentForce
    return body.force
end

function calculateForceHelper(node::Node, body::Body)
    currentForce = [0.0,0.0,0.0]
    if node.hasChildren && node.numOfChild > 1
        # calculate the s/d ratio
        d = dis(node.centerOfMass,body.position)
        delx = (node.maxBounds[1] - node.minBounds[1])*distanceRate
        dely = (node.maxBounds[2] - node.minBounds[2])*distanceRate
        delz = (node.maxBounds[3] - node.minBounds[3])*distanceRate
        s = max(delx,dely,delz)
        if s/d > THETA
            # this means that the node is considered close to the body
            for i in 1:8
                currentForce = currentForce .+ calculateForce(node.children[i], body)
            end
        else
            # this means that the node is far away from the body
            currentForce = currentForce .+ G*node.mass*body.mass./(dis(node.centerOfMass,body.position)^3 )*distanceRate.*(node.centerOfMass .- body.position)
        end
    elseif node.hasChildren && node.numOfChild == 1 && node.body.id!=body.id
        if node.body.position == body.position
            println("crap")
        end
        currentForce = currentForce .+ G*node.mass*body.mass./(dis(node.centerOfMass,body.position)^3)*distanceRate .*(node.body.position .- body.position) # try switching direction of force
    end
    return currentForce
end

# this function first calls apply force to get updates on the new net force acted on the body,
# and then call creepHelper to update the new position and velocity for each body
function creep(node::Node)
    applyForce(node)
    creepHelper(node)
end

function creepHelper(node::Node)
    if node.hasChildren && node.numOfChild > 1
        for i in 1:8
            creepHelper(node.children[i])
        end
    elseif node.hasChildren && node.numOfChild == 1
        acceleration = node.body.force/node.body.mass
        time = frameRate*timeScalar
        newPosition = node.body.position*distanceRate + time .* node.body.velocity + 0.5*time^2 .*acceleration
        node.body.position = newPosition/distanceRate
        velocity = node.body.velocity .+ time .* acceleration
        speed = sqrt(velocity[1]^2+velocity[2]^2+velocity[3]^2)
        velocityDir = velocity / speed
        if speed >= 3e8
            velocity = 3e8 .* velocityDir
            println("speed of light!")
        end
        node.body.velocity = velocity
    end
end

# this function traverse the node and calculate and apply the new force of each body
function applyForce(node::Node)
    applyForceHelper(node, node)

end

# this function is a recursive helper for applyForce
function applyForceHelper(originalNode::Node, node::Node)
    if node.hasChildren && node.numOfChild > 1
        for i in 1:8
            applyForceHelper(originalNode, node.children[i])
        end
    elseif node.hasChildren && node.numOfChild == 1
        calculateForce(originalNode, node.body)
    end
end

# this function takes a node write the positions of the children nodes into a file
function writeStatsToFrame(node::Node)
    df = DataFrame( mass = Float64[], positionX = Float64[], positionY = Float64[],positionZ = Float64[],
    velocityX = Float64[], velocityY = Float64[], velocityZ = Float64[],
    forceX = Float64[], forceY = Float64[], forceZ = Float64[] )
    # first entry is the center of mass
    push!(df, [node.mass, node.centerOfMass[1],node.centerOfMass[2],node.centerOfMass[3],
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    writeStatsToFrameHelper(node, df)
    return df
end

# this function is a recursive helper function for writeStatsToFile
function writeStatsToFrameHelper(node::Node, df::DataFrame)
    if node.hasChildren && node.numOfChild > 1
        for i in 1:8
            writeStatsToFrameHelper(node.children[i],df)
        end
    elseif node.hasChildren && node.numOfChild == 1
        # push!(df, [node.body.position,node.body.velocity,node.body.force])
        push!(df, [node.body.mass, node.body.position[1],node.body.position[2],node.body.position[3],
         node.body.velocity[1],node.body.velocity[2],node.body.velocity[3],
          node.body.force[1],node.body.force[2],node.body.force[3]])
    end
end



# this funciton makes a body given the data from the matrix from the DataFrame from the CVS file
function makeBodyFromStats(id,stats::Array{Float64,1})
    body = Body(id, stats[1],
    [stats[2],stats[3], stats[4]],
    [stats[8],stats[9], stats[10]],
    [stats[5],stats[6], stats[7]])
    return body
end


function detectBounds(stats::Matrix{Float64})
    minimumX = Inf
    minimumY = Inf
    minimumZ = Inf
    maximumX = 0
    maximumY = 0
    maximumZ = 0

    row = size(stats)[1]
    for i in 1:row
        if stats[i,2] < minimumX
            minimumX = stats[i,2]
        end
        if stats[i,2] > maximumX
            maximumX = stats[i,2]
        end
        if stats[i,3] < minimumY
            minimumY = stats[i,3]
        end
        if stats[i,3] > maximumY
            maximumY = stats[i,3]
        end
        if stats[i,4] < minimumZ
            minimumZ = stats[i,4]
        end
        if stats[i,4] > maximumZ
            maximumZ = stats[i,4]
        end
    end
    return minimumX-100, minimumY-100, minimumZ-100, maximumX+100, maximumY+100,maximumZ+100
end

# this function make a tree from the csv file
function makeTreeFromFrame(df)
    stats = convert(Matrix{Float64}, df)
    numOfBodies = size(stats)[1]
    # make the first body (essential to initialize the tree)
    firstBody = makeBodyFromStats(0,stats[1,1:10])
    minimumX, minimumY, minimumZ, maximumX, maximumY,maximumZ = detectBounds(stats)

    # initialize the tree
    C= Node(
       false, # hasChildren
       Array{Node,1}(undef, 8), # children
       0, # numOfChild
       Body(), # body
       0.0, # mass
       SVector{3, Float64}(zeros(3)), # centerOfMass
       SVector{3, Float64}(minimumX,minimumY,minimumZ), # minBounds
       SVector{3, Float64}(maximumX,maximumY,maximumZ)) # maxBounds

    # add bodies, skip the center of mass column
    for i in 2:numOfBodies
        body = makeBodyFromStats(i,stats[i,1:10])
        insertBody(C, body)
    end

    return C
end


# this function update the tree of bodies with stable orbit velocity
function caluculateInitialSpeed(node::Node)
    caluculateInitialSpeedHelper(node, node.centerOfMass, node.mass)
end

function caluculateInitialSpeedHelper(node::Node, centerOfMass::SVector{3, Float64}, totalMass::Float64)
    if node.hasChildren && node.numOfChild > 1
        for i in 1:8
            caluculateInitialSpeedHelper(node.children[i], centerOfMass, totalMass)
        end
    elseif node.hasChildren && node.numOfChild == 1
        node.body.velocity = sqrt(G*totalMass*100/dis(centerOfMass,node.body.position)) .*getDirectionOfVelocity(centerOfMass,node.body.position)
    end
end

# this function update the tree of bodies with stable orbit velocity
function caluculateInitialSpeedPlus(node::Node, centerOfMassVelocity::SVector{3, Float64})
    caluculateInitialSpeedPlusHelper(node, node.centerOfMass, node.mass, centerOfMassVelocity)
end

function caluculateInitialSpeedPlusHelper(node::Node, centerOfMass::SVector{3, Float64}, totalMass::Float64, centerOfMassVelocity::SVector{3, Float64})
    if node.hasChildren && node.numOfChild > 1
        for i in 1:8
            caluculateInitialSpeedPlusHelper(node.children[i], centerOfMass, totalMass,centerOfMassVelocity )
        end
    elseif node.hasChildren && node.numOfChild == 1 && node.body.velocity ==[0.0,0.0,0.0]
        node.body.velocity = centerOfMassVelocity + sqrt(G*totalMass/dis(centerOfMass,node.body.position)) .*getDirectionOfVelocity(centerOfMass,node.body.position)
    end
end


# this function update the tree of bodies with stable orbit velocity
function caluculateInitialSpeedForEnvironment(node::Node, vector::Array{Float64,1})
    caluculateInitialSpeedForEnvironmentHelper(node, vector)
end

function caluculateInitialSpeedForEnvironmentHelper(node::Node, vector::Array{Float64,1})
    if node.hasChildren && node.numOfChild > 1
        for i in 1:8
            caluculateInitialSpeedForEnvironmentHelper(node.children[i], vector)
        end
    elseif node.hasChildren && node.numOfChild == 1 && node.body.velocity == [0.0,0.0,0.0]
        println("olaolaoa")
        x = rand()
        println(x)
        if x>0.5

            node.body.velocity =  (rand(0:1000) + rand(Float64)) .* vector
        end
    end
end
