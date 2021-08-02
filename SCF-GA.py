import random

my_file = open("Test.txt", "r")
content = my_file.read()

new_list = content.split("\n")
test_cases = int(new_list[0])
del new_list[0]

new_list2 = []
for i in new_list:
    Str = i.split(" ")
    new_list2.append(Str)
xy_points = []
for i in range(len(new_list2)):
    mat = [float(new_list2[i][0]), float(new_list2[i][1])]
    xy_points.append(mat)


def Random(degreeOf):
    randChrom = []
    for deg in range(degreeOf + 1):
        rand = random.uniform(10, -10)
        randChrom.append(rand)
    return randChrom


def ValidMut(offsp1, offsp2, cross1, cross2):
    for gene in range(len(offsp1)):
        if 10 < offsp1[gene] < -10:
            offsp1[gene] = cross1[gene]
        if 10 < offsp2[gene] < -10:
            offsp2[gene] = cross2[gene]
    return offsp1, offsp2


def Solutions(degreeOf, popSize):
    sol = []
    while len(sol) < popSize:
        randChromosome = Random(degreeOf)
        sol.append(randChromosome)
    return sol


def Predict(x, y, sol, degreeOf):
    deg = degreeOf
    y2 = 0
    for a in sol:
        y2 += a * (x ** deg)
        deg -= 1
    return (y2 - y) ** 2


def MSE(sol, pointsXY, degreeOf):
    Fit = []
    for chrome in sol:
        y = 0
        for point in pointsXY:
            y += Predict(point[0], point[1], chrome, degreeOf)
        Fit.append(len(pointsXY) / y)
    return Fit


def Selection(pointsXY, degreeOf, sol):
    lFit = MSE(sol, pointsXY, degreeOf)
    lPer = []
    lRange = []
    totFit = sum(lFit)
    z = random.random()
    # i = 0
    for indv in lFit:
        per = indv / totFit
        lPer.append(per)

    start = 0
    end = lPer[0]
    for p in lPer:
        lRange.append([start, end])
        start = end
        end = end + p
    count = 0
    for inde, r in enumerate(lRange):
        if r[0] <= z <= r[1]:
            count = inde
    return sol[count]


def Pc():
    x = random.uniform(0.4, 0.7)
    return x


def CrossOver(p1, p2, degreeOf):
    pc = Pc()
    Pos1 = random.randint(1, degreeOf)
    rand = random.random()
    OfSp1 = []
    OfSp2 = []
    if rand <= pc:
        OfSp1 = p1[:Pos1] + p2[Pos1:]
        OfSp2 = p2[:Pos1] + p1[Pos1:]
    else:
        OfSp1 = p1
        OfSp2 = p2
    return [OfSp1, OfSp2]


def Mutation(off, t, T, b):
    offs = off
    lower = -10
    upper = 10
    for g in range(len(offs)):
        deltaL = offs[g] - lower
        deltaU = upper - offs[g]
        r = random.random()
        y = 0
        if r <= 0.5:
            y = deltaL
        else:
            y = deltaU
        x = y * (1 - (r ** ((1 - (t / T)) ** b)))
        if y == deltaL:
            offs[g] = offs[g] - x
        else:
            offs[g] = x - offs[g]
    return offs


def Replacement(aftMut, pointsXY, sol, degreeOf):
    replace = []
    fsol = MSE(sol, pointsXY, degreeOf)
    faftMut = MSE(aftMut, pointsXY, degreeOf)
    replacement = []
    for i in range(len(fsol)):
        for j in range(len(faftMut)):
            if faftMut[j] >= fsol[i]:
                replacement = aftMut[j]
            else:
                replacement = sol[i]
        replace.append(replacement)
    return replace


def maximum(pointsXY, degreeOf, sol):
    fsol = MSE(sol, pointsXY, degreeOf)
    x = max(fsol)
    z = fsol.index(max(fsol))
    return x, sol[z]


def GA(sol, pointsXY, degreeOf, tG):
    coSol = sol
    mini = maximum(pointsXY, degreeOf, coSol)
    b = random.randint(1, 5)
    for cG in range(tG):
        newGen = []
        sel1 = []
        sel2 = []
        cross = []
        for ind in range(0, len(coSol), 2):
            sel1 = Selection(pointsXY, degreeOf, coSol)
            sel2 = Selection(pointsXY, degreeOf, coSol)
            while sel1 == sel2:
                sel2 = Selection(pointsXY, degreeOf, sol)
            cross = CrossOver(sel1, sel2, degreeOf)
            mut1 = Mutation(cross[0], cG, tG, b)
            mut2 = Mutation(cross[1], cG, tG, b)
            aftCheck = ValidMut(mut1, mut2, cross[1], cross[0])
            newGen.append(aftCheck[0])
            newGen.append(aftCheck[1])
        rep = Replacement(newGen, pointsXY, coSol, degreeOf)
        coSol = rep
        minx = maximum(pointsXY, degreeOf, coSol)
        if mini[0] <= minx[0]:
            mini = minx
    return mini[1], mini[0]


GenNum = 100
PopSize = 700

file1 = open("out.txt", "w")
for i in range(test_cases):
    print("case: ", i + 1)
    file1.write("case: ")
    file1.write(str((i + 1)))
    file1.write("\n")
    points = int(xy_points[0][0])
    print("points: ", points)
    degree = int(xy_points[0][1])
    print("degree: ", degree)
    del xy_points[0]
    XY = []
    for k in range(points):
        XY.append(xy_points[k])
    del xy_points[0:points]
    solFin = Solutions(degree, PopSize)
    rslt = GA(solFin, XY, degree, GenNum)
    file1.write(str(rslt[0]))
    file1.write(" , ")
    file1.write(str(1.0 / rslt[1]))
    file1.write("\n")

file1.close()
