import random
import math


def anglepoint(rad, theta):
    return [round(rad * math.cos(math.radians(theta)), 2), round(rad * math.sin(math.radians(theta)), 2)]


inx = 40
a = 10
r = .5

drawlines = False
drawpoints = False
drawcircle = False
drawblob = True
dashedlines = False

print('\\documentclass{article}\n\\usepackage{tikz}\\begin{document}')

for p in range(45):
    angles = []
    for i in range(360 // inx):
        theta = inx * i + random.gauss(0, a)
        angles += [theta]

    radangles = [[random.gauss(2, r), angles[i]] for i in range(360 // inx)]

    points = [[tuple(anglepoint(*point)), round(90 + point[1]), round(point[1])] for point in radangles]

    points += [points[0]]

    xscale = round(1 + (random.uniform(0, 1)) ** 2, 2)
    xscale = 1
    dashed = '[dashed]' if dashedlines else ''

    stringout = '\\begin{tikzpicture}[yscale=1,xscale=' + str(xscale) + ']'
    if drawblob:
        stringout = stringout + '\n\draw' + dashed + ' '
        for i in range(360 // inx):
            stringout = stringout + str(points[i][0]) + ' to[out=' + str(points[i][1]) + ',in=' + str(
                points[i + 1][1] + 180) + '] '
        stringout = stringout + str(points[0][0]) + ';\n'
    if drawlines:
        for x in points:
            stringout = stringout + '\n\\draw (0,0) -- ' + str(x[0]) + ';'
    if drawpoints:
        stringout += '\n\\filldraw (0,0) circle(1pt);'
        for x in points:
            stringout = stringout + '\n\\filldraw' + str(x[0]) + ' circle(1pt);'
    if drawcircle:
        stringout += '\n\\draw[dashed] (0,0) circle[radius=2];'
    stringout += '\end{tikzpicture}'

    print(stringout)

print('\end{document}')