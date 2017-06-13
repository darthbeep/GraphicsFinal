from display import *
from matrix import *
from math import *
from gmath import *
from random import random

def shading(polygons, ind, color, colormod, mod, normal):
	red = (ind * 18) % 255;
	green = (ind * 32) % 255;
	blue = (ind * 90) % 255;
	color = [red, green, blue]
	ka = [colormod[2][mod][0], colormod[2][mod][3], colormod[2][mod][6]]
	kd = [colormod[2][mod][1], colormod[2][mod][4], colormod[2][mod][7]]
	ks = [colormod[2][mod][2], colormod[2][mod][5], colormod[2][mod][8]]
	mag = sqrt((normal[0] * normal[0]) + (normal[1] * normal[1]) + (normal[2] * normal[2]))
	nn = [normal[0]/mag, normal[1]/mag, normal[2]/mag]

	for i in range(3):

		ambient = colormod[0]['c'][i] * ka[i]
		for l in colormod[1]:

			col = colormod[1][l][0:3]
			source = colormod[1][l][3:6]
			#Sorry
			L = [source[0]-polygons[ind][0], source[1]-polygons[ind][1], source[2]-polygons[ind][2]]
			m = sqrt((L[0] * L[0]) + (L[1] * L[1]) + (L[2] * L[2]))
			nl = [L[0]/m, L[1]/m, L[2]/m]
			dp = (nl[0] * nn[0]) + (nl[1] * nn[1]) + (nl[2] * nn[2])
			diffuse = col[i]*kd[i]*dp
			sca = [nn[0] * dp * 2, nn[1] * dp * 2, nn[2] * dp * 2]
			su = [sca[0] - nl[0], sca[1] - nl[1], sca[2] - nl[2]]
			mm = sqrt((su[0] * su[0]) + (su[1] * su[1]) + (su[2] * su[2]))
			nr = [su[0]/mm, su[1]/mm, su[2]/mm]
			nv = [0,0,1]
			dp2 = (nr[0] * nv[0]) + (nr[1] * nv[1]) + (nr[2] * nv[2])
			specular = col[i]*ks[i]*(dp2)**3

			color[i] += ambient + diffuse+ specular
			#print color[i]

	color = [int(min(max(c,0),255)) for c in color]

	return color

def scanline_convert(polygons, i, screen, zbuffer, color, colormod, mod, normal):
	if mod != 'nope':
		color = shading(polygons, i, color, colormod, mod, normal)
	A = polygons[i]
	B = polygons[i+1]
	C = polygons[i+2]
	top = []
	middle = []
	bottom = []
	#There is probably a much more efficient way to do this but I don't really care
	if A[1] > B[1] and A[1] > C[1]:
		top = A
		if B[1] > C[1]:
			middle = B
			bottom = C
		else:
			middle = C
			bottom = B
	elif B[1] > A[1] and B[1] > C[1]:
		top = B
		if A[1] > C[1]:
			middle = A
			bottom = C
		else:
			middle = C
			bottom = A
	else:
		top = C
		if A[1] > B[1]:
			middle = A
			bottom = B
		else:
			middle = B
			bottom = A
	#Due to a lot of failed testing variable names are mixed up and confusing
	tb = float(top[1] - bottom[1])
	mb = float(middle[1]-bottom[1])
	tm = float(top[1]-middle[1])
	T0 = 0
	T1 = 0
	T0Z = 0
	T1Z = 0
	B0 = 0
	B1 = 0
	B0Z = 0
	B1Z = 0
	if tb != 0 and tm != 0:
		T0 = ((top[0]-middle[0])/tm)
		T1 = ((top[0]-bottom[0])/tb)
		T0Z = (top[2]-middle[2])/tm
		T1Z = (top[2]-bottom[2])/tb
	if mb != 0 and tb != 0:
		B0 = (middle[0]-bottom[0])/mb
		B1 = (top[0]-bottom[0])/tb
		B0Z = (middle[2]-bottom[2])/mb
		B1Z = (top[2]-bottom[2])/tb
	zAv = 0#abs(top[2]-bottom[2])
	C2=0
	#if hh != 0:
	#	C2 = (x2)/hh
	for n in range(int(math.ceil(tm))+1):
		draw_line(int(top[0]-(T0*n)), int(top[1] - n), top[2]-(T0Z * n), int(top[0]-(T1*n)) -0 , int(top[1] - n), top[2]-(T1Z*n), screen, zbuffer, color)
	for n in range(int(math.ceil(mb))):
		draw_line(int(bottom[0]+(B0*n)), int(bottom[1] + n), bottom[2]+(B0Z*n), int(bottom[0]+(B1*n)) +0, int(bottom[1] + n), bottom[2]+(B1Z*n), screen, zbuffer, color)
	draw_line(int(bottom[0]), int(bottom[1]), int(bottom[2]), int(top[0]), int(top[1]), int(top[2]), screen, zbuffer, color)
	draw_line(int(middle[0]), int(middle[1]), int(middle[2]), int(top[0]), int(top[1]), int(top[2]), screen, zbuffer, color)
	draw_line(int(bottom[0]), int(bottom[1]), int(bottom[2]), int(middle[0]), int(middle[1]), int(middle[2]), screen, zbuffer, color)

	#draw_line(int(bottom[0]+(B0 * mb)), int(math.ceil(middle[1]))-1, middle[2], int(bottom[0]+(B1 * mb)), int(math.ceil(middle[1]))-1, middle[2], screen, zbuffer, color)

	#if i == 3:
	#	print 'start'
	#	print [T, M, B]


def add_polygon( polygons, x0, y0, z0, x1, y1, z1, x2, y2, z2 ):
    add_point(polygons, x0, y0, z0);
    add_point(polygons, x1, y1, z1);
    add_point(polygons, x2, y2, z2);

def draw_polygons( matrix, screen, zbuffer, color, colormod, mod='nope' ):
    #print mod
    if len(matrix) < 2:
        print 'Need at least 3 points to draw'
        return

    point = 0
    while point < len(matrix) - 2:

        normal = calculate_normal(matrix, point)[:]
        #print normal
        if normal[2] > 0:
            #color = [0, 255, 0]
            #c = [int((point * 10) % 255), int((point * 20) % 255), int((point * 30) % 255)]
            scanline_convert(matrix, point, screen, zbuffer, color, colormod, mod, normal)
            #print point
            #draw_line( int(matrix[point][0]),  int(matrix[point][1]), matrix[point][2],int(matrix[point+1][0]), int(matrix[point+1][1]), matrix[point+1][2],screen, zbuffer, color)
            #draw_line( int(matrix[point+2][0]), int(matrix[point+2][1]),matrix[point+2][2],int(matrix[point+1][0]),int(matrix[point+1][1]),matrix[point+1][2] screen, zbuffer, color)
            #draw_line( int(matrix[point][0]),  int(matrix[point][1]),matrix[point][2], int(matrix[point+2][0]), int(matrix[point+2][1]), matrix[point+2][2],screen, zbuffer, color)
            color = [0, 0, 0]
        point+= 3


def add_box( polygons, x, y, z, width, height, depth ):
    x1 = x + width
    y1 = y - height
    z1 = z - depth

    #front
    add_polygon(polygons, x, y, z, x1, y1, z, x1, y, z);
    add_polygon(polygons, x, y, z, x, y1, z, x1, y1, z);

    #back
    add_polygon(polygons, x1, y, z1, x, y1, z1, x, y, z1);
    add_polygon(polygons, x1, y, z1, x1, y1, z1, x, y1, z1);

    #right side
    add_polygon(polygons, x1, y, z, x1, y1, z1, x1, y, z1);
    add_polygon(polygons, x1, y, z, x1, y1, z, x1, y1, z1);
    #left side
    add_polygon(polygons, x, y, z1, x, y1, z, x, y, z);
    add_polygon(polygons, x, y, z1, x, y1, z1, x, y1, z);

    #top
    add_polygon(polygons, x, y, z1, x1, y, z, x1, y, z1);
    add_polygon(polygons, x, y, z1, x, y, z, x1, y, z);
    #bottom
    add_polygon(polygons, x, y1, z, x1, y1, z1, x1, y1, z);
    add_polygon(polygons, x, y1, z, x, y1, z1, x1, y1, z1);

def add_sphere( edges, cx, cy, cz, r, step ):
    points = generate_sphere(cx, cy, cz, r, step)
    num_steps = int(1/step+0.1)

    lat_start = 0
    lat_stop = num_steps
    longt_start = 0
    longt_stop = num_steps

    num_steps+= 1
    for lat in range(lat_start, lat_stop):
        for longt in range(longt_start, longt_stop):

            p0 = lat * (num_steps) + longt
            p1 = p0+1
            p2 = (p1+num_steps) % (num_steps * (num_steps-1))
            p3 = (p0+num_steps) % (num_steps * (num_steps-1))

            if longt != num_steps - 2:
	        add_polygon( edges, points[p0][0],
		             points[p0][1],
		             points[p0][2],
		             points[p1][0],
		             points[p1][1],
		             points[p1][2],
		             points[p2][0],
		             points[p2][1],
		             points[p2][2])
            if longt != 0:
	        add_polygon( edges, points[p0][0],
		             points[p0][1],
		             points[p0][2],
		             points[p2][0],
		             points[p2][1],
		             points[p2][2],
		             points[p3][0],
		             points[p3][1],
		             points[p3][2])

def generate_sphere( cx, cy, cz, r, step ):
    points = []
    num_steps = int(1/step+0.1)

    rot_start = 0
    rot_stop = num_steps
    circ_start = 0
    circ_stop = num_steps

    for rotation in range(rot_start, rot_stop):
        rot = step * rotation
        for circle in range(circ_start, circ_stop+1):
            circ = step * circle

            x = r * math.cos(math.pi * circ) + cx
            y = r * math.sin(math.pi * circ) * math.cos(2*math.pi * rot) + cy
            z = r * math.sin(math.pi * circ) * math.sin(2*math.pi * rot) + cz

            points.append([x, y, z])
            #print 'rotation: %d\tcircle%d'%(rotation, circle)
    return points

def add_torus( edges, cx, cy, cz, r0, r1, step ):
    points = generate_torus(cx, cy, cz, r0, r1, step)
    num_steps = int(1/step+0.1)

    lat_start = 0
    lat_stop = num_steps
    longt_start = 0
    longt_stop = num_steps

    for lat in range(lat_start, lat_stop):
        for longt in range(longt_start, longt_stop):

            p0 = lat * (num_steps) + longt;
            if (longt == num_steps - 1):
	        p1 = p0 - longt;
            else:
	        p1 = p0 + 1;
            p2 = (p1 + num_steps) % (num_steps * num_steps);
            p3 = (p0 + num_steps) % (num_steps * num_steps);

            add_polygon(edges,
                        points[p0][0],
                        points[p0][1],
                        points[p0][2],
                        points[p3][0],
                        points[p3][1],
                        points[p3][2],
                        points[p2][0],
                        points[p2][1],
                        points[p2][2] )
            add_polygon(edges,
                        points[p0][0],
                        points[p0][1],
                        points[p0][2],
                        points[p2][0],
                        points[p2][1],
                        points[p2][2],
                        points[p1][0],
                        points[p1][1],
                        points[p1][2] )

def generate_torus( cx, cy, cz, r0, r1, step ):
    points = []
    num_steps = int(1/step+0.1)

    rot_start = 0
    rot_stop = num_steps
    circ_start = 0
    circ_stop = num_steps

    for rotation in range(rot_start, rot_stop):
        rot = step * rotation
        for circle in range(circ_start, circ_stop):
            circ = step * circle

            x = math.cos(2*math.pi * rot) * (r0 * math.cos(2*math.pi * circ) + r1) + cx;
            y = r0 * math.sin(2*math.pi * circ) + cy;
            z = -1*math.sin(2*math.pi * rot) * (r0 * math.cos(2*math.pi * circ) + r1) + cz;

            points.append([x, y, z])
    return points

def add_circle( points, cx, cy, cz, r, step ):
    x0 = r + cx
    y0 = cy
    t = step

    while t <= 1.00001:
        x1 = r * math.cos(2*math.pi * t) + cx;
        y1 = r * math.sin(2*math.pi * t) + cy;

        add_edge(points, x0, y0, cz, x1, y1, cz)
        x0 = x1
        y0 = y1
        t+= step

def add_curve( points, x0, y0, x1, y1, x2, y2, x3, y3, step, curve_type ):

    xcoefs = generate_curve_coefs(x0, x1, x2, x3, curve_type)[0]
    ycoefs = generate_curve_coefs(y0, y1, y2, y3, curve_type)[0]

    t = step
    while t <= 1.00001:
        x = xcoefs[0] * t*t*t + xcoefs[1] * t*t + xcoefs[2] * t + xcoefs[3]
        y = ycoefs[0] * t*t*t + ycoefs[1] * t*t + ycoefs[2] * t + ycoefs[3]

        add_edge(points, x0, y0, 0, x, y, 0)
        x0 = x
        y0 = y
        t+= step

def draw_lines( matrix, screen, zbuffer, color ):
    if len(matrix) < 2:
        print 'Need at least 2 points to draw'
        return

    point = 0
    while point < len(matrix) - 1:
        draw_line( int(matrix[point][0]),
                   int(matrix[point][1]),
                   matrix[point][2],
                   int(matrix[point+1][0]),
                   int(matrix[point+1][1]),
                   matrix[point+1][2],
                   screen, zbuffer, color)
        point+= 2

def add_edge( matrix, x0, y0, z0, x1, y1, z1 ):
    add_point(matrix, x0, y0, z0)
    add_point(matrix, x1, y1, z1)

def add_point( matrix, x, y, z=0 ):
    matrix.append( [x, y, z, 1] )




def draw_line( x0, y0, z0, x1, y1, z1, screen, zbuffer, color ):

    #swap points if going right -> left
    if x0 > x1:
        xt = x0
        yt = y0
        zt = z0
        x0 = x1
        y0 = y1
        z0 = z1
        x1 = xt
        y1 = yt
        z1 = zt

    x = x0
    y = y0
    z = z0
    A = 2 * (y1 - y0)
    B = -2 * (x1 - x0)

    wide = False
    tall = False

    if ( abs(x1-x0) >= abs(y1 - y0) ): #octants 1/8
        wide = True
        loop_start = x
        loop_end = x1
        dx_east = dx_northeast = 1
        dy_east = 0
        d_east = A
        distance = x1 - x
        if ( A > 0 ): #octant 1
            d = A + B/2
            dy_northeast = 1
            d_northeast = A + B
        else: #octant 8
            d = A - B/2
            dy_northeast = -1
            d_northeast = A - B

    else: #octants 2/7
        tall = True
        dx_east = 0
        dx_northeast = 1
        distance = abs(y1 - y)
        if ( A > 0 ): #octant 2
            d = A/2 + B
            dy_east = dy_northeast = 1
            d_northeast = A + B
            d_east = B
            loop_start = y
            loop_end = y1
        else: #octant 7
            d = A/2 - B
            dy_east = dy_northeast = -1
            d_northeast = A - B
            d_east = -1 * B
            loop_start = y1
            loop_end = y

    dz = (z1-z0)
    if loop_end != loop_start:
        dz = (z1-z)/(loop_end-loop_start)

    while ( loop_start < loop_end ):
        plot( screen, zbuffer, color, x, y, z )
        if ( (wide and ((A > 0 and d > 0) or (A < 0 and d < 0))) or
             (tall and ((A > 0 and d < 0) or (A < 0 and d > 0 )))):
            x+= dx_northeast
            y+= dy_northeast
            d+= d_northeast
            z+=dz
        else:
            x+= dx_east
            y+= dy_east
            z+= dz
            d+= d_east
        loop_start+= 1

    plot( screen, zbuffer, color, x, y, z )
