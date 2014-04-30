###############################################################################
## tstPin_Cell.py
## 9te
## Mon Dec 13 11:30:07 2010
## $Id: template.py,v 1.3 2008/01/02 17:18:47 9te Exp $
###############################################################################
## Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
###############################################################################

import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt

##---------------------------------------------------------------------------##

class Pin(object):

    def __init__(self, dx, dy, height=14.28):
        self.d      = np.array([dx, dy, height])
        self.vessel = False

    def add_vessel(self, R, xoff, yoff):
        self.R      = R
        self.vessel = False
        self.off    = np.array([xoff, yoff, 0.0])

        (nearR, farR) = self.nfp(xoff, yoff)

        if self.R > nearR and self.R < farR:
            self.vessel = True

    def angle(self, rp, ro):

        o    = np.zeros((3))
        o[:] = rp[:] - ro[:]
        o    = o / math.sqrt(np.dot(o, o))
        return o

    def d2r(self, r, omega):

        r0 = np.zeros((3))
        np.copyto(r0, r)
        r0 = self.l2g(r0)

        a = omega[0]**2 + omega[1]**2
        b = 2.0 * (omega[1]*r0[1] + omega[0]*r0[0])
        c = r0[0]**2 + r0[1]**2 - self.R**2

        d1 = -1.0
        d2 = -1.0

        discriminant = b**2 - 4.0 * a * c

        if discriminant >= 0.0:
            d1 = (-b + math.sqrt(discriminant))/(2.0 * a)
            d2 = (-b - math.sqrt(discriminant))/(2.0 * a)

        return np.array([d1, d2])

    def d2b(self, r, omega):

        db = np.zeros((3))

        for d in xrange(3):
            if (omega[d] > 0.0):
                db[d] = (self.hi(d) - r[d]) / omega[d]
            else:
                db[d] = (self.lo(d) - r[d]) / omega[d]

        return db

    def l2g(self, point):
        gbl    = np.zeros((3))
        gbl[:] = point[:] + self.off[:]
        return gbl

    def pitch(self, dir):
        return self.d[dir]

    def lo(self, dir):
        if dir == 2:
            return 0.0
        return -self.d[dir] * 0.5

    def hi(self, dir):
        if dir == 2:
            return self.d[2]
        return self.d[dir] * 0.5

    @property
    def has_vessel(self):
        return self.vessel

    @property
    def height(self):
        return self.d[2]

    def output(self, r, omega, dr, db):
        print "r     = Vector(%6.2f, %6.2f, %6.2f);" % (r[0], r[1], r[2])
        print "omega = Vector(%16.12f, %16.12f, %16.12f);" \
            % (omega[0], omega[1], omega[2])
        print "d1    = %16.10e" % dr[0]
        print "d2    = %16.10e" % dr[1]
        print "dx    = %16.10e" % db[0]
        print "dy    = %16.10e" % db[1]
        print "dz    = %16.10e" % db[2]
        print "---------------"

    def plotit(self, r0, rp):
        y  = np.zeros((1000))
        dx = self.R / 1000.0
        x  = np.zeros((1000))

        for i in xrange(1, 1000):
            x[i] = x[i-1] + dx

        y[:] = (self.R**2 - x[:]**2)**0.5

        lob = self.l2g(np.array([self.lo(0), self.lo(1), self.lo(2)]))
        hib = self.l2g(np.array([self.hi(0), self.hi(1), self.lo(2)]))

        r   = np.array([0, self.R])
        lox = np.array([lob[0], lob[0]])
        hix = np.array([hib[0], hib[0]])
        loy = np.array([lob[1], lob[1]])
        hiy = np.array([hib[1], hib[1]])

        xp = np.array([r0[0], rp[0]])
        yp = np.array([r0[1], rp[1]])

        fig   = plt.figure()
        axes  = fig.add_subplot(111)
        qmesh = axes.plot(x, y, xp, yp, lox, r, hix, r, r, loy, r, hiy)
        axes.set_ylim(bottom = 0.0)
        plt.show()

    def plot2(self, r0, rp, R1):
        y  = np.zeros((1000))
        dx = self.R / 1000.0
        x  = np.zeros((1000))

        for i in xrange(1, 1000):
            x[i] = x[i-1] + dx

        y[:] = (self.R**2 - x[:]**2)**0.5

        y1 = np.zeros((1000))
        dx = R1 / 1000.0
        x1 = np.zeros((1000))

        for i in xrange(1, 1000):
            x1[i] = x1[i-1] + dx

        y1[:] = (R1**2 - x1[:]**2)**0.5

        lob = self.l2g(np.array([self.lo(0), self.lo(1), self.lo(2)]))
        hib = self.l2g(np.array([self.hi(0), self.hi(1), self.lo(2)]))

        r   = np.array([0, self.R * 1.5])
        lox = np.array([lob[0], lob[0]])
        hix = np.array([hib[0], hib[0]])
        loy = np.array([lob[1], lob[1]])
        hiy = np.array([hib[1], hib[1]])

        xp = np.array([r0[0], rp[0]])
        yp = np.array([r0[1], rp[1]])

        fig   = plt.figure()
        axes  = fig.add_subplot(111)
        qmesh = axes.plot(x, y, xp, yp, lox, r, hix, r, r, loy, r, hiy,
                          x1, y1)
        axes.set_ylim(bottom = 0.0)
        plt.show()

    def nfp(self, xoff, yoff):

        near = np.zeros((2))
        far  = np.zeros((2))
        off  = np.array([xoff, yoff])

        near[:] = off[:] - self.d[0:2] * 0.5
        far[:]  = off[:] + self.d[0:2] * 0.5

        for dir in xrange(2):
            if off[dir] < 0.0:
                near[dir] = off[dir] + self.d[dir] * 0.5
                far[dir]  = off[dir] - self.d[dir] * 0.5

        nearR = math.sqrt(near[0]**2 + near[1]**2)
        farR  = math.sqrt(far[0]**2 + far[1]**2)

        return (nearR, farR)

##---------------------------------------------------------------------------##
## PIN TESTS
##---------------------------------------------------------------------------##

def init():

    xoff = 21.62 * 5.5
    yoff = 30.00 * 2.5

    pin = Pin(21.62, 30.00)
    (nearR, farR) = pin.nfp(xoff, yoff)
    print nearR, farR

##---------------------------------------------------------------------------##

def track_lor():

    xoff = 21.62 * 5.5
    yoff = 30.00 * 2.5

    pin = Pin(21.62, 30.00)
    (nearR, farR) = pin.nfp(xoff, yoff)
    print nearR, farR

    pin.add_vessel(135.0, xoff, yoff)

    r0    = np.array([-0.52, -0.25, 11.8])
    rp    = np.array([21.62*0.5, 0.01, 7.8])
    omega = pin.angle(rp, r0)

    dr = pin.d2r(r0, omega)
    db = pin.d2b(r0, omega)

    pin.output(r0, omega, dr, db)
    pin.plot2(pin.l2g(r0), pin.l2g(rp), 140.0)

    r0    = np.array([-0.52, -0.25, 11.8])
    rp    = np.array([-21.62*0.5, -12.0, 7.8])
    omega = pin.angle(rp, r0)

    dr = pin.d2r(r0, omega)
    db = pin.d2b(r0, omega)

    pin.output(r0, omega, dr, db)
    pin.plot2(pin.l2g(r0), pin.l2g(rp), 140.0)

    r0    = np.array([-2.52, -4.2, 11.8])
    rp    = np.array([-1.0, -15.0, 7.8])
    omega = pin.angle(rp, r0)

    dr = pin.d2r(r0, omega)
    db = pin.d2b(r0, omega)

    pin.output(r0, omega, dr, db)
    pin.plot2(pin.l2g(r0), pin.l2g(rp), 140.0)

    r0    = np.array([-2.52, -4.2, 11.8])
    rp    = np.array([-1.0, 15.0, 7.8])
    omega = pin.angle(rp, r0)

    dr = pin.d2r(r0, omega)
    db = pin.d2b(r0, omega)

    pin.output(r0, omega, dr, db)
    pin.plot2(pin.l2g(r0), pin.l2g(rp), 140.0)

##---------------------------------------------------------------------------##

def track_hir():

    xoff = 21.62 * 5.5
    yoff = 30.00 * 2.5

    pin = Pin(21.62, 30.00)
    (nearR, farR) = pin.nfp(xoff, yoff)
    print nearR, farR

    pin.add_vessel(140.0, xoff, yoff)

    r0    = np.array([-0.52, -0.25, 11.8])
    rp    = np.array([21.62*0.5, 0.01, 7.8])
    omega = pin.angle(rp, r0)

    dr = pin.d2r(r0, omega)
    db = pin.d2b(r0, omega)

    pin.output(r0, omega, dr, db)
    pin.plotit(pin.l2g(r0), pin.l2g(rp))

    r0    = np.array([-0.52, -0.25, 11.8])
    rp    = np.array([-21.62*0.5, -12.0, 7.8])
    omega = pin.angle(rp, r0)

    dr = pin.d2r(r0, omega)
    db = pin.d2b(r0, omega)

    pin.output(r0, omega, dr, db)
    pin.plotit(pin.l2g(r0), pin.l2g(rp))

    r0    = np.array([-2.52, -4.2, 11.8])
    rp    = np.array([-1.0, -15.0, 7.8])
    omega = pin.angle(rp, r0)

    dr = pin.d2r(r0, omega)
    db = pin.d2b(r0, omega)

    pin.output(r0, omega, dr, db)
    pin.plotit(pin.l2g(r0), pin.l2g(rp))

    r0    = np.array([-2.52, -4.2, 11.8])
    rp    = np.array([-1.0, 15.0, 7.8])
    omega = pin.angle(rp, r0)

    dr = pin.d2r(r0, omega)
    db = pin.d2b(r0, omega)

    pin.output(r0, omega, dr, db)
    pin.plotit(pin.l2g(r0), pin.l2g(rp))

##---------------------------------------------------------------------------##

def track_hi2lo():

    xoff = 21.62 * 5.5
    yoff = 30.00 * 2.5

    pin = Pin(21.62, 30.00)

    pin.add_vessel(140.0, xoff, yoff)

    r0    = np.array([-0.52, -0.25, 11.8])
    rp    = np.array([-21.62*0.5, -12.0, 7.8])
    omega = pin.angle(rp, r0)

    dr = pin.d2r(r0, omega)
    db = pin.d2b(r0, omega)

    pin.output(r0, omega, dr, db)
    pin.plot2(pin.l2g(r0), pin.l2g(rp), 135.0)

    r0[:] += omega[:] * dr[1]

    pin.add_vessel(135, xoff, yoff)

    dr = pin.d2r(r0, omega)
    db = pin.d2b(r0, omega)

    pin.output(r0, omega, dr, db)
    pin.plot2(pin.l2g(r0), pin.l2g(rp), 140.0)

    r0[:] += omega[:] * dr[1]

    dr = pin.d2r(r0, omega)
    db = pin.d2b(r0, omega)

    pin.output(r0, omega, dr, db)
    pin.plot2(pin.l2g(r0), pin.l2g(rp), 140.0)

##---------------------------------------------------------------------------##

def track_lo2hi():

    xoff = 21.62 * 5.5
    yoff = 30.00 * 2.5

    pin = Pin(21.62, 30.00)

    pin.add_vessel(135.0, xoff, yoff)

    r0    = np.array([-10.0, 0.25, 11.8])
    rp    = np.array([6.0, 15.0, 7.8])
    omega = pin.angle(rp, r0)

    dr = pin.d2r(r0, omega)
    db = pin.d2b(r0, omega)

    pin.output(r0, omega, dr, db)
    pin.plot2(pin.l2g(r0), pin.l2g(rp), 140.0)

    r0[:] += omega[:] * dr[0]

    pin.add_vessel(140.0, xoff, yoff)

    dr = pin.d2r(r0, omega)
    db = pin.d2b(r0, omega)

    pin.output(r0, omega, dr, db)
    pin.plot2(pin.l2g(r0), pin.l2g(rp), 135.0)

    r0[:] += omega[:] * dr[0]

    dr = pin.d2r(r0, omega)
    db = pin.d2b(r0, omega)

    pin.output(r0, omega, dr, db)
    pin.plot2(pin.l2g(r0), pin.l2g(rp), 135.0)

##---------------------------------------------------------------------------##

def track_assembly():

    r0 = np.array([1.4, 2.25, 5.0])
    rp = np.array([1.4, -2.25, 5.0])

    pin   = Pin(4.5, 4.5)
    omega = pin.angle(rp, r0)
    pin.add_vessel(6.5, 4.5, 4.5)

    dr = pin.d2r(r0, omega)
    db = pin.d2b(r0, omega)

    pin.output(r0, omega, dr, db)
    pin.plot2(pin.l2g(r0), pin.l2g(rp), 6.0)

    r0[:] += omega[:] * dr[1]

    dr = pin.d2r(r0, omega)
    db = pin.d2b(r0, omega)

    pin.output(r0, omega, dr, db)
    pin.plot2(pin.l2g(r0), pin.l2g(rp), 6.0)

    ##------

    r0 = np.array([1.4, 0.75, 5.0])
    rp = np.array([1.4, -0.75, 5.0])

    pin = Pin(4.5, 1.5)
    omega = pin.angle(rp, r0)
    pin.add_vessel(6.0, 4.5, 1.5)

    dr = pin.d2r(r0, omega)
    db = pin.d2b(r0, omega)

    pin.output(r0, omega, dr, db)
    pin.plot2(pin.l2g(r0), pin.l2g(rp), 6.5)

    r0[:] += omega[:] * dr[1]

    dr = pin.d2r(r0, omega)
    db = pin.d2b(r0, omega)

    pin.output(r0, omega, dr, db)
    pin.plot2(pin.l2g(r0), pin.l2g(rp), 6.5)


##---------------------------------------------------------------------------##

# init()
#track_hi2lo()
#track_lo2hi()

track_assembly()

#track_lor()
#track_hir()

###############################################################################
##                            end of tstPin_Cell.py
###############################################################################
