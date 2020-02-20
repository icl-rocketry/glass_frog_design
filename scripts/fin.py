import numpy as np
import matplotlib.pyplot as plt

import pint

ureg = pint.UnitRegistry()

def si(v):
    return v.to_base_units().m

class Fin:
    # for a single fin, accurate if using 3 or 4 fins
    def __init__(self, root_chord, tip_chord, span, tube_diameter, mid_chord_sweep):

        self.c_r = root_chord
        self.c_t = tip_chord
        self.b = span
        self.r_t = tube_diameter
        self.gamma_c = mid_chord_sweep

        self.sonic_speed = 340*ureg.m/ureg.s #assumed constant

    def Aref(self):
        return np.pi*(self.r_t)**2

    def aspect_ratio(self):

        return 2*self.b / (self.c_r + self.c_t)

    def fin_body_correction(self):

        return 1 + 1/self.tau()

    def tau(self):

        return (self.b + self.r_t)/self.r_t

    def x_t(self):

        # down the root chord
        x = self.c_r/2
        # across the fin
        x += self.b * np.tan(self.gamma_c)

        #up the tip chord
        x -= self.c_t/2

        return x

    def CNa(self, mach=0.0):

        K = self.fin_body_correction()
        AR = self.aspect_ratio()
        cr = self.c_r
        ct = self.c_t
        b = self.b
        rt = self.r_t
        cosgc = np.cos(self.gamma_c)

        beta = np.sqrt(1-mach**2)

        CNa = (1/beta) * K * (AR * (cr+ct)*(b)/(rt**2))/(2 + np.sqrt(4+AR/cosgc))

        return CNa

    def dynPressure(self, Vinf, rho=1.225*ureg.kg/ureg.m**3):
        q = 0.5*rho*Vinf**2

        return q.to(ureg.Pa)

    def lift(self, Vinf, V_gust=9*ureg.m/ureg.s, rho=1.225*ureg.kg/ureg.m**3):

        q = self.dynPressure(Vinf=Vinf, rho=rho)

        Aref = self.Aref()
        alpha = np.arctan(V_gust/Vinf)
        mach = Vinf/self.sonic_speed

        N = q * self.CNa(mach=mach)*alpha*Aref

        return N.to(ureg.N)

    def root_bending_moment(self, Vinf, V_gust=9*ureg.m/ureg.s, rho=1.225*ureg.kg/ureg.m**3):

        N = self.lift(Vinf, V_gust=V_gust, rho=rho)

        moment_arm = self.ycp() - self.r_t

        return (N*moment_arm).to(ureg.N*ureg.mm)

    def equivalent_tip_load(self, Vinf, V_gust=9*ureg.m/ureg.s, rho=1.225*ureg.kg/ureg.m**3):

        M = self.root_bending_moment(Vinf=Vinf, V_gust=V_gust, rho=rho)

        b = self.b

        F = M/b

        return F.to(ureg.N)

    def zcp(self):

        xt = self.x_t()
        cr = self.c_r
        ct = self.c_t

        zcp = (xt/3) * ((cr + 2*ct)/(cr + ct)) + (1/6)*(cr + ct - cr*ct / (cr + ct))

        return zcp

    def ycp(self):

        rt = self.r_t
        cr = self.c_r
        ct = self.c_t
        b = self.b

        ycp = rt + (b/3)*(cr + 2*ct)/(cr+ct)

        return ycp

    def area(self):

        return 0.5*(self.c_r + self.c_t) * self.b


    def plot(self):

        rt = si(self.r_t)
        cr = si(self.c_r)
        ct = si(self.c_t)
        xt = si(self.x_t())
        b = si(self.b)

        coords = np.array([[rt, 0], [rt, -cr], [rt+b, -xt - ct], [rt+b, -xt], [rt, 0]])

        plt.plot(coords[:,0], coords[:,1])

        ax = plt.gca()
        ax.set_aspect(1)
        plt.grid(True)
        plt.xlabel('Radial [m]')
        plt.ylabel('Longitudinal [m]')

        ycp = si(self.ycp())
        zcp = si(self.zcp())

        plt.plot(ycp, -zcp, 'x')
