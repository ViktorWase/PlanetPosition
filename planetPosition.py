 # This code uses data provided by JPL to calculate the approximate positions
 # of all of the major planets in out solar system.

from math import sin, cos, sqrt, pi, fmod, fabs, atan2
import sys


def check_error(cond, msg):
    if not cond:
        raise PlanetPositionError(msg)


class PlanetPositionError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


def keplers_eq(x, c, tol=1.0e-6):
    """Solves x=y-c*sin(y) numerically"""
    diff = 100000.0
    y = x + sin(x)
    while fabs(diff)>tol:
        d_x = x - (y-c*sin(y))
        d_y = d_x/(1.0-c*cos(y))
        y += d_y
        diff = d_y #This is the diff of the previous iteration, but still.

    return y


def get_julian_from_unix(unixMillisecs):
    """
        Converts from unix time (in milliseconds) to modified Julian Date.
    """
    return unixMillisecs / 86400000.0 + 2440587.5


def get_pos(planet, t, unix_time=False, icrf=False, circular=False):
    """
        Calculates the position of a planet. The position is measured in AU (the mean
        distance from the sun to the earth) and the time in years.
        The position is given in an ecliptic coordinate system, and t=2000 means
        the turn of the latest millenium.

        - planet the name or index of the planet. The index should be 0-indexed.

        - t time measured in years since the birth of Christ. (or milliseconds since
          1 Jan 1970 0:00 if unix_time is used)

       - circular=True uses an approximation of the elliptical orbit using a circle. This
         is slightly faster.

        - unix_time=True means that the input time is given in illiseconds since 1 Jan 1970 0:00.

        - icrf=True changes the coordinat system to an equatoral coordinate system.

        returns: An array of length 3, with the position of the planet.
   """
    parameters = [[0.38709927, 0.20563593, 7.00497902, 252.25032350, 77.45779628, 48.33076593, 0.00000037, \
        0.00001906, -0.00594749, 149472.67411175, 0.16047689, -0.12534081], [0.72333566, 0.00677672,\
        3.39467605, 181.97909950, 131.60246718, 76.67984255, 0.00000390, -0.00004107, -0.00078890, \
        58517.81538729, 0.00268329, -0.27769418], [1.00000261, 0.01671123, -0.00001531, 100.46457166,\
        102.93768193, 0.0, 0.00000562, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0.0], \
        [1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891, 0.00001847, \
        0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343], [5.20288700, 0.04838624, \
        1.30439695, 34.39644051, 14.72847983, 100.47390909, -0.00011607, -0.00013253, -0.00183714, \
        3034.74612775, 0.21252668, 0.20469106],[9.53667594, 0.05386179, 2.48599187, 49.95424423, \
        92.59887831, 113.66242448, -0.00125060, -0.00050991, 0.00193609, 1222.49362201, -0.41897216, \
        -0.28867794], [19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.95427630, 74.01692503, \
        -0.00196176, -0.00004397, -0.00242939, 428.48202785, 0.40805281, 0.04240589], [30.06992276, \
        0.00859048, 1.77004347, -55.12002969, 44.96476227, 131.78422574, 0.00026291, 0.00005105, \
        0.00035372, 218.45945325, -0.32241464, -0.00508664]]

    planet_names = ['mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune']
    planet_name = None
    planet_idx = None
    if (sys.version_info >= (3,0) and isinstance(planet, (int))) or (isinstance(planet, (int))):
        check_error(planet>=0, "The planet index must be positive.")
        check_error(planet<8, "The planet index must be less than 8. It is 0-indexed.")
        planet_name = planet_names[planet]
        planet_idx = planet
    else:
        check_error(isinstance(planet, (str)), "The planet must be an integer or string.")
        is_in_list = False
        planet_lower = planet.lower()
        for i, p in enumerate(planet_names):
            if p.lower() == planet_lower:
                planet_idx = i
                is_in_list = True
                break
        check_error(is_in_list==True, "The planet-string must be the name of a planet!")
    check_error(planet_idx != None, """I'm not sure how you got this bug. This shouldn't happen.
    It can't really happen. If you're reading this it probably did happen, so why don't you do
     me a favor and email me at viktorwase@gmail.com and tell me which input-data you used,
     which version of python you're running and maybe even operating system? I promise to include
     you in my speech when I accept my Academy Award.""")

    if sys.version_info >= (3,0):
        check_error(isinstance(t, (float, int)), "The variable t needs to be a number.")
    else:
        check_error(isinstance(t, (float, int, long)), "The variable t needs to be a number.")

    if unix_time:
        mjd = getJulianFromUnix(t*1000.0)
        T = (mjd-2451544.5)/36525.0
    else:
        T = (t-2000.0)/100.0
    pla_par = parameters[planet_idx]
    conv = pi/180.0
    a = pla_par[0]+T*pla_par[6]
    e = pla_par[1]+T*pla_par[7]
    I = (pla_par[2]+T*pla_par[8])*conv
    L = (pla_par[3]+T*pla_par[9])*conv
    w_line = (pla_par[4]+T*pla_par[10])*conv
    om = (pla_par[5]+T*pla_par[11])*conv

    w = w_line - om
    M = L - w_line #TODO: Add the corrections for jupiter to pluto here.
    M = fmod(M + pi, 2*pi) - pi

    if not circular:
        E = keplers_eq(M, e)
    else:
        E = M

    x = a*(cos(E)-e)
    y = a*sqrt(1.0-e*e)*sin(E)

    cos_I = cos(I)
    cos_om = cos(om)
    sin_om = sin(om)
    sin_w = sin(w)
    pos_ecl = [(cos(w)*cos_om-sin_w*sin_om*cos_I)*x - (sin_w*cos_om+cos(w)*sin_om*cos_I)*y,\
        (cos(w)*sin_om+sin_w*cos_om*cos_I)*x + (-sin_w*sin_om+cos(w)*cos_om*cos_I)*y,
        sin(I)*(sin(w)*x + cos(w)*y)
    ]

    if icrf == False:
        return pos_ecl
    else:
        eps_rad = 23.43928*conv
        pos_icrf = [pos_ecl[0],\
            cos(eps_rad)*pos_ecl[1]-sin(eps_rad)*pos_ecl[2],\
            sin(eps_rad)*pos_ecl[1]+cos(eps_rad)*pos_ecl[2]
            ]
        return pos_icrf

if __name__ == '__main__':
    #Examples
    print(get_pos(2, 2000, icrf=True))
    print(get_pos('jupiter', 2017.101, circular=True))
    print(get_pos("earth", 2100.5, icrf=True, circular=True))
    print(get_pos(0, 2000.7))
