__author__ = 'matthias'

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle

linear = lambda m, x, x1, y1 : m * (x - x1) + y1

def calcTangentAngleGeneric(LineOrigin, CircleCenter, CircleRadius):
    """ This function calculates the two angles under which a line starting at LineOrigin touches a circle.
    Input:
     - LineOrigin: [x y]
     - CircleCenter: [x y]
     - CircleDiameter: r
     Output:
      [higherAngle lowerAngle] in degrees """
    distance = np.sqrt(np.sum((LineOrigin - CircleCenter)**2))
    angle_touching = np.degrees(np.arcsin(CircleRadius/distance))
    angle_to_ring = calcAngle(LineOrigin, CircleCenter)
    touchingAngles = np.hstack([angle_to_ring - angle_touching, angle_to_ring + angle_touching])
    return touchingAngles

def calcAngle(LineOrigin, CircleCenter):
    """ Calculate the angle between two points in the cartesian coordinate system in degrees """
    deltaXY = (LineOrigin - CircleCenter)
    angle = np.arctan(deltaXY[1]/deltaXY[0])
    return np.degrees(angle)

def calcTangentAngles(LineOrigin, CircleCenter, CircleRadius):
    """ Returns an array of angle pairs under which a line going trough LineOrigin touches circle(s) with center at
    CircleCenter and a radius of CircleRadius """
    LineOrigin = np.array(LineOrigin)
    CircleCenter = np.array(CircleCenter)

    nRings = CircleCenter.shape[0]
    angles = np.zeros([nRings, 2])

    for i in range(nRings):
        angles[i, :] = calcTangentAngleGeneric(LineOrigin, CircleCenter[i, :], CircleRadius)

    return angles


def findView(touchingAngles):
    """" Returns an array of True or False values: True if there is a gap to view trough and False if there is no gap"""
    touchingAngles = np.array(touchingAngles)
    n_rings = touchingAngles.shape[0]

    view = np.zeros([n_rings-1])
    for i in range(0, n_rings-1):
        openingAngle = touchingAngles[i, 1] - touchingAngles[i, 0]
        if openingAngle > 0.:
            view[i] = 1
    return view

def getOpeningAngle(viewableAngles):
    """" Calculates the opening angle of an array of angle pairs and returns them as a list  """
    openingAngle = []
    for i in range(viewableAngles.shape[0]):
        openingAngle.append(viewableAngles[i, 1] - viewableAngles[i, 0])

    openingAngle = np.array(openingAngle)
    return openingAngle

def calcTPCIntersectionPoints(LaserPosition, LaserAngles, endpoint):

    m = np.tan(np.radians(LaserAngles)).flatten()

    n_angles = m.shape[0]
    intersection = linear(m, endpoint, LaserPosition[0], LaserPosition[1])
    x = np.ones([n_angles, 1])*endpoint

    intersectionXY = np.vstack([x.flatten(), intersection])

    return intersectionXY

def produceLaserTracks(trackResolution, viewableAngles, start=0, stop=99):

    openingAngle = getOpeningAngle(viewableAngles)

    if stop == 99:
        stop = len(openingAngle)

    shotsInOpening = np.floor(openingAngle/trackResolution)
    stepsInOpeining = openingAngle/(shotsInOpening+1)
    tracks = []
    for i in range(start, stop):
        for k in range(int(shotsInOpening[i])):
            tracks.append(viewableAngles[i, 0] + (k+1)*stepsInOpeining[i])


    return np.array(tracks)

def plotRings(ax, RingPositions, RingDiameter):

    n_rings = RingPositions.shape[0]
    circles = []
    for i in range(n_rings):
        circle = Circle(RingPositions[i,:], RingDiameter)
        circles.append(circle)

    RingCollection = PatchCollection(circles)
    ax.add_collection(RingCollection)

def plotLaserTracks(ax, LaserPosition, LaserAngles, intersection):
    lines = []
    intersectionPointsEnd = calcTPCIntersectionPoints(LaserPosition, LaserAngles, intersection)
    n_angles = intersectionPointsEnd.shape[1]
    for i in range(n_angles):
        lines.append([LaserPosition, intersectionPointsEnd[:, i]])

    LaserTracksCollection = LineCollection(lines)
    ax.add_collection(LaserTracksCollection)


# basic TPC definitions
ring_diameter = 2.5
ring_radius = ring_diameter/2
center_to_center = 4
gap = center_to_center - ring_diameter
tpc_length = 1000
tpc_width = 256

# defining the field cage rings
n_fieldcages = tpc_width/center_to_center

RingPositionsX = np.zeros([1, 64])
RingPositionsY = np.linspace(1, 64, 64)*center_to_center - gap
RingPositionsUpstream = np.vstack([RingPositionsX, RingPositionsY]).transpose()
RingPositionsDownstream = np.vstack([RingPositionsX, RingPositionsY]).transpose() + [tpc_length, 0]

# defining the position of the laser
LaserPositionUpstreamX = -25
LaserPositionUpstreamY = 127
LaserPositionUpstream = np.hstack([LaserPositionUpstreamX, LaserPositionUpstreamY])

LaserPositionDownstreamX = 25
LaserPositionDownstreamY = 127
LaserPositionDownstream = np.hstack([LaserPositionDownstreamX, LaserPositionDownstreamY]) + [tpc_length, 0]

# calculations for the upstream laser mirror
anglesUpstream = calcTangentAngles(LaserPositionUpstream, RingPositionsUpstream, ring_radius)
anglesUpstream = np.roll(anglesUpstream, 1)
viewableUpstream = findView(anglesUpstream)
viewableAnglesUpstream = anglesUpstream[np.flatnonzero(viewableUpstream)]

# calculations for the upstream laser mirror
anglesDownstream = calcTangentAngles(LaserPositionDownstream, RingPositionsDownstream, ring_radius)
anglesDownstream = np.roll(np.fliplr(anglesDownstream), 1)
viewableDownstream = np.logical_not(findView(anglesDownstream))
viewableAnglesDownstream = np.fliplr(anglesDownstream[np.flatnonzero(viewableDownstream)])

# laser angles generation
trackResolution = 0.1

laser_tracks_upstream = produceLaserTracks(trackResolution, viewableAnglesUpstream + 360, 0, 9)
laser_tracks_downstream = produceLaserTracks(trackResolution, viewableAnglesDownstream - 180 + 360, 0, 9)

fig1, ax1 = plt.subplots()
plotRings(ax1, RingPositionsUpstream, ring_radius)
plotRings(ax1, RingPositionsDownstream, ring_radius)
#plotLaserTracks(ax, LaserPositionUpstream, viewableAnglesUpstream, 1000)
#plotLaserTracks(ax, LaserPositionDownstream, viewableAnglesDownstream, 0)
plotLaserTracks(ax1, LaserPositionUpstream, laser_tracks_upstream, 1000)
plotLaserTracks(ax1, LaserPositionDownstream, laser_tracks_downstream, 0)

plt.xlim(-30,1030)
plt.ylim(0, 256)
plt.xlabel("z [cm]")
plt.ylabel("x [cm]")
plt.title("top view of TPC")
ax1.set_aspect(1)
plt.show()


# vertical calculations
tpc_height = 242.6
LaserPositionUpstreamZ = tpc_height/2
halfOpeningAngleVertical = np.degrees((np.pi/2 - np.arctan(LaserPositionDownstreamX/(tpc_height/2)))) - 4
verticalTrackRes = 1

NstepsVertical = np.floor(halfOpeningAngleVertical) / verticalTrackRes
stepsVertical = np.linspace(-halfOpeningAngleVertical, halfOpeningAngleVertical, NstepsVertical)


fig2, ax2 = plt.subplots()

plotLaserTracks(ax2, LaserPositionUpstream, stepsVertical, 1000)
plotLaserTracks(ax2, LaserPositionDownstream, stepsVertical+180, 0)

plt.xlim(-30, 1030)
plt.ylim(0, tpc_height)
plt.xlabel("z [cm]")
plt.ylabel("y [cm]")
plt.title("side view of TPC")
ax2.set_aspect(1)
plt.show()


# save angles to text file

UpstramLaserTracks = []
for verticalIdx in range(len(stepsVertical)):
    for horizontalIdx in range(len(laser_tracks_upstream)):
        pos = [LaserPositionDownstream[1], 0, -LaserPositionUpstream[0], np.mod(stepsVertical[verticalIdx] + 360, 360), laser_tracks_upstream[horizontalIdx]]
        UpstramLaserTracks.append(pos)

DownStreamLaserTracks = []
for verticalIdx in range(len(stepsVertical)):
    for horizontalIdx in range(len(laser_tracks_upstream)):
        pos = [LaserPositionDownstream[1], 0, LaserPositionDownstream[0], stepsVertical[verticalIdx]+180, laser_tracks_downstream[horizontalIdx]]
        DownStreamLaserTracks.append(pos)


format = "%.1f %.1f %.1f %.4f %.4f"
header = "Generated Laser Tracks with the following resolutions: " \
         + str(trackResolution) + " deg (horizontal) / " + str(verticalTrackRes) + " deg (vertical)\n" + \
         "angles measured from beam direction: horizontal (phi): counter-clockwise from beam direction\n " + \
         "                                    vertical (thetha): counter-clockwise (towards top) from beam direction\n" \
         " X[cm], Y[cm], Z[cm], thetha[deg], phi[deg] "

np.savetxt("upstream_laser.txt", UpstramLaserTracks, header=header, fmt=format)
np.savetxt("downstream_laser.txt", DownStreamLaserTracks, header=header, fmt=format)

