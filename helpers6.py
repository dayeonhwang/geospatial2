import math
import re
import csv
import copy
from datetime import datetime
from geopy.distance import geodesic
from geopy.distance import great_circle

class Probe:
    def __init__(self, sampleID, dateTime, sourceCode, latitude, longitude, altitude, speed, heading):
        self.sampleID = sampleID
        self.dateTime = dateTime
        self.sourceCode = sourceCode
        self.latitude = latitude
        self.longitude = longitude
        self.altitude = altitude
        self.speed = speed
        self.heading = heading
        self.linkPVID = 0
        self.direction = 'B'
        self.distFromRef = 0
        self.distFromLink = 0

class Link:
    def __init__(self, linkPVID, refNodeID, nrefNodeID, length, directionOfTravel, shapeInfo, curvatureInfo, slopeInfo):
        self.linkPVID = linkPVID
        self.refNodeID = refNodeID
        self.nrefNodeID = nrefNodeID
        self.length = length
        self.directionOfTravel = directionOfTravel
        self.shapeInfo = shapeInfo
        self.curvatureInfo = curvatureInfo
        self.slopeInfo = slopeInfo

class Params:
    def __init__(self, beta, sigma):
        self.beta = beta
        self.sigma = sigma

#compute great-circle distance between two points where each point is (latitude,longitude)
def compute_great_circle_distance(lat1, lon1, lat2, lon2):
    p1 = (lat1, lon1)
    p2 = (lat2, lon2)
    return great_circle(p1,p2).meters

#get candidate link nodes that are within 200 meters distance of probe point
def get_candidate_nodes(probe,link_lst, prev_lst, first = True):
    candidates = []
    links = []

    # Only take the links that are connected/are the current set of links
    if not first:
        linkPVIDs = []
        refNodeIDs = []
        nrefNodeIDs = []
        for l in prev_lst:
            linkPVIDs.append(l.linkPVID)
            refNodeIDs.append(l.refNodeID)
            nrefNodeIDs.append(l.nrefNodeID)
        possible_links = [x for x in link_lst if (x.refNodeID in nrefNodeIDs)
            or (x.nrefNodeID in refNodeIDs) or (x.linkPVID in linkPVIDs)]
    else:
        possible_links = copy.deepcopy(link_lst)

    for i in range(0, len(possible_links)):
        link = possible_links[i]
        distances = []
        min_shape_idx = 0
        min_shape_idx2 = 0
        min_distance = 200
        within_range = False
        # Find shape point closest to probe point
        for j in range(0, len(link.shapeInfo)):
            ref = link.shapeInfo[j][:]
            ref_lat = float(ref[0])
            ref_lon = float(ref[1])
            distances.append(compute_great_circle_distance(probe.latitude,probe.longitude,ref_lat,ref_lon))
            if (distances[j]  <= 200):
                within_range = True
                if (distances[j] < min_distance):
                    min_distance = distances[j]
                    min_shape_idx = j

        # Compute closest point on link to probe point
        if (within_range):
            # Find second closest shape point
            if (min_shape_idx > 0 and min_shape_idx < len(link.shapeInfo)-1):
                if (distances[min_shape_idx - 1] < distances[min_shape_idx + 1]):
                    min_shape_idx2 = min_shape_idx - 1
                else:
                    min_shape_idx2 = min_shape_idx + 1
            elif (min_shape_idx == 0):
                min_shape_idx2 = min_shape_idx + 1
            elif (min_shape_idx == len(link.shapeInfo)):
                min_shape_idx2 = min_shape_idx - 1
            # Use Heron's formula to compute area
            pt1 = link.shapeInfo[min_shape_idx][:]
            pt2 = link.shapeInfo[min_shape_idx2][:]
            dist1p = compute_great_circle_distance(probe.latitude, probe.longitude, float(pt1[0]), float(pt1[1]))
            dist2p = compute_great_circle_distance(probe.latitude, probe.longitude, float(pt2[0]), float(pt2[1]))
            dist12 = compute_great_circle_distance(pt1[0],pt1[1],pt2[0],pt2[1])
            S = (dist1p + dist2p + dist12)/2.0
            A = math.sqrt(S*(S-dist1p)*(S-dist2p)*(S-dist12))
            # Calculate altitude from probe point using computed compute_great_circle_distance
            height = 2*A/dist12
            candidate = [link.linkPVID, min_shape_idx, pt1[0], pt1[1], height]
            candidates.append(candidate)
            links.append(link)

    return candidates, links

# convert time string (e.g. "6/12/2009  6:12:49 AM) to miliseconds
def compute_probe_time(raw_time):
    date_obj = datetime.strptime(raw_time, '%m/%d/%Y  %I:%M:%S %p')
    milisec = date_obj.timestamp()*1000
    return int(milisec) # avoid overflow

def process_probe_point(probe_pt):
    probe_obj = Probe(int(probe_pt[0]), compute_probe_time(probe_pt[1]),
        int(probe_pt[2]), float(probe_pt[3]), float(probe_pt[4]),
        int(probe_pt[5]), int(probe_pt[6]), int(probe_pt[7], ))
    return probe_obj

def process_link(link):
    shapeInfo = []
    shape_pts = link[14]
    shape_lst = shape_pts.split('|')
    for i in range(0,len(shape_lst)):
        shape = shape_lst[i].split('/')
        shapeInfo.append(shape)
    curvatureInfo = []
    curvature_pts = link[15]
    curvature_lst = curvature_pts.split('|')
    for i in range(0,len(curvature_lst)):
        curvature = curvature_lst[i].split('/')
        curvatureInfo.append(curvature)
    slopeInfo = []
    slope_pts = link[16]
    slope_lst = slope_pts.split('|')
    for i in range(0,len(slope_lst)):
        slope = slope_lst[i].split('/')
        slopeInfo.append(slope)
    link_obj = Link(int(link[0]),int(link[1]),int(link[2]),float(link[3]),link[5],shapeInfo,curvatureInfo,slopeInfo)
    return link_obj

# compute slope between probe points
def derive_road_slope(probe1, probe2):
    # change in alt
    # change in dist
    probe_alt1 = probe1.altitude
    probe_alt2 = probe2.altitude
    dist = compute_great_circle_distance(probe1.latitude, probe1.longitude, probe2.latitude, probe2.longitude)
    derived_slope = math.abs(probe_alt1 - probe_alt2)/dist

    return derived_slope

# make trajectory dict {probe_id: list of probe points}
def make_trajectory(probe_obj):
    trajectory = {} # {probe ID: list of probe points}
    currID = probe_obj[0].sampleID
    for i in range(0,len(probe_obj)):
        currID = probe_obj[i].sampleID
        if currID in trajectory:
            trajectory[currID].append(probe_obj[i])
        else:
            trajectory[currID] = [probe_obj[i]]
    return trajectory

# return a trajectory of probe data which has same probe ID
def get_trajectory(traj,probe_id):
    return traj[probe_id]
