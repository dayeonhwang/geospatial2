import numpy as np
import copy
import math
from helpers6 import *
from pprint import pprint
from collections import deque

def get_closest_shape(probe, link):
    shape = link.shapeInfo
    minDist = minLat = minLon = 0
    for i in range(0,len(shape)):
        lat = float(shape[i][0])
        lon = float(shape[i][1])
        currDist = compute_great_circle_distance(probe.latitude,probe.longitude,lat,lon)
        if (minDist > currDist):
            minDist = currDist
            shapeIndex = i
    return [shapeIndex, minLat, minLon]

def get_dist_curvature(curvatureInfo, idx):
    curvature = curvatureInfo[idx-1]
    print (curvature[0])
    dist = float(curvature[0])
    return dist

def emission(link_pt, probe_pt, params):
    tmp = 1/(math.sqrt(2*math.pi)*params.sigma)
    return tmp*math.exp(-0.5*(link_pt[4]/params.sigma)**2)

def compute_curve_dist(link, idx):
    if (idx == 0):
        curve_dist = 0
    elif (idx == len(link.shapeInfo)-1):
        curve_dist = link.length
    else:
        # If no curvature info, assume straight line
        if (link.curvatureInfo[0][0]== ''):
            curve_dist = compute_great_circle_distance(link.shapeInfo[0][0],
                link.shapeInfo[0][1], link.shapeInfo[idx][0],
                link.shapeInfo[idx][1])
        else:
            curve_dist = get_dist_curvature(link.curvatureInfo, idx)
    return curve_dist

def transition(probe1, probe2, candidate1, candidate2, link1, link2, params):
    # x_t,i: lat/lon pt on link1 nearest probe1 (shapeIndex1, minLat1, minLon1)
    # x_t+1,j: lat/lon pt on link2 nearest probe2 (shapeIndex2, minLat2, minLon2)
    # dist2: distance between x_t,j and x_t+1,j
    dist1 = compute_great_circle_distance(probe1.latitude,probe1.longitude,probe2.latitude,probe2.longitude)
    # shapeIndex1, minLat1, minLon1 = get_closest_shape(probe1,link1)
    # shapeIndex2, minLat2, minLon2 = get_closest_shape(probe2,link2)
    notFound = False
    # print ("entire cand1", candidate1)
    # print ("entire cand2", candidate2)
    # print ("link 1 curvature", link1.curvatureInfo)
    # print ("link 2 curvature", link2.curvatureInfo)
    shape_idx1 = candidate1[1]
    shape_idx2 = candidate2[1]
    # print ("shape info 1", link1.shapeInfo)
    # print ("shape info 2", link2.shapeInfo)

    curve_dist1 = compute_curve_dist(link1, shape_idx1)
    curve_dist2 = compute_curve_dist(link2, shape_idx2)

    # Check if it is along the same link
    if (link1.linkPVID == link2.linkPVID):
        dist2 = math.fabs(curve_dist1 - curve_dist2)
        raw_dist = curve_dist2 - curve_dist1
        if raw_dist < 0:
            direction = 'T'
            dist2 = raw_dist * -1
        else:
            direction = 'F'
            dist2 = raw_dist

    else:
        if link1.nrefNodeID == link2.refNodeID:
            # traveling link1 -> link2
            dist2 = (link1.length - curve_dist1) + curve_dist2
            direction = 'F'
        elif link1.refNodeID == link2.nrefNodeID:
            # traveling link2 -> link1
            dist2 = (link2.length - curve_dist2) + curve_dist1
            direction = 'T'
        else:
            notFound = True
            direction = 'F'

    if notFound:
        return 0
    else:
        dist = math.fabs(dist1 - dist2)
        print("Distance", dist)
        beta = (1/math.log(2))*dist
        print("Beta", beta)
        return (1/params.beta)*math.exp(-dist/params.beta)

def find_direction(link1, link2, candidate1, candidate2):
    # F = from ref node, T = towards ref node
    shape_idx1 = candidate1[1]
    shape_idx2 = candidate2[1]

    curve_dist1 = compute_curve_dist(link1, shape_idx1)
    curve_dist2 = compute_curve_dist(link2, shape_idx2)

    if (link1.linkPVID == link2.linkPVID):
        raw_dist = curve_dist2 - curve_dist1
        if raw_dist < 0:
            direction = 'T'
            dist2 = raw_dist * -1
        else:
            direction = 'F'
            dist2 = raw_dist

    else:
        if link1.nrefNodeID == link2.refNodeID:
            # traveling link1 -> link2
            dist2 = (link1.length - curve_dist1) + curve_dist2
            direction = 'F'
        elif link1.refNodeID == link2.nrefNodeID:
            # traveling link2 -> link1
            dist2 = (link2.length - curve_dist2) + curve_dist1
            direction = 'T'
        else:
            notFound = True
            direction = 'F'

    return direction

def find_dist_from_ref(probe, link):
    # distance from the reference node to the map-matched probe point location on the link in decimal meters
    ref = (link.shapeInfo).split('|')[0].split('/')
    dist = compute_great_circle_distance(probe.latitude, probe.longitude, lat2, lon2)
    return dist

def find_dist_from_link(probe,link):
    # perpendicular distance from the map-matched probe point location on the link to the probe point in decimal meters


def MapMatchHMM(params, trajectory, links):
    T = len(trajectory)
    # Create matrix that stores joint probability up to that state (states x time)
    score = []
    sequence = np.zeros(T-1)
    backpointers = np.empty((0,T), int)
    # Create index dictionary for road links in score matrix
    linkIdx = {}
    idx = 0 # running index for dictionary
    # Sequence of directions
    direction = [[]]
    # Sequence of distFromRef
    distFromRef = [[]]
    # Sequence of distFromLink
    distFromLink = [[]]

    # Initialize current candidates list
    R = []
    R_links = []
    for t in range(T-1):
        if t == 0:
            R, R_links = get_candidate_nodes(trajectory[t], links, [])
            L, L_links = get_candidate_nodes(trajectory[t+1], links, R_links, False)
        else:
            L, L_links = get_candidate_nodes(trajectory[t+1], links, R_links, False)

        print ("Step", t)
        R_links_IDs = []
        L_links_IDs = []
        for x in R_links:
            R_links_IDs.append(x.linkPVID)
        for x in L_links:
            L_links_IDs.append(x.linkPVID)
        print ("R_links", R_links_IDs)
        print ("L_links", L_links_IDs)


        if len(R) == 0:
            print ("no candidates")
        for l in range(0, len(L)):
            eProb = emission(L[l], trajectory[t+1], params)
            tProb = np.zeros(len(R))
            fProb = np.zeros(len(R))
            for r in range(0, len(R)):
                # Initialize scores
                if t == 0 and l == 0:
                    entry = np.zeros(T-1)
                    entry[t] = emission(R[r], trajectory[t], params)
                    score.append(entry)
                    linkIdx[R_links[r].linkPVID] = idx
                    idx += 1

                tProb[r] = transition(trajectory[t], trajectory[t+1], R[r], L[l], R_links[r], L_links[l], params)
                # print ("Link Index", linkIdx)
                fProb[r] = tProb[r] * score[linkIdx[R_links[r].linkPVID]][t]

            maxr, maxri = np.max(fProb), np.argmax(fProb)
            sequence[t] = R_links[maxri].linkPVID
            finalScore = emission(L[l], trajectory[t+1], params)*maxr
            # Add state to the score matrix
            if (L_links[l].linkPVID not in linkIdx):
                entry = np.zeros(T-1)
                entry[t] = finalScore
                score.append(entry)
                linkIdx[L_links[l].linkPVID] = idx
                idx += 1
            else:
                score[linkIdx[L_links[l].linkPVID]][t] = finalScore

        # Transfer next candidates to current candidates
        R = copy.deepcopy(L)
        R_links = copy.deepcopy(L_links)
    return score, sequence

if __name__ == "__main__":
    probe_csv = open("Partition6467ProbePoints.csv", "r")
    link_csv = open("Partition6467LinkData.csv", "r")
    output_csv = open("Partition6467MatchedPoints.csv", "w")
    probe_reader = csv.reader(probe_csv)
    link_reader = csv.reader(link_csv)
    writer = csv.writer(output_csv,lineterminator='\n')
    probe_rows = list(probe_reader)
    link_rows = list(link_reader)

    # make each probe row as object of class probe
    probe_obj = []
    for i in range (950, 1000):
        curr_obj = process_probe_point(probe_rows[i])
        probe_obj.append(curr_obj)
    print ("done making probe obj")

    # make each link row as object of class link
    link_obj = []
    for i in range (0, len(link_rows)):
        curr_obj = process_link(link_rows[i])
        link_obj.append(curr_obj)
    print ("done making link_obj")

    #unsorted_probe_obj = copy.deepcopy(probe_obj)
    #probe_obj.sort(key=lambda x: x.dateTime)
    first_probe = probe_obj[0]
    print ("sorting dataTime")

    traj = make_trajectory(probe_obj) # dict of probe IDs and probe pts
    #TODO: sort probe pts by time

    probe_ids = deque() # queue list of all probe IDs
    for k in traj.keys():
        probe_ids.append(k)
    print ("got traj")
    print ("probe ids:", probe_ids)

    while (probe_ids):
        curr_id = probe_ids[0] # next probe ID to process
        curr_probe = [p for p in probe_obj if p.sampleID==curr_id][0] # finding probe with next probe ID
        probe_ids.popleft()
        print ("current probe being processed: ")
        pprint(vars(curr_probe))

        probe_traj = traj[curr_probe.sampleID] # all probe pts with same ID
        params = Params(4.07, 4)
        print ("starting HMM with", curr_id)
        scores, sequences = MapMatchHMM(params, probe_traj, link_obj)
        print ("scores", scores)
        print ("sequences", sequences)