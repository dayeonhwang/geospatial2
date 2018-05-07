import numpy as np
import copy
import math
from helpers6 import *
from pprint import pprint
from collections import deque

def getLinkID(linkIdx, idx):
    for key, value in linkIdx.items():
        if value == idx:
            return key

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
    # print (curvature[0])
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
        raw_dist = curve_dist2 - curve_dist1
        if (raw_dist < 0):
            # To the reference node
            # if (link1.directionOfTravel == 'F'):
            #     return 0
            raw_dist *= -1
        # else:
            # From the reference node
            # if (link1.directionOfTravel == 'T'):
            #     return 0

        dist2 = raw_dist

    else:
        if link1.nrefNodeID == link2.refNodeID:
            # traveling link1 -> link2
            # if (link1.directionOfTravel == 'T'):
            #     return 0
            dist2 = (link1.length - curve_dist1) + curve_dist2
        elif link1.refNodeID == link2.nrefNodeID:
            # traveling link2 -> link1
            # if (link1.directionOfTravel == 'F'):
            #     return 0
            dist2 = (link2.length - curve_dist2) + curve_dist1
            #direction = 'T'
        else:
            notFound = True

    if notFound:
        return 0
    else:
        dist = math.fabs(dist1 - dist2)
        # print("Distance", dist)
        beta = (1/math.log(2))*dist
        # print("Beta", beta)
        return (1/params.beta)*math.exp(-dist/params.beta)

def find_direction(link1, link2, shape_idx1, shape_idx2):
    # F = from ref node, T = towards ref node

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

def find_dist_from_ref_from_link(probe,matched_link):
    # distance from the reference node to the map-matched probe point location on the link in decimal meters

    link = matched_link
    distances = []
    min_shape_idx = 0
    min_shape_idx2 = 0

    # Find shape point closest to probe point
    for j in range(0, len(link.shapeInfo)):
        ref = link.shapeInfo[j][:]
        ref_lat = float(ref[0])
        ref_lon = float(ref[1])
        distances.append(compute_great_circle_distance(probe.latitude,probe.longitude,ref_lat,ref_lon))
        if (distances[j] < min_distance):
            min_distance = distances[j]
            min_shape_idx = j

    # Compute closest point on link to probe point
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
    shape_idx = min_shape_idx

    return compute_curve_dist(link, shape_idx), height, shape_idx

def MapMatchHMM(params, trajectory, links):
    T = len(trajectory)
    # Create matrix that stores joint probability up to that state (states x time)
    score = np.empty((0,T), float)
    #score = []
    backpointers = np.empty((0,T), int)
    # Create index dictionary for road links in score matrix
    linkIdx = {}

    # Initialize current candidates list
    R = []
    R_links = []
    for t in range(0,T-1):
        if t == 0:
            R, R_links = get_candidate_nodes(trajectory[t], links, [])
            L, L_links = get_candidate_nodes(trajectory[t+1], links, R_links, False)
        else:
            L, L_links = get_candidate_nodes(trajectory[t+1], links, R_links, False)

        # print ("Step", t)
        R_links_IDs = []
        L_links_IDs = []
        for x in R_links:
            R_links_IDs.append(x.linkPVID)
        for x in L_links:
            L_links_IDs.append(x.linkPVID)
        # print ("R_links", R_links_IDs)
        # print ("L_links", L_links_IDs)


        if len(R) == 0:
            print ("no candidates")
            return [], 0
        for l in range(0, len(L)):
            eProb = emission(L[l], trajectory[t+1], params)
            tProb = np.zeros(len(R))
            fProb = np.zeros(len(R))
            for r in range(0, len(R)):
                # Initialize scores and backpointer states
                if t == 0 and l == 0:
                    score_entry = np.zeros((1,T), float)
                    score_entry[0,t] = emission(R[r], trajectory[t], params)
                    score = np.append(score, score_entry, axis = 0)
                    back_entry = np.zeros((1,T), int)
                    backpointers = np.append(backpointers, back_entry, axis = 0)
                    linkIdx[R_links[r].linkPVID] = (np.shape(score)[0])-1

                tProb[r] = transition(trajectory[t], trajectory[t+1], R[r], L[l], R_links[r], L_links[l], params)
                fProb[r] = float(tProb[r] * score[linkIdx[R_links[r].linkPVID],t])
                # print("tProb", tProb[r])
                # print("Score Matrix", score)
                # print("Score Index", linkIdx[R_links[r].linkPVID], t)
                # print("Score", score[linkIdx[R_links[r].linkPVID],t])
                # print("fProb", fProb[r])
            maxr, maxri = np.max(fProb), np.argmax(fProb)
            finalScore = emission(L[l], trajectory[t+1], params)*maxr
            # print("Max score", maxr)
            # Add state to the score and backpointer matrix
            if (L_links[l].linkPVID not in linkIdx):
                score_entry = np.zeros((1,T), float)
                score_entry[0, t+1] = finalScore
                score = np.append(score, score_entry, axis = 0)
                back_entry = np.zeros((1,T), int)
                backpointers = np.append(backpointers, back_entry, axis = 0)
                linkIdx[L_links[l].linkPVID] = (np.shape(score)[0])-1
            else:
                score[linkIdx[L_links[l].linkPVID], t+1] = finalScore

            backpointers[linkIdx[L_links[l].linkPVID], t+1] = linkIdx[R_links[maxri].linkPVID]

        # Transfer next candidates to current candidates
        R = copy.deepcopy(L)
        R_links = copy.deepcopy(L_links)

    # Follow backpointers to resolve sequence
    # print("Score Matrix", score)
    # print ("Backpointers", backpointers)
    sequence = []
    sequence_idx = []
    sequence_idx.append(np.argmax(score[:, T-1]))
    sequence.append(getLinkID(linkIdx, sequence_idx[0]))
    for i in range(T-1, 0, -1):
        sequence_idx.append(backpointers[sequence_idx[-1], i])
        sequence.append(getLinkID(linkIdx, sequence_idx[-1]))

    return list(reversed(sequence)), np.max(score[:, T-1])

if __name__ == "__main__":
    print ("Start Time", datetime.now().time())
    probe_csv = open("Partition6467ProbePoints.csv", "r")
    link_csv = open("Partition6467LinkData.csv", "r")
    output_csv = open("Partition6467MatchedPoints.csv", "w")
    probe_reader = csv.reader(probe_csv)
    link_reader = csv.reader(link_csv)
    writer = csv.writer(output_csv)
    probe_rows = list(probe_reader)
    link_rows = list(link_reader)

    # make each probe row as object of class probe
    probe_obj = []
    for i in range (0, 5):
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
    probe_obj.sort(key=lambda x: x.dateTime)
    print ("sorting dataTime")

    traj = make_trajectory(probe_obj) # dict of probe IDs and probe pts
    print ("got traj")

    probe_ids = deque() # queue list of all probe IDs
    for k in traj.keys():
        probe_ids.append(k)
    print ("probe ids:", probe_ids)

    while (probe_ids):
        curr_id = probe_ids[0] # next probe ID to process
        probe_traj = traj[curr_id] # all probe pts with same ID
        curr_probe = probe_traj[0]
        probe_ids.popleft()
        print ("current probe being processed: ")
        pprint(vars(curr_probe))

        params = Params(4.07, 10)
        print ("starting HMM with", curr_id)
        sequence, prob = MapMatchHMM(params, probe_traj, link_obj)
        print ("sequence", sequence)
        print ("probability", prob)
        for i in range (0, len(sequence)):
            distFromRef1, distFromLink1, shape_idx1 = find_dist_from_ref_from_link(curr_probe,sequence[i]) # current link
            distFromRef2, distFromLink2, shape_idx2 = find_dist_from_ref_from_link(curr_probe,sequence[i+1]) # next link
            if i+1 == len(sequence):
                curr_probe.direction = sequence[i].directionOfTravel
            else:
                curr_probe.direction = find_direction(sequence[i], sequence[i+1], shape_idx1, shape_idx2)

            curr_probe.linkPVID = sequence[i].linkPVID
            curr_probe.distFromRef = distFromRef
            curr_probe.distFromLink = distFromLink
            writer.writerows(curr_probe)

    print ("End Time", datetime.now().time())
