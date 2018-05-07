import numpy as np
import copy
import math
from helpers6 import *
from pprint import pprint
from collections import deque
import time

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
    dist = float(curvature[0])
    return dist

def emission(candidate, probe_pt, params):
    tmp = 1/(math.sqrt(2*math.pi)*params.sigma)
    return tmp*math.exp(-0.5*(candidate[2]/params.sigma)**2)

def compute_curve_dist(link, candidate):
    # Compute route distance from reference node of link
    idx = candidate[1]
    delta = candidate[3]
    side = candidate[4]
    if (idx == 0):
        curve_dist = 0 + delta
    elif (idx == len(link.shapeInfo)-1):
        curve_dist = link.length - delta
        if (curve_dist < 0):
            curve_dist = 0
    else:
        # If no curvature info, assume straight line
        if (link.curvatureInfo[0][0]== ''):
            curve_dist = compute_great_circle_distance(link.shapeInfo[0][0],
                link.shapeInfo[0][1], link.shapeInfo[idx][0],
                link.shapeInfo[idx][1])
        else:
            curve_dist = get_dist_curvature(link.curvatureInfo, idx)

        if (side == 1):
            curve_dist += delta
        else:
            curve_dist -= delta
            if (curve_dist < 0):
                if (idx > 1) and (link.curvatureInfo[0][0] !=  ''):
                    curve_dist = get_dist_curvature(link.curvatureInfo, idx-1)
                else:
                    curve_dist = 0
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
    # print ("shape info 1", link1.shapeInfo)
    # print ("shape info 2", link2.shapeInfo)

    curve_dist1 = compute_curve_dist(link1, candidate1)
    curve_dist2 = compute_curve_dist(link2, candidate2)

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

def find_direction(link1, link2, candidate1, candidate2):
    # F = from ref node, T = towards ref node

    curve_dist1 = compute_curve_dist(link1, candidate1)
    curve_dist2 = compute_curve_dist(link2, candidate2)

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

def find_matched_params(probe,matched_link):
    # distance from the reference node to the map-matched probe point location on the link in decimal meters

    links = [matched_link]
    candidates, link_results = get_candidate_nodes(probe, links, [], True, 1000)
    candidate = candidates[0]
    distFromRef = compute_curve_dist(matched_link, candidate)
    distFromLink = candidate[2]

    return distFromRef, distFromLink, candidate

def MapMatchHMM(params, trajectory, links):
    min_dist = 200
    T = len(trajectory)
    # Create matrix that stores joint probability up to that state (states x time)
    score = np.empty((0,T), float)
    #score = []
    backpointers = np.empty((0,T), int)
    # Create index dictionary for road links in score matrix
    linkIdx = {}
    # Create running list of total possible_links whose indices correspond to linkIdx
    candidate_links = []

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
        # R_links_IDs = []
        # L_links_IDs = []
        # for x in R_links:
        #     R_links_IDs.append(x.linkPVID)
        # for x in L_links:
        #     L_links_IDs.append(x.linkPVID)
        # print ("R_links", R_links_IDs)
        # print ("L_links", L_links_IDs)


        while (len(R) == 0):
            # Loosen distance requirement
            min_dist += 200
            print ("no candidates, trying", min_dist)
            R, R_links = get_candidate_nodes(trajectory[t], links, [], True, min_dist)
            if (min_dist > 500):
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
                    candidate_links.append(R_links[r])
                tProb[r] = transition(trajectory[t], trajectory[t+1], R[r], L[l], R_links[r], L_links[l], params)
                fProb[r] = float(tProb[r] * score[linkIdx[R_links[r].linkPVID],t])
                # print("tProb", tProb[r])
                # print("Score Matrix", score)
                # print("Score Index", linkIdx[R_links[r].linkPVID], t)
                # print("Score", score[linkIdx[R_links[r].linkPVID],t])
                # print("fProb", fProb[r])
            maxr, maxri = np.max(fProb), np.argmax(fProb)
            finalScore = eProb*maxr
            # print("Max score", maxr)
            # Add state to the score and backpointer matrix
            if (L_links[l].linkPVID not in linkIdx):
                score_entry = np.zeros((1,T), float)
                score_entry[0, t+1] = finalScore
                score = np.append(score, score_entry, axis = 0)
                back_entry = np.zeros((1,T), int)
                backpointers = np.append(backpointers, back_entry, axis = 0)
                linkIdx[L_links[l].linkPVID] = (np.shape(score)[0])-1
                candidate_links.append(L_links[l])
            else:
                score[linkIdx[L_links[l].linkPVID], t+1] = finalScore

            backpointers[linkIdx[L_links[l].linkPVID], t+1] = linkIdx[R_links[maxri].linkPVID]

        # Transfer next candidates to current candidates
        R = copy.deepcopy(L)
        R_links = copy.deepcopy(L_links)

    print ("Length of candidates_links", len(candidate_links))
    print ("Length of score matrix", np.shape(score)[0])
    # Follow backpointers to resolve sequence
    # print("Score Matrix", score)
    # print ("Backpointers", backpointers)
    if (np.shape(score)[0] == 0):
        return [], 0
    sequence = []
    sequence_idx = []
    sequence_idx.append(np.argmax(score[:, T-1]))
    sequence.append(candidate_links[sequence_idx[0]])
    for i in range(T-1, 0, -1):
        sequence_idx.append(backpointers[sequence_idx[-1], i])
        sequence.append(candidate_links[sequence_idx[-1]])

    return list(reversed(sequence)), np.max(score[:, T-1])

if __name__ == "__main__":
    print ("Start Time", datetime.now().time())
    probe_csv = open("Partition6467ProbePoints.csv", "r")
    link_csv = open("Partition6467LinkData.csv", "r")
    output_csv = open("Partition6467MatchedPoints.csv", "a")
    slope_csv = open("Partition6467SlopeData.csv", "a")
    probe_reader = csv.reader(probe_csv)
    link_reader = csv.reader(link_csv)
    writer = csv.writer(output_csv)
    slope_writer = csv.writer(slope_csv)
    probe_rows = list(probe_reader)
    link_rows = list(link_reader)

    # make each probe row as object of class probe
    probe_obj = []
    for i in range (0, 498):
        curr_obj = process_probe_point(probe_rows[i])
        probe_obj.append(curr_obj)
    print ("done making probe obj")

    #make each link row as object of class link
    link_obj = []
    for i in range (0, len(link_rows)):
        curr_obj = process_link(link_rows[i])
        link_obj.append(curr_obj)
    print ("done making link_obj")

    #make a running list of links that have derived slopes computed
    link_slopes = []
    link_slope_idx = {}

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
        first_probe = probe_traj[0]
        probe_ids.popleft()
        print ("First probe being processed: ")
        pprint(vars(first_probe))

        params = Params(4.07, 20)
        print ("starting HMM with", curr_id)
        sequence, prob = MapMatchHMM(params, probe_traj, link_obj)
        print ("Sequence")
        candidate1 = []
        distFromRef1 = 0.0
        distFromLink1 = 0.0
        for i in range (0, len(sequence)):
            link1 = sequence[i]
            pprint(vars(link1))
            if (i+1 == len(sequence)):
                link2 = sequence[i]
            else:
                link2 = sequence[i+1]

            curr_probe = probe_traj[i]
            # Determine distance from reference node and distance from link
            if (i == 0):
                distFromRef1, distFromLink1, candidate1 = find_matched_params(curr_probe,link1) # current link
            if (i != len(sequence)-1):
                # Don't need to compute parameters for next link when already at the end
                distFromRef2, distFromLink2, candidate2 = find_matched_params(curr_probe,link2) # next link

            # Determine direction of travel of probe
            if (sequence[i].directionOfTravel != 'B'):
                curr_probe.direction = sequence[i].directionOfTravel
            else:
                if (i == len(sequence)-1):
                    curr_probe.direction = probe_traj[i-1].direction
                else:
                    curr_probe.direction = find_direction(link1, link2, candidate1, candidate2)

            # Derive the slope if link ID is the same
            if (link1.linkPVID == link2.linkPVID and i != len(sequence)-1):
                # Check if its a new link in the running list
                if link1.linkPVID not in link_slope_idx:
                    matching_links = [x for x in link_obj if (x.linkPVID == link1.linkPVID)]
                    matching_link = matching_links[0]
                    link_slopes.append(matching_link)
                    link_slope_idx[link1.linkPVID] = len(link_slopes)-1
                slope = derive_road_slope(curr_probe, probe_traj[i+1])
                previous_avg = link_slopes[link_slope_idx[link1.linkPVID]].derivedSlope
                n = link_slopes[link_slope_idx[link1.linkPVID]].numDerivedSlopes
                new_avg = (slope + (n * previous_avg))/(n+1)
                link_slopes[link_slope_idx[link1.linkPVID]].derivedSlope = new_avg
                link_slopes[link_slope_idx[link1.linkPVID]].numDerivedSlopes = n + 1

            curr_probe.linkPVID = link1.linkPVID
            curr_probe.distFromRef = distFromRef1
            curr_probe.distFromLink = distFromLink1
            # Update the current probe parameters
            candidate1 = candidate2
            distFromRef1 = distFromRef2
            distFromLink1 = distFromLink2
            probe_traj[i] = curr_probe
            row = []
            for attr, value in curr_probe.__dict__.items():
                if (attr == 'dateTime'):
                    time_str = time.strftime('%m/%d/%Y  %I:%M:%S %p', time.gmtime(value/1000.0))
                    row.append(time_str)
                else:
                    row.append(value)
            writer.writerow(row)
        print ("Probability", prob)

    print ("Processing derived slopes")
    # Process links with derived slopes
    for i in range(0, len(link_slopes)):
        link_slope = link_slopes[i]
        avg_slope_info = 0.0
        # Skip links without slope info
        if (link_slope.slopeInfo[0][0] == ''):
            continue
        # Compute average slopeInfo
        for j in range(0, len(link_slope.slopeInfo)):
            avg_slope_info = (float(link_slope.slopeInfo[j][1]) + j*avg_slope_info)/(j+1)
        # Compute absolute difference
        link_slope.absoluteDiff = math.fabs(avg_slope_info - link_slope.derivedSlope)
        # Compute relative error
        link_slope.percentDiff = link_slope.absoluteDiff/avg_slope_info
        row = []
        for attr, value in link_slope.__dict__.items():
            if (attr == 'linkPVID') or (attr == 'derivedSlope') or (attr == 'absoluteDiff') or (attr == 'percentDiff'):
                row.append(value)
        slope_writer.writerow(row)

    print ("End Time", datetime.now().time())
