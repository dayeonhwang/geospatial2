import numpy as np
import math
from helpers6 import *

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

def transition(probe1, probe2, candidate1, candidate2, link1, link2, params):
    # x_t,i: lat/lon pt on link1 nearest probe1 (shapeIndex1, minLat1, minLon1)
    # x_t+1,j: lat/lon pt on link2 nearest probe2 (shapeIndex2, minLat2, minLon2)
    # dist2: distance between x_t,j and x_t+1,j
    dist1 = compute_great_circle_distance(probe1.latitude,probe1.longitude,probe2.latitude,probe2.longitude)
    # shapeIndex1, minLat1, minLon1 = get_closest_shape(probe1,link1)
    # shapeIndex2, minLat2, minLon2 = get_closest_shape(probe2,link2)
    notFound = False
    print ("entire cand1", candidate1)
    print ("entire cand2", candidate2)
    print ("candidate", candidate1[1])
    print ("link 1 curvature", link1.curvatureInfo)
    shape_idx1 = candidate1[1]
    shape_idx2 = candidate2[1]
    print ("shape idx", shape_idx1)
    print ("shape info 1", link1.shapeInfo)
    print ("shape info 2", link1.shapeInfo)
    if (len(link1.shapeInfo) > 2):
        curve_dist1 = get_dist_curvature(link1.curvatureInfo, shape_idx1)
    elif (shape_idx1 == 0):
        curve_dist1 = 0
    else:
        curve_dist1 = link1.length

    if (len(link2.shapeInfo) > 2):
        curve_dist2 = get_dist_curvature(link2.curvatureInfo, shape_idx2)
    elif (shape_idx2 == 0):
        curve_dist2 = 0
    else:
        curve_dist2 = link2.length

    # Check if its along the same link
    if (link1.linkPVID == link2.linkPVID):
        dist2 = math.fabs(curve_dist1 - curve_dist2)
    else:
        if link1.nrefNodeID == link2.refNodeID:
            # traveling link1 -> link2
            dist2 = (link1.length - curve_dist1) + curve_dist2
        elif link1.refNodeID == link2.nrefNodeID:
            # traveling link2 -> link1
            dist2 = (link2.length - curve_dist2) + curve_dist1
        else:
            notFound = True
    if notFound:
        return 0
    else:
        dist = math.fabs(dist1 - dist2)
        beta = (1/math.log(2))*dist
        return (1/params.beta)*math.exp(-dist/params.beta)

def MapMatchHMM(params, trajectory, links):
    T = len(trajectory)
    # Create matrix that stores joint probability up to that state (states x time)
    #score = np.empty([0,T-1])
    score = []
    sequence = np.zeros(T-1)
    # Create index dictionary for road links in score matrix
    linkIdx = {}
    idx = 0 # running index for dictionary
    for t in range(T):
        R, R_links = get_candidate_nodes(trajectory[t], links)
        L, L_links= get_candidate_nodes(trajectory[t+1], links)
        if len(R) == 0:
            print ("no candidates")
        for l in range(0, len(L)):
            eProb = emission(L[l], trajectory[t+1], params)
            tProb = np.zeros(len(R))
            fProb = np.zeros(len(R))
            for r in range(0, len(R)):
                # Initialize scores
                if t == 0:
                    entry = np.zeros(T)
                    entry[t] = emission(R[r], trajectory[t], params)
                    score.append(entry)
                    linkIdx[R_links[r].linkPVID] = idx
                    idx += 1
                tProb[r] = transition(trajectory[t], trajectory[t+1], R[r], L[l], R_links[r], L_links[l], params)
                fProb[r] = tProb[r] * score[linkIdx[R_links[r].linkPVID]][t]
            maxr, maxri = np.max(fProb), np.argmax(fProb)
            sequence[t] = R_links[maxri].linkPVID
            finalScore = emission(L[l], trajectory[t+1], params)*maxr
            # Add state to the score matrix
            if (L_links[l].linkPVID not in linkIdx):
                entry = np.zeros(T)
                entry[t] = finalScore
                score.append(entry)
                linkIdx[L_links[l].linkPVID] = idx
                idx += 1
            else:
                score[linkIdx[L_links[l].linkPVID]][t] = finalScore
    return score, sequence

if __name__ == "__main__":
    probe_csv = open("Partition6467ProbePoints.csv", "r")
    link_csv = open("Partition6467LinkData.csv", "r")
    probe_reader = csv.reader(probe_csv)
    link_reader = csv.reader(link_csv)
    probe_rows = list(probe_reader)
    link_rows = list(link_reader)

    #candidate_nodes = get_candidate_nodes(probe_pt, link_rows)

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

    unsorted_probe_obj = probe_obj
    probe_obj.sort(key=lambda x: x.dateTime)
    print ("sorting dataTime")
    probe1 = probe_obj[0]
    traj = get_trajectory(probe_obj)
    #traj.sort(key=lambda x: x.dateTime)
    print ("got traj")
    probe1 = traj[probe1.sampleID]
    params = Params(4.07, 8.0)
    print ("starting HMM")
    hmm = MapMatchHMM(params, probe1, link_obj[0:100])
    print (hmm)