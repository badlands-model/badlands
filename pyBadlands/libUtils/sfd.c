//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~//
//                                                                                   //
//  This file forms part of the Badlands surface processes modelling application.    //
//                                                                                   //
//  For full license and copyright information, please refer to the LICENSE.md file  //
//  located at the project root, or contact the authors.                             //
//                                                                                   //
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~//

// This module computes the Single Flow Direction for any given surface.

#include <stdio.h>

#define MAX_NEIGHBOURS 20

void directions(double pyElev[], double pyZ[], int pyNgbs[][MAX_NEIGHBOURS], double pyEdge[][MAX_NEIGHBOURS],
    double pyDist[][MAX_NEIGHBOURS], int pyGIDs[], double sealimit, int pyBase[],
    int pyRcv[], double pyMaxh[], double pyMaxDep[],
    int pylocalNb, int pyglobalNb)
{
    int i;

    for (i = 0; i < pyglobalNb; i++) {
        pyBase[i] = -1;
        pyRcv[i] = -1;
        pyMaxh[i] = 1.e6;
        pyMaxDep[i] = 0.;
    }

    int k;
    for (k = 0; k < pylocalNb; k++) {
        int gid = pyGIDs[k];
        int lowestID = gid;
        double diffH = 1.e6;
        double diffD = 0.;
        int p;

        for (p = 0; p < MAX_NEIGHBOURS; p++) {
            int ngbid = pyNgbs[gid][p];
            if (ngbid < 0) {
                break;
            }

            if (pyElev[ngbid] < pyElev[lowestID]) {
                lowestID = ngbid;
            }

            double dh = pyZ[ngbid] - pyZ[gid];

            if (dh >= 0. && dh < diffH) {
                diffH = dh;
            }

            if (dh > diffD) {
                diffD = dh;
            }
        }

        pyRcv[gid] = lowestID;
        if (pyZ[gid] < sealimit) {
            pyRcv[gid] = gid;
        }

        if (gid == pyRcv[gid] && pyZ[gid] + diffD > sealimit) {
            pyBase[gid] = gid;
        }

        if (diffH > 9.99e5) {
            diffH = 0.;
        }

        pyMaxh[gid] = diffH;
        pyMaxDep[gid] = diffD;
    }
}

void directions_base(double pyZ[], int pyNgbs[][MAX_NEIGHBOURS], double pyEdge[][MAX_NEIGHBOURS],
    double pyDist[][MAX_NEIGHBOURS], int pyGIDs[], double sealimit, int pyBase[],
    int pyRcv[], int pylocalNb, int pyglobalNb)
{
    int i;

    for (i = 0; i < pyglobalNb; i++) {
        pyBase[i] = -1;
        pyRcv[i] = -1;
    }

    int k;
    for (k = 0; k < pylocalNb; k++) {
        int gid = pyGIDs[k];
        int lowestID = gid;
        double diffD = 0.;
        int p;

        for (p = 0; p < MAX_NEIGHBOURS; p++) {
            int ngbid = pyNgbs[gid][p];
            if (ngbid < 0) {
                break;
            }

            if (pyZ[ngbid] < pyZ[lowestID]) {
                lowestID = ngbid;
            }

            double dh = pyZ[ngbid] - pyZ[gid];

            if (dh > diffD) {
                diffD = dh;
            }
        }

        pyRcv[gid] = lowestID;
        if (pyZ[gid] < sealimit) {
            pyRcv[gid] = gid;
        }

        if (gid == pyRcv[gid] && pyZ[gid] + diffD > sealimit) {
            pyBase[gid] = gid;
        }
    }
}

void diffusion(double pyZ[], int pyNgbs[][MAX_NEIGHBOURS], double pyEdge[][MAX_NEIGHBOURS],
    double pyDist[][MAX_NEIGHBOURS], int pyGIDs[], double pyDiff[], int pylocalNb, int pyglobalNb)
{
    int i;

    for (i = 0; i < pyglobalNb; i++) {
        pyDiff[i] = 0.;
    }

    int k;
    for (k = 0; k < pylocalNb; k++) {
        int gid = pyGIDs[k];
        int p;

        for (p = 0; p < MAX_NEIGHBOURS; p++) {
            int ngbid = pyNgbs[gid][p];
            if (ngbid < 0) {
                break;
            }

            pyDiff[gid] += pyEdge[gid][p] * (pyZ[ngbid] - pyZ[gid]) / pyDist[gid][p];
        }
    }
}
