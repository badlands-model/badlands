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

void dirview(double pyElev[], double pyZ[], int pyNgbs[][MAX_NEIGHBOURS], double pyEdge[][MAX_NEIGHBOURS],
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
        int p;

        for (p = 0; p < MAX_NEIGHBOURS; p++) {
            int ngbid = pyNgbs[gid][p];
            if (ngbid < 0) {
                break;
            }

            if (pyElev[ngbid] < pyElev[lowestID]) {
                lowestID = ngbid;
            }
        }

        pyRcv[gid] = lowestID;
        if (pyZ[gid] < sealimit) {
            pyRcv[gid] = gid;
        }

        if (gid == pyRcv[gid]) {
            pyBase[gid] = gid;
        }
    }
}

void directions(double pyElev[], double pyZ[], int pyNgbs[][MAX_NEIGHBOURS], double pyEdge[][MAX_NEIGHBOURS],
    double pyDist[][MAX_NEIGHBOURS], int pyGIDs[], int pyBase[],
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

        if (gid == pyRcv[gid]) {
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
    double pyDist[][MAX_NEIGHBOURS], int pyGIDs[], int pyBase[],
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
        if (gid == pyRcv[gid]) {
            pyBase[gid] = gid;
        }
    }
}

void diffusion(double pyZ[], int pyBord[], int pyNgbs[][MAX_NEIGHBOURS], double pyEdge[][MAX_NEIGHBOURS],
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
        if (pyBord[gid]>0) {
          for (p = 0; p < MAX_NEIGHBOURS; p++) {
              int ngbid = pyNgbs[gid][p];
              if (ngbid < 0) {
                  break;
              }
              if (pyBord[ngbid]>0 && pyDist[gid][p] > 0.){
                pyDiff[gid] += pyEdge[gid][p] * (pyZ[ngbid] - pyZ[gid]) / pyDist[gid][p];
              }
              if (pyBord[ngbid]<1){
                if (pyZ[ngbid] < pyZ[gid] && pyDist[gid][p] > 0.){
                  pyDiff[gid] += pyEdge[gid][p] * (pyZ[ngbid] - pyZ[gid]) / pyDist[gid][p];
                }
              }
          }
        }
    }
}

void diffusionero(double pyZ[], int pyBord[], int pyNgbs[][MAX_NEIGHBOURS], double pyEdge[][MAX_NEIGHBOURS],
    double pyDist[][MAX_NEIGHBOURS], int pyGIDs[], double pyEro[], int pylocalNb, int pyglobalNb)
{
    int i;

    for (i = 0; i < pyglobalNb; i++) {
        pyEro[i] = 0.;
    }

    int k;
    for (k = 0; k < pylocalNb; k++) {
        int gid = pyGIDs[k];
        int p;
        if (pyBord[gid]>0) {
          for (p = 0; p < MAX_NEIGHBOURS; p++) {
              int ngbid = pyNgbs[gid][p];
              if (ngbid < 0) {
                  break;
              }
              if (pyBord[ngbid]>0 && pyZ[gid] > pyZ[ngbid] && pyDist[gid][p] > 0.){
                pyEro[gid] += pyEdge[gid][p] * (pyZ[ngbid] - pyZ[gid]) / pyDist[gid][p];
              }
              if (pyBord[ngbid]<1){
                if (pyZ[ngbid] < pyZ[gid] && pyDist[gid][p] > 0. ){
                  pyEro[gid] += pyEdge[gid][p] * (pyZ[ngbid] - pyZ[gid]) / pyDist[gid][p];
                }
              }
          }
        }
    }
}

void diffusionmarine(double pyZ[], int pyBord[], int pyDep[], int pyNgbs[][MAX_NEIGHBOURS], double pyEdge[][MAX_NEIGHBOURS],
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
        if (pyBord[gid]>0.) {
          for (p = 0; p < MAX_NEIGHBOURS; p++) {
              int ngbid = pyNgbs[gid][p];
              if (ngbid < 0) {
                  break;
              }
              if (pyBord[ngbid]>0){
                if(pyDep[gid] > 0 && pyZ[gid] > pyZ[ngbid] && pyDist[gid][p] > 0.) {
                  pyDiff[gid] += pyEdge[gid][p] * (pyZ[ngbid] - pyZ[gid]) / pyDist[gid][p];
                }
                if(pyDep[ngbid] > 0 && pyZ[gid] < pyZ[ngbid] && pyDist[gid][p] > 0.) {
                  pyDiff[gid] += pyEdge[gid][p] * (pyZ[ngbid] - pyZ[gid]) / pyDist[gid][p];
                }
              }
              if (pyBord[ngbid]<1){
                if (pyZ[ngbid] < pyZ[gid] && pyDep[gid] > 0 && pyDist[gid][p] > 0.){
                  pyDiff[gid] += pyEdge[gid][p] * (pyZ[ngbid] - pyZ[gid]) / pyDist[gid][p];
                }
              }
          }
        }
    }
}
