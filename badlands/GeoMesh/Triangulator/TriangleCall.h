/* ============================================================================
 *  	This file calls  Triangle to create the top surface for SPModel.
 *  	Jonathan Richard Shewchuk
 *  	2360 Woolsey #H
 *  	Berkeley, California  94705-1927
 *  	jrs@cs.berkeley.edu
 * ============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#define NUL '\0'
#define VOID int
#define REAL double
#define FILENAMESIZE 2048

// TYPES

typedef REAL **triangle;
typedef REAL **subseg;
typedef REAL *vertex;

// STRUCTURES

struct otri {
  triangle *tri;
  int orient;
};

struct badsubseg {
  subseg encsubseg;
  vertex subsegorg, subsegdest;
};

struct badtriang {
  triangle poortri;
  REAL key;
  vertex triangorg, triangdest, triangapex;
  struct badtriang *nexttriang;
};

struct flipstacker {
  triangle flippedtri;
  struct flipstacker *prevflip;
};

struct memorypool {
  VOID **firstblock, **nowblock;
  VOID *nextitem;
  VOID *deaditemstack;
  VOID **pathblock;
  VOID *pathitem;
  int alignbytes;
  int itembytes;
  int itemsperblock;
  int itemsfirstblock;
  long items, maxitems;
  int unallocateditems;
  int pathitemsleft;
};

struct mesh {
  struct memorypool triangles;
  struct memorypool subsegs;
  struct memorypool vertices;
  struct memorypool viri;
  struct memorypool badsubsegs;
  struct memorypool badtriangles;
  struct memorypool flipstackers;
  struct memorypool splaynodes;
  struct badtriang *queuefront[4096];
  struct badtriang *queuetail[4096];
  int nextnonemptyq[4096];
  int firstnonemptyq;
  struct flipstacker *lastflip;
  REAL xmin, xmax, ymin, ymax;
  REAL xminextreme;
  int invertices;
  int inelements;
  int insegments;
  int holes;
  int regions;
  int undeads;
  long edges;
  int mesh_dim;
  int nextras;
  int eextras;
  long hullsize;
  int steinerleft;
  int vertexmarkindex;
  int vertex2triindex;
  int highorderindex;
  int elemattribindex;
  int areaboundindex;
  int checksegments;
  int checkquality;
  int readnodefile;
  long samples;
  long incircledcount;
  long counterclockcount;
  long oriented3dcount;
  long hyperbolacount;
  long circumcentercount;
  long circletopcount;
  vertex infvertex1, infvertex2, infvertex3;
  triangle *dummytri;
  triangle *dummytribase;
  subseg *dummysub;
  subseg *dummysubbase;
  struct otri recenttri;
};

struct behavior {
  int poly, refine, quality, vararea, fixedarea, usertest;
  int regionattrib, convex, weighted, jettison;
  int firstnumber;
  int edgesout, voronoi, neighbors, geomview;
  int nobound, nopolywritten, nonodewritten, noelewritten, noiterationnum;
  int noholes, noexact, conformdel;
  int incremental, sweepline, dwyer;
  int splitseg;
  int docheck;
  int quiet, verbose;
  int usesegments;
  int order;
  int nobisect;
  int steiner;
  REAL minangle, goodangle, offconstant;
  REAL maxarea;
  char innodefilename[FILENAMESIZE];
  char inelefilename[FILENAMESIZE];
  char inpolyfilename[FILENAMESIZE];
  char areafilename[FILENAMESIZE];
  char outnodefilename[FILENAMESIZE];
  char outelefilename[FILENAMESIZE];
  char outpolyfilename[FILENAMESIZE];
  char edgefilename[FILENAMESIZE];
  char vnodefilename[FILENAMESIZE];
  char vedgefilename[FILENAMESIZE];
  char neighborfilename[FILENAMESIZE];
  char offfilename[FILENAMESIZE];
};

// FUNCTIONS

char *F2C_trim( char *str ) ;

void trianglegen_( char copt[], char fsurf[] ) ;

void triangulation_fct( int nbcmd, char *argument[] ) ;




