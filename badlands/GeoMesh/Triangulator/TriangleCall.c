/* ============================================================================
 * 		TriangleCall.c is based on J. R. Shewchuk 
 *
 *  	This file calls  Triangle to create the top surface for SPModel.
 *  	Jonathan Richard Shewchuk
 *  	2360 Woolsey #H
 *  	Berkeley, California  94705-1927
 *  	jrs@cs.berkeley.edu
 * ============================================================================*/

#include "TriangleCall.h"
#include "triangle.h"

#include <stdio.h>

#ifndef NULL
#define NULL   ((void *) 0)
#endif

/* ============================================================================
 * Function F2C_trim()
 * When passing string from fortran to C some end character are mismatching.
 * This function ensures the consistency between string for the 2 languages
 * ============================================================================*/
char *F2C_trim( char *str )
{
  char *ibuf = str, *obuf = str;
  int i = 0, cnt = 0;
  if (str) {
    for (ibuf = str; *ibuf && isspace(*ibuf); ++ibuf) ;
    if (str != ibuf) memmove(str, ibuf, ibuf - str);
    while (*ibuf) {
      if (isspace(*ibuf) && cnt)
	ibuf++;
      else {
	if (!isspace(*ibuf))
	  cnt = 0;
	else {
	  *ibuf = ' ';
	  cnt = 1;
	}
	obuf[i++] = *ibuf++;
      }
    }
    obuf[i] = NUL;
    while (--i >= 0) {
      if (!isspace(obuf[i])) break;
    }
    obuf[++i] = NUL;
  }

  return str;
}

/* ============================================================================
 * Function trianglegen_()
 * This function ensures the interface between SPModel and Triangle. It takes the input nodes
 * from SPModel and passes it over to Triangle which creates the appropriate TIN surface.
 * ============================================================================*/
void trianglegen_( char copt[], char fsurf[] )
{
  int chosen ;
  int nbcmd = 4 ;
  char *argument[4] ;
  char *prog = "triangle" ;
  argument[0] = prog ;
  argument[1] = "-QpevD" ;
  argument[2] = F2C_trim(copt) ;
  argument[3] = F2C_trim(fsurf) ;
  // Call adapted function for running Triangle
  triangulation_fct( nbcmd, argument );

}

/* ============================================================================
 * Function triangulation_fct()
 * This function is derived from the triangulate() function from Triangle. We use this modified
 * version because we want to pass the nodes file to Triangle and we want to have the elements
 * from Triangle directly output in a file. Instead of having to create out own reporting function
 * it was making more sense to let Triangle do it and just adapting the calling function..
 * ============================================================================*/
void triangulation_fct(int nbcmd, char *argument[])
{
  struct mesh m;
  struct behavior b;
  REAL *holearray;
  REAL *regionarray;
  FILE *polyfile;
  triangleinit(&m);
  parsecommandline(nbcmd, argument, &b);

  m.steinerleft = b.steiner;

  readnodes(&m, &b, b.innodefilename, b.inpolyfilename, &polyfile);

  if (b.refine) {
    m.hullsize = reconstruct(&m, &b, b.inelefilename, b.areafilename,
			     b.inpolyfilename, polyfile);
  } else {
    m.hullsize = delaunay(&m, &b);
  }

  m.infvertex1 = (vertex) NULL;
  m.infvertex2 = (vertex) NULL;
  m.infvertex3 = (vertex) NULL;

  if (b.usesegments) {
    m.checksegments = 1;
    if (!b.refine) {
      formskeleton(&m, &b, polyfile, b.inpolyfilename);
    }
  }

  if (b.poly && (m.triangles.items > 0)) {
    readholes(&m, &b, polyfile, b.inpolyfilename, &holearray, &m.holes,
	      &regionarray, &m.regions);
    if (!b.refine) {
      carveholes(&m, &b, holearray, m.holes, regionarray, m.regions);
    }
  } else {
    m.holes = 0;
    m.regions = 0;
  }

  if (b.quality && (m.triangles.items > 0)) enforcequality(&m, &b);
  m.edges = (3l * m.triangles.items + m.hullsize) / 2l;
  if (b.order > 1) highorder(&m, &b);
  if (b.nonodewritten || (b.noiterationnum && m.readnodefile)) {
    numbernodes(&m, &b);
  } else {
    writenodes(&m, &b, b.outnodefilename, nbcmd, argument);
  }
  if (b.noelewritten) {;}
  else {
    writeelements(&m, &b, b.outelefilename, nbcmd, argument);
  }
  if (b.poly || b.convex) {
    if (b.nopolywritten || b.noiterationnum) {;}
    else {
      writepoly(&m, &b, b.outpolyfilename, holearray, m.holes, regionarray,
		m.regions, nbcmd, argument);
    }
  }
  if (m.regions > 0) trifree((VOID *) regionarray);
  if (m.holes > 0) trifree((VOID *) holearray);
  if (b.geomview) writeoff(&m, &b, b.offfilename, nbcmd, argument);
  if (b.edgesout) writeedges(&m, &b, b.edgefilename, nbcmd, argument);
  if (b.voronoi) writevoronoi(&m, &b, b.vnodefilename, b.vedgefilename, nbcmd, argument);
  if (b.neighbors) writeneighbors(&m, &b, b.neighborfilename, nbcmd, argument);
  if (!b.quiet) statistics(&m, &b);
  if (b.docheck) {
    checkmesh(&m, &b);
    checkdelaunay(&m, &b);
  }
  triangledeinit(&m, &b);

}
