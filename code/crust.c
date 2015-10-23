#include "delauany.h"
#include "crust.h"
#include "shape.h"

/*-------------------------------------------------------------------*/
//headers from qhull
#include "qhull.h"
#include "poly.h"
#include "qset.h"
#include "geom.h"
#include "alpha-shape.h"

#define TRUE 1
#define FALSE 0
#define MAX_DIST 10000
#define IS_ON_CH 3

/*-------------------------------------------------------------------*/
//defined in shape.h/.c
extern tVertex vertices;
extern tEdge edges;
extern tFace faces;
extern tTetra tetras;

tVertex voronoi_vertices = NULL;

tVertex MakeTempVertex(double *vector) {
	tVertex v;
	NEW(v, tsVertex);
	v->v[0] = vector[0];
	v->v[1] = vector[1];
	v->v[2] = vector[2];
	v->vvlist = NULL;
	v->ispole = false;
}

tVertex	MakeNullVoronoiVertex(void)
{
	tVertex  v;
	NEW(v, tsVertex);
	ADD(voronoi_vertices, v);
	v->vvlist = NULL;
	v->ispole = false;

	return v;
}

tList MakeNullVoronoiVertexList(void) {
	tList  vvl;
	NEW(vvl, tsList);
	vvl->next = NULL;
	vvl->prev = NULL;
	vvl->p = NULL;

	return vvl;
}

/*
get angle between two vertex
*/
double getAngle(tVertex v1, tVertex v2) {
	double pt1[3] = { v1->v[0], v1->v[1], v1->v[2] };
	double pt2[3] = { v2->v[0], v2->v[1], v2->v[2] };

	return qh_getangle(pt1, pt2);
}

/*
compute the inner product of the vectors, (n+, n-)
*/
double innerProduct(tVertex site, tVertex pole, tVertex v) {
	double dist1, dist2;
	double angle = 0.0;
	double result = 0.0;

	dist1 = pointDist(site, pole);
	dist2 = pointDist(site, v);
	angle = getAngle(pole, v);
	result = dist1 * dist2* cos(angle);

	return result;
}

void findPoleAntiPole(int vsize) {
	tVertex site;
	double *pole_vector;
	double avg_normal[3] = { 0, };
	facetT *neighbor, **neighborp;
	tVertex pole_voronoi_vertex = NULL;
	tVertex antipole_voronoi_vertex = NULL;
	vertexT *temp_voronoi_vertexT = NULL;
	tVertex temp_voronoi_vertex;
	tList site_voronoi_vertices;
	double temp_dist = 0;
	double max_dist = 0;
	int neighbor_size = 0;
	int i;
	//int count = 1;

	site = vertices;
	// for all vertices
	do {
		site_voronoi_vertices = site->vvlist;
		if (!site_voronoi_vertices)	{
			site = site->next;
			continue;
		}
		temp_voronoi_vertexT = (vertexT*)site_voronoi_vertices->p;

		if (!temp_voronoi_vertexT) {	// lies on the CH: compute the average of the outer nomals of the adjacents.
			pole_vector = (double *)calloc(3, sizeof(double));
			neighbor_size = 0;
			FOREACHneighbor_(temp_voronoi_vertexT) {
				neighbor_size++;
				avg_normal[X] += neighbor->normal[X];
				avg_normal[Y] += neighbor->normal[Y];
				avg_normal[Z] += neighbor->normal[Z];
			}

			avg_normal[X] /= neighbor_size;
			avg_normal[Y] /= neighbor_size;
			avg_normal[Z] /= neighbor_size;
			for (i = 0; i < 3; i++) {
				pole_vector[i] = avg_normal[i];
			}
		}
		else {
			site_voronoi_vertices = site->vvlist;
			do {			// Find the farthest Voronoi vertex from s.
				temp_voronoi_vertex = (tVertex)(site_voronoi_vertices->p);
				temp_dist = pointDist(site, temp_voronoi_vertex);
				if (temp_dist > max_dist) {
					pole_voronoi_vertex = temp_voronoi_vertex;
					max_dist = temp_dist;
					pole_voronoi_vertex->ispole = true;
				}
			} while (site_voronoi_vertices != site->vvlist);

			if (pole_voronoi_vertex){
				if (max_dist > MAX_DIST) {
					pole_voronoi_vertex = NULL;
				}
				else {
					pole_voronoi_vertex->ispole = TRUE;
					pole_voronoi_vertex->vnum = vsize++;

					ADD(vertices, pole_voronoi_vertex);
					pole_vector = (double *)calloc(3, sizeof(double));
					for (i = 0; i < 3; i++) {
						pole_vector[i] = pole_voronoi_vertex->v[i];
					}
				}
			}

		}

		max_dist = 0;
		temp_dist = 0;

		// find a antipole
		site_voronoi_vertices = site->vvlist;
		assert(pole_vector != NULL);

		if (pole_voronoi_vertex != NULL) {
			do {
				temp_voronoi_vertex = (tVertex)(site_voronoi_vertices->p);
				if (temp_voronoi_vertex->ispole) {
					continue;
				}
				temp_dist = innerProduct(site, pole_voronoi_vertex, temp_voronoi_vertex);
				if (temp_dist > max_dist) {
					antipole_voronoi_vertex = temp_voronoi_vertex;
					max_dist = temp_dist;
				}
				site_voronoi_vertices = site_voronoi_vertices->next;
			} while (site_voronoi_vertices != site->vvlist);
		}
		else {	// the site is on CH
			do {
				temp_voronoi_vertex = (tVertex)(site_voronoi_vertices->p);
				double dist1, dist2;
				tVertex temp_vertex = MakeTempVertex(pole_vector);
				temp_dist = innerProduct(site, temp_vertex, temp_voronoi_vertex);
				if (temp_dist > max_dist) {
					antipole_voronoi_vertex = temp_voronoi_vertex;
					max_dist = temp_dist;
				}
				site_voronoi_vertices = site_voronoi_vertices->next;
			} while (site_voronoi_vertices != site->vvlist);
		}
		if (antipole_voronoi_vertex){
			if (max_dist > MAX_DIST) {
				antipole_voronoi_vertex = NULL;
			}
			else {
				antipole_voronoi_vertex->ispole = TRUE;
				antipole_voronoi_vertex->vnum = vsize++;

				ADD(vertices, antipole_voronoi_vertex);
			}
		}
		max_dist = 0;
		pole_vector[0] = 0;
		pole_vector[1] = 0;
		pole_vector[2] = 0;
		pole_voronoi_vertex = NULL;
		antipole_voronoi_vertex = NULL;
		site = site->next;
	} while (site != vertices);

}




/*-------------------------------------------------------------------*/
void	Crust(void)
{
	// implement 3d crust alogorithm here
	// put results (vertices, edges, faces, and tetrahedra)
	// in the global structures above
	// see main.c for detail on how qhull does it

	tVertex  ptr_v;
	tVertex * all_v = NULL;
	int vsize = 0;
	int id = 0;

	/* global varibles for qhull */

	//	Qhull options
	//	Qz   - add point-at-infinity to Voronoi diagram
	static char * options = (char *)"qhull Qbb QJ Pp";
	int curlong, totlong;
	coordT * pt = NULL;
	facetT *facet = NULL;
	vertexT *vertex = NULL;
	vertexT **vertexp = NULL;
	facetT *neighbor, **neighborp;
	int vid = 0;

	// for crust
	tTetra tetra;
	double facetArea = 0.0;

	tVertex voronoi_vertex;
	tVertex site;
	int vvlnum = 0;
	tList vvl;

	//count number of points
	ptr_v = vertices;
	do {
		vsize++;
		ptr_v = ptr_v->next;
	} while (ptr_v != vertices);

	//allocate memory
	pt = (coordT*)calloc(vsize * 4, sizeof(coordT)); //each point will have four coord
	all_v = (tVertex*)calloc(vsize, sizeof(tVertex));
	assert(pt && all_v);

	//copy points
	ptr_v = vertices;
	do {
		pt[id++] = ptr_v->v[0];
		pt[id++] = ptr_v->v[1];
		pt[id++] = ptr_v->v[2];
		pt[id++] = (ptr_v->v[0] * ptr_v->v[0]) + (ptr_v->v[1] * ptr_v->v[1]) + (ptr_v->v[2] * ptr_v->v[2]);
		all_v[ptr_v->vnum] = ptr_v;
		ptr_v = ptr_v->next;
	} while (ptr_v != vertices);

	//using qhull
	qh_init_A(stdin, stdout, stderr, 0, NULL);

	//qh DELAUNAY= True;     /* 'd'   */
	//qh SCALElast= True;    /* 'Qbb' */
	//qh KEEPcoplanar= True; /* 'Qc', to keep coplanars in 'p' */

	qh_initflags(options);
	qh_init_B(pt, vsize, 4, false);
	qh_qhull();
	qh_check_output();

	//compute Voronoi centers for all facets
	//includes upperDelaunay facets if qh.UPPERdelaunay ('Qu')
	qh_setvoronoi_all();

	FORALLfacets
	{
		if (facet->center == NULL)
		facet->center = qh_facetcenter(facet->vertices);
		facetArea = qh_facetarea(facet);
		FOREACHvertex_(facet->vertices) {
			//if (facet->upperdelaunay)	continue;
			if (facet->upperdelaunay || !facetArea)	continue;
			site = all_v[qh_pointid(vertex->point)];

			if ((site->vvlist) == NULL) {	// avoid from checking redundantly

				FOREACHneighbor_(vertex) {
					if ((neighbor->upperdelaunay)) {
						//is_on_convexhull = TRUE;
					}
					else {
						voronoi_vertex = MakeNullVoronoiVertex();
						voronoi_vertex->v[X] = neighbor->center[X];
						voronoi_vertex->v[Y] = neighbor->center[Y];
						voronoi_vertex->v[Z] = neighbor->center[Z];
						voronoi_vertex->vnum = vvlnum++;
						vvl = MakeNullVoronoiVertexList();
						vvl->p = (void*)voronoi_vertex;
						ADD(site->vvlist, vvl);

					}	//else
				}	//FOREACHneighbor_

			}	// if ((site->vvlist) == NULL)
		}	//FOREACHvertex_
	}	//FORALLfacets

		// find a pole
	findPoleAntiPole(vsize);
	printf("vsize = %d\n", vsize);

	Delaunay();

	//not used
	free(pt);
	free(all_v);
	pt = NULL;
	all_v = NULL;

	//free mem
	qh_freeqhull(!qh_ALL);
	qh_memfreeshort(&curlong, &totlong);


}
