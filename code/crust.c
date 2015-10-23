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
	check the vertex lies on the convex hull
*/
int isOnConvexHull(tVertex vertex) {
	tEdge ptr_e;

	ptr_e = edges;
	do {
		if (ptr_e->endpts[0] == vertex || ptr_e->endpts[1] == vertex) {
			return TRUE;
		}
		ptr_e = ptr_e->next;
	} while (ptr_e != edges);

	return False;
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

/*
	For sample point s, find a pole

	vertexT *vertices, tVertex site
*/
void findPoleAntiPole() {
	tVertex site = NULL;
	tList site_voronoi_vertices = NULL;
	tVertex voronoi_vertex = NULL;
	tVertex pole_voronoi_vertex = NULL;
	tVertex anti_voronoi_vertex = NULL;
	vertexT *neighborVertex = NULL;
	facetT *neighbor, **neighborp;

	int isOnCH = 0;
	int vid = 0;
	double temp_norm = 0.0;
	double vector_norm = 0.0;

	site = vertices;
	isOnCH = isOnConvexHull(site);

	do {
		site_voronoi_vertices = site->vvlist;

		if (isOnCH) {	// If site s lies on the convex hull
			// n+ be the aveerage of the outer normals of the adjacent triangles;

			double facet_normal[3];
			
			//FOREACHneighbor_(site) {	// assign 'neighbor' to each neighbor in vertices->neighbors
			//	neighbor = (facetT*)site->neighbors;
			//	facet_normal[X] += neighbor->normal[X];
			//	facet_normal[Y] += neighbor->normal[Y];
			//	facet_normal[Z] += neighbor->normal[Z];
			//	vid++;
			////}
			//facet_normal[X] /= vid;
			//facet_normal[Y] /= vid;
			//facet_normal[Z] /= vid;
			
		}
		else {				// If site s does not lie on the cconvex hull,
			do {			// Find the farthest Voronoi vertex from s.
				voronoi_vertex = (tVertex)(site_voronoi_vertices->p);
				temp_norm = pointDist(site, voronoi_vertex);
				if (temp_norm > vector_norm) {
					pole_voronoi_vertex = voronoi_vertex;
					vector_norm = temp_norm;
				}
				site_voronoi_vertices = site_voronoi_vertices->next;
			} while (site_voronoi_vertices != site->vvlist);

			pole_voronoi_vertex->ispole = true;
			ADD(vertices, pole_voronoi_vertex);
			//findAntiPole(pole_voronoi_vertex, site);
		}
		site = site->next;
	} while (site != vertices);
	ADD(vertices, pole_voronoi_vertex);
}


void findAntiPole(tVertex pole, tVertex site) {
	tList site_voronoi_vertices;
	tVertex voronoi_vertex;
	tVertex antipole_voronoi_vertex;
	double temp_norm = 0.0;
	double vector_norm = 0.0;

	site_voronoi_vertices = site->vvlist;

	do {
		voronoi_vertex = (tVertex)(site_voronoi_vertices->p);
		temp_norm = innerProduct(site, pole, voronoi_vertex);
		if (temp_norm > vector_norm) {
			antipole_voronoi_vertex = voronoi_vertex;
			vector_norm = temp_norm;
		}
		site_voronoi_vertices = site_voronoi_vertices->next;
	} while (site_voronoi_vertices != site->vvlist);
	antipole_voronoi_vertex->ispole = TRUE;
	ADD(vertices, antipole_voronoi_vertex);
	
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
	//tsFace temp_face;
	double facetArea = 0.0;

	tVertex voronoi_vertex;
	
	int vvlnum = 0;
	tList vvl;
	int i = 0;
	int is_on_convexhull = 0;
	tVertex site;


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

		
		FOREACHvertex_(facet->vertices) {
			facetArea = qh_facetarea(facet);
			if (facet->upperdelaunay)	continue;
			site = all_v[qh_pointid(vertex->point)];

			if ((site->vvlist) == NULL) {	// avoid from checking redundantly
				int neighbor_size = 0;
				double avg_normal[3] = { 0, };

				FOREACHneighbor_(vertex) {
					if ((neighbor->upperdelaunay)) {	// is on convexhull
						is_on_convexhull = TRUE;
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

						if (neighbor->toporient) {
							neighbor_size++;
							avg_normal[X] += neighbor->normal[X];
							avg_normal[Y] += neighbor->normal[Y];
							avg_normal[Z] += neighbor->normal[Z];
						}

					}	//else
				}	//FOREACHneighbor_
				
				double *pole_vector;
				double antipole_vector[3] = { 0, };
				tVertex pole_voronoi_vertex = NULL;
				tVertex antipole_voronoi_vertex = NULL;
				tVertex temp_voronoi_vertex = NULL;
				tList site_voronoi_vertices;
				double temp_dist = 0;
				double max_dist = 0;

				if (is_on_convexhull) {
					pole_vector = (double *)calloc(3, sizeof(double));
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
						temp_dist = pointDist(site, voronoi_vertex);
						if (temp_dist > max_dist) {
							pole_voronoi_vertex = voronoi_vertex;
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
						voronoi_vertex = (tVertex)(site_voronoi_vertices->p);
						if (voronoi_vertex->ispole) {
							continue;
						}
						temp_dist = innerProduct(site, pole_voronoi_vertex, voronoi_vertex);
						if (temp_dist > max_dist) {
							antipole_voronoi_vertex = voronoi_vertex;
							max_dist = temp_dist;
						}
						site_voronoi_vertices = site_voronoi_vertices->next;
					} while (site_voronoi_vertices != site->vvlist);
				}
				else {	// the site is on CH
					do {
						voronoi_vertex = (tVertex)(site_voronoi_vertices->p);
						double dist1, dist2;
						tVertex temp_vertex = MakeTempVertex(pole_vector);
						temp_dist = innerProduct(site, temp_vertex, voronoi_vertex);
						if (temp_dist > max_dist) {
							antipole_voronoi_vertex = voronoi_vertex;
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
						ADD(vertices, antipole_voronoi_vertex);
					}
				}
				
			}	// if ((site->vvlist) == NULL)
		}	//FOREACHvertex_
	}	//FORALLfacets
				
				//site = vertices;
				//do {
				//	double max_dist = 0;
				//	double *pole_vector;
				//	double antipole_vector[3] = { 0, };
				//	int neighbor_size = 0;
				//	tVertex pole_voronoi_vertex = NULL;
				//	tVertex temp_voronoi_vertex = NULL;
				//	tList site_voronoi_vertices;
				//	double temp_dist = 0;
				//	double max_dist = 0;

				//	site_voronoi_vertices = site->vvlist;
				//	if (site_voronoi_vertices == NULL) {
				//		site = site->next;
				//		continue;
				//	}

				//	if (site->) {	// If site s lies on the convex hull
				//		// n+ be the aveerage of the outer normals of the adjacent triangles;


				//		double avg_normal[3] = { 0, };

				//		//FOREACHneighbor_(site) {	// assign 'neighbor' to each neighbor in vertices->neighbors
				//		//	neighbor = (facetT*)site->neighbors;
				//		//	facet_normal[X] += neighbor->normal[X];
				//		//	facet_normal[Y] += neighbor->normal[Y];
				//		//	facet_normal[Z] += neighbor->normal[Z];
				//		//	vid++;
				//		////}
				//		//facet_normal[X] /= vid;
				//		//facet_normal[Y] /= vid;
				//		//facet_normal[Z] /= vid;

				//	}
				//	else {				// If site s does not lie on the cconvex hull,
				//		do {			// Find the farthest Voronoi vertex from s.
				//			temp_voronoi_vertex = (tVertex)(site_voronoi_vertices->p);
				//			temp_dist = pointDist(site, voronoi_vertex);
				//			if (temp_dist > max_dist) {
				//				pole_voronoi_vertex = voronoi_vertex;
				//				max_dist = temp_dist;
				//			}
				//		} while (site_voronoi_vertices != site->vvlist);

				//		pole_voronoi_vertex->ispole = true;
				//		ADD(vertices, pole_voronoi_vertex);
				//		pole_vector = (double *)calloc(3, sizeof(double));
				//		for (i = 0; i < 3; i++) {
				//			pole_vector[i] = pole_voronoi_vertex->v[i];
				//		}
				//	}
				//	site = site->next;
				//} while (site != vertices);




	//site = vertices;
	//do {
	//	double max_dist = 0;
	//	double *pole_vector;
	//	double antipole_vector[3] = { 0, };
	//	int neighbor_size = 0;
	//	tVertex pole_voronoi_vertex = NULL;
	//	tVertex temp_voronoi_vertex = NULL;
	//	tList site_voronoi_vertices;
	//	double temp_dist = 0;
	//	double max_dist = 0;

	//	site_voronoi_vertices = site->vvlist;
	//	if (site_voronoi_vertices == NULL) {
	//		site = site->next;
	//		continue;
	//	}
	//	
	//	if (1) {	// If site s lies on the convex hull
	//		// n+ be the aveerage of the outer normals of the adjacent triangles;


	//		double avg_normal[3] = { 0, };

	//		//FOREACHneighbor_(site) {	// assign 'neighbor' to each neighbor in vertices->neighbors
	//		//	neighbor = (facetT*)site->neighbors;
	//		//	facet_normal[X] += neighbor->normal[X];
	//		//	facet_normal[Y] += neighbor->normal[Y];
	//		//	facet_normal[Z] += neighbor->normal[Z];
	//		//	vid++;
	//		////}
	//		//facet_normal[X] /= vid;
	//		//facet_normal[Y] /= vid;
	//		//facet_normal[Z] /= vid;

	//	}
	//	else {				// If site s does not lie on the cconvex hull,
	//		do {			// Find the farthest Voronoi vertex from s.
	//			temp_voronoi_vertex = (tVertex)(site_voronoi_vertices->p);
	//			temp_dist = pointDist(site, voronoi_vertex);
	//			if (temp_dist > max_dist) {
	//				pole_voronoi_vertex = voronoi_vertex;
	//				max_dist = temp_dist;
	//			}
	//		} while (site_voronoi_vertices != site->vvlist);

	//		pole_voronoi_vertex->ispole = true;
	//		ADD(vertices, pole_voronoi_vertex);
	//		pole_vector = (double *)calloc(3, sizeof(double));
	//		for (i = 0; i < 3; i++) {
	//			pole_vector[i] = pole_voronoi_vertex->v[i];
	//		}
	//	}
	//	site = site->next;
	//} while (site != vertices);
	//

	// find a pole

		/*if (is_on_convexhull) {
		pole_vector = (double *)calloc(3, sizeof(double));
		avg_normal[X] /= neighbor_size;
		avg_normal[Y] /= neighbor_size;
		avg_normal[Z] /= neighbor_size;
		for (i = 0; i < 3; i++) {
		pole_vector[i] = avg_normal[i];
		}
		if (pole_voronoi_vertex) {
		pole_vector = (double *)calloc(3, sizeof(double));
		pole_voronoi_vertex->ispole = true;
		ADD(vertices, pole_voronoi_vertex);
		for (i = 0; i < 3; i++) {
		pole_vector[i] = pole_voronoi_vertex->v[i];
		}
		}
		vvl = MakeNullVoronoiVertexList();
		vvl->p = (void *)vertex;
		ADD(site->vvlist, vvl);
		}*/
/*
	double max_dist = 0;
	double avg_normal[3] = { 0, };
	double *pole_vector;
	double antipole_vector[3] = { 0, };
	int neighbor_size = 0;
	tVertex pole_voronoi_vertex = NULL;
	tVertex temp_voronoi_vertex = NULL;
	if (neighbor->upperdelaunay){
		is_on_convexhull = TRUE;
	}
	else {
		voronoi_vertex = MakeNullVoronoiVertex();
		voronoi_vertex = neighbor->center;

		if (neighbor->toporient) {
			neighbor_size++;
			avg_normal[X] += neighbor->normal[X];
			avg_normal[Y] += neighbor->normal[Y];
			avg_normal[Z] += neighbor->normal[Z];
		}

		double temp_dist = pointDist(vertex, voronoi_vertex);

		if (temp_dist > max_dist) {
			pole_voronoi_vertex = voronoi_vertex;
			pole_voronoi_vertex->ispole = true;
			max_dist = temp_dist;
		}
	}

	}*/ //FOREACHneighbor_
	//	if (is_on_convexhull) {
	//		pole_vector = (double *)calloc(3, sizeof(double));
	//		avg_normal[X] /= neighbor_size;
	//		avg_normal[Y] /= neighbor_size;
	//		avg_normal[Z] /= neighbor_size;
	//		for (i = 0; i < 3; i++) {
	//			pole_vector[i] = avg_normal[i];
	//		}
	//	}
	//if (pole_voronoi_vertex) {
	//	pole_vector = (double *)calloc(3, sizeof(double));
	//	pole_voronoi_vertex->ispole = true;
	//	ADD(vertices, pole_voronoi_vertex);
	//	for (i = 0; i < 3; i++) {
	//		pole_vector[i] = pole_voronoi_vertex->v[i];
	//	}
	//}

	//// fine a antipole
	//if (pole_vector) {
	//	/*temp_voronoi_vertex = voronoi_vertices;
	//	do {
	//	if (!(temp_voronoi_vertex->ispole))
	//	temp_dist = innerProduct(vertex->, pole, voronoi_vertex);
	//	if (temp_norm > vector_norm) {
	//	antipole_voronoi_vertex = voronoi_vertex;
	//	vector_norm = temp_norm;
	//	}
	//	site_voronoi_vertices = site_voronoi_vertices->next;
	//	} while (temp_voronoi_vertex != voronoi_vertices);
	//	antipole_voronoi_vertex->ispole = TRUE;
	//	ADD(vertices, antipole_voronoi_vertex);*/
	//}

				//FOREACHneighbor_(site) {	// assign 'neighbor' to each neighbor in vertices->neighbors
				//	neighbor = (facetT*)site->neighbors;
				//	facet_normal[X] += neighbor->normal[X];
				//	facet_normal[Y] += neighbor->normal[Y];
				//	facet_normal[Z] += neighbor->normal[Z];
				//	vid++;
				////}
				//facet_normal[X] /= vid;
				//facet_normal[Y] /= vid;
				//facet_normal[Z] /= vid;
	
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
