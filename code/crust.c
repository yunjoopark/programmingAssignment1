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
/*-------------------------------------------------------------------*/
//defined in shape.h/.c
extern tVertex vertices;
extern tEdge edges;
extern tFace faces;
extern tTetra tetras;

tVertex voronoi_vertices = NULL;

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
void findPole(int vsize) {
	tVertex site;
	tList site_voronoi_vertices;
	tVertex voronoi_vertex;
	tVertex pole_voronoi_vertex;
	facetT *neighbor, **neighborp;
	vertexT *neighborVertex = NULL;

	int isOnCH = 0;
	double temp_norm = 0.0;
	double vector_norm = 0.0;

	site = vertices;
	isOnCH = isOnConvexHull(site);

	//printf("is convex hull ? %d\n", isOnCH);

	do {
		site_voronoi_vertices = site->vvlist;

		if (isOnCH) {	// If site s lies on the convex hull
			// n+ be the aveerage of the outer normals of the adjacent triangles;

			//FOREACHneighbor_(vertices) {	// assign 'neighbor' to each neighbor in vertices->neighbors
			//	neighbor = (facetT*)vertices->neighbors;
			//	neighbor->normal;

			//}

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
		}

		site = site->next;
	} while (site != vertices);

}

void findAntiPole(tVertex site) {
	tList site_voronoi_vertices;
	tVertex voronoi_vertex;
	tVertex pole_voronoi_vertex;
	tVertex antipole_voronoi_vertex;
	double temp_norm = 0.0;
	double vector_norm = 0.0;
	
	// seek a pole
	site_voronoi_vertices = site->vvlist;
	do {
		pole_voronoi_vertex = (tVertex)(site_voronoi_vertices->p);
		if (!pole_voronoi_vertex->ispole){
			site_voronoi_vertices = site_voronoi_vertices->next;
			continue;
		}
		else {
			break;
		}
	} while (site_voronoi_vertices != site->vvlist);

	site_voronoi_vertices = site->vvlist;
	do {
		voronoi_vertex = (tVertex)(site_voronoi_vertices->p);
		temp_norm = innerProduct(site, pole_voronoi_vertex, voronoi_vertex);
		if (temp_norm > vector_norm) {
			antipole_voronoi_vertex = voronoi_vertex;
			vector_norm = temp_norm;
		}
		site_voronoi_vertices = site_voronoi_vertices->next;
	} while (site_voronoi_vertices != site->vvlist);
}

/*-------------------------------------------------------------------*/
void	Crust( void )
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

	// for delaunay triangulation and crust
	tTetra tetra = NULL;
	tsFace temp_face;
	double volume;

	tVertex voronoi_vertex;
	int vvlnum = 0;
	tList vvl;
	int i = 0;

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
	//qh_setvoronoi_all();
	

	//loop through all faces
	FORALLfacets
	{
		tetra = MakeNullTetra();

		//get vertices of facet
		//loop through each vertex
		vid = 0;
		FOREACHvertex_(facet->vertices)
		{
			//get the id of the vertex
			tetra->vertex[vid++] = all_v[qh_pointid(vertex->point)];
		}

		// for checking the volume of a face.
		temp_face.vertex[0] = tetra->vertex[0];
		temp_face.vertex[1] = tetra->vertex[1];
		temp_face.vertex[2] = tetra->vertex[2];

		volume = fabs(Volumei(&temp_face, tetra->vertex[3]));
		/*temp = qh_facetarea(facet);*/
		if (!volume) {
			DELETE(tetras, tetra);
			continue;
		}

		
		//// for facing-down triangles and abs(volume) > 0 
		//if (facet->normal[3] < 0.0) {
		//	printf("%f	%f	   %f	  %f	 %f\n", facet->normal[0], facet->normal[1], facet->normal[2], facet->normal[3], facet->normal[4]);

		//	MakeFace(tetra->vertex[0], tetra->vertex[1], tetra->vertex[2], NULL);
		//	MakeFace(tetra->vertex[3], tetra->vertex[1], tetra->vertex[0], NULL);
		//	MakeFace(tetra->vertex[2], tetra->vertex[3], tetra->vertex[0], NULL);
		//	MakeFace(tetra->vertex[1], tetra->vertex[2], tetra->vertex[3], NULL);
		
			voronoi_vertex = MakeNullVoronoiVertex();

			voronoi_vertex->v[X] = qh_facetcenter(facet->vertices)[X];
			voronoi_vertex->v[Y] = qh_facetcenter(facet->vertices)[Y];
			voronoi_vertex->v[Z] = qh_facetcenter(facet->vertices)[Z];
			voronoi_vertex->vnum = vvlnum++;

			vvl = MakeNullVoronoiVertexList();
			vvl->p = (void*)voronoi_vertex;

			for (i = 0; i < 4; i++) {
				tVertex tetra_vertex = tetra->vertex[i];
				ADD(tetra_vertex->vvlist, vvl);
			}
			//findPole(vsize, facet);

			
		
		FOREACHneighbor_(facet)
		{
			if (facet < neighbor) {
				double facet_normal[3];
				vid = 0;
				FOREACHvertex_(neighbor->vertices)
				{
					vid++;
					// get normal
					facet_normal[X] += facet->normal[X];
					facet_normal[Y] += facet->normal[Y];
					facet_normal[Z] += facet->normal[Z];

				}
				facet_normal[X] /= vid;
				facet_normal[Y] /= vid;
				facet_normal[Z] /= vid;

			}
		}//FOREACHneighbor_

		
		


		//compute Voronoi centers for all facets
		//includes upperDelaunay facets if qh.UPPERdelaunay ('Qu')
		//qh_setvoronoi_all();	
		// qh_ASvoronoi;

	}
	
		
		//not used
	free(pt);
	free(all_v);
	pt = NULL;
	all_v = NULL;

	//free mem
	qh_freeqhull(!qh_ALL);
	qh_memfreeshort(&curlong, &totlong);


}



