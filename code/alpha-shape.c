#include "alpha-shape.h"
#include "shape.h"

/*-------------------------------------------------------------------*/
//headers from qhull
#include "qhull.h"
#include "poly.h"
#include "qset.h"
#include "geom.h"


/*-------------------------------------------------------------------*/
//defined in shape.h/.c
extern tVertex vertices;
extern tEdge edges;
extern tFace faces;
extern tTetra tetras;

double pointDist(tVertex v1, tVertex v2) {
	double pt1[3] = { v1->v[X], v1->v[Y], v1->v[Z]};
	double pt2[3] = { v2->v[X], v2->v[Y], v2->v[Z]};

	return qh_pointdist(pt1, pt2, 3);
}

/* 
	calculate the area of a triangle whose sides are products of the opposite edge of the tetra
	There is a tetra.
	Let the length of edges be p, q, r, a, b, c.
	Triangle efg, e = ap, f = bg, g = cr
*/
double areaTriangle(tTetra tetra) {
	double area = 0.0;
	double p, q, r, a, b, c;
	double e, f, g;
	double s;

	p = pointDist(tetra->vertex[0], tetra->vertex[1]);
	q = pointDist(tetra->vertex[0], tetra->vertex[2]);
	r = pointDist(tetra->vertex[0], tetra->vertex[3]);
	a = pointDist(tetra->vertex[2], tetra->vertex[3]);
	b = pointDist(tetra->vertex[3], tetra->vertex[1]);
	c = pointDist(tetra->vertex[1], tetra->vertex[2]);
	
	e = a * p;
	f = b * q;
	g = c * r;
	//printf("%f, %f, %f, %f, %f, %f\n", p, q, r, a, b, c);
	s = (e + f + g) * 0.5;
	area = sqrt(s * (s - e) * (s - f) * (s - g));

	return area;
}

/*-------------------------------------------------------------------*/
void AlphaShape( unsigned int alpha )
{
	//implement 3d Alpha Shape alogorithm here
	//
	// put results (vertices, edges, faces, and tetrahedra)
	// in the global structures above
	// see main.c for detail on how qhull does it
	//
	tVertex  ptr_v;
	tVertex * all_v = NULL;
	int vsize = 0;
	int id = 0;

	//global varibles for qhull
	static char * options = (char *)"delaunay QJ Pp";
	int curlong, totlong;
	coordT * pt = NULL;
	facetT *facet = NULL;
	vertexT *vertex = NULL;
	vertexT **vertexp = NULL;
	facetT *neighbor, **neighborp;
	int vid = 0;

	// for delaunay triangulation
	tTetra tetra;
	tsFace face;
	double volume;
	double area;
	double radius;
	double squared_radius;

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

	// qh DELAUNAY= True;     /* 'd'   */
	//qh SCALElast= True;    /* 'Qbb' */
	//qh KEEPcoplanar= True; /* 'Qc', to keep coplanars in 'p' */

	qh_initflags(options);
	qh_init_B(pt, vsize, 4, false);
	qh_qhull();
	qh_check_output();

	//loop through all faces
	FORALLfacets
	{
		tetra = MakeNullTetra();	// make a tetra

		//get vertices of facet
		//loop through each vertex
		vid = 0;
		FOREACHvertex_(facet->vertices)
		{
			//get the id of the vertex
			tetra->vertex[vid++] = all_v[qh_pointid(vertex->point)];
		}

		face.vertex[0] = tetra->vertex[0];
		face.vertex[1] = tetra->vertex[1];
		face.vertex[2] = tetra->vertex[2];
		
		volume = fabs(Volumei(&face, tetra->vertex[3]));

		if (facet->normal[3] < 0.0 && volume) {
			area = areaTriangle(tetra);	// same withqh_facetarea(tetra->face);
			radius = area / (6 * volume);	// same with qh_pointdist(tetra->vertex[0], qh_facetcenter(facet->vertices), 3);
			//radius = qh_pointdist(tetra->vertex[0], qh_facetcenter(facet->vertices), 3);
			
			squared_radius = SQR(radius);
			if (squared_radius > alpha) {// squared_radius <= alpha
				DELETE(tetras, tetra);
				continue;
			}
			tetra->face[0] = MakeFace(tetra->vertex[0], tetra->vertex[1], tetra->vertex[2], NULL);
			tetra->face[1] = MakeFace(tetra->vertex[3], tetra->vertex[1], tetra->vertex[0], NULL);
			tetra->face[2] = MakeFace(tetra->vertex[2], tetra->vertex[3], tetra->vertex[0], NULL);
			tetra->face[3] = MakeFace(tetra->vertex[1], tetra->vertex[2], tetra->vertex[3], NULL);
		}
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

