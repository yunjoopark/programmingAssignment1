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

/*******************************************/
/* tetracircumcenter() Find the circumcentre of the tetra */
/*******************************************/

double areaTriangle(tTetra tetra) {
	double area = 0.0;
	double p, q, r, a, b, c;
	double e, f, g;
	double s;

	p = qh_pointdist(tetra->vertex[0], tetra->vertex[1], 4);
	q = qh_pointdist(tetra->vertex[0], tetra->vertex[2], 4);
	r = qh_pointdist(tetra->vertex[0], tetra->vertex[3], 4);
	a = qh_pointdist(tetra->vertex[2], tetra->vertex[3], 4);
	b = qh_pointdist(tetra->vertex[3], tetra->vertex[1], 4);
	c = qh_pointdist(tetra->vertex[1], tetra->vertex[2], 4);
	
	e = a * p;
	f = b * q;
	g = c * r;
	
	s = (e + f + g) *0.5;
	area = sqrt(s * (s - e) * (s - f) * (s - g));

	return area;
}

void drawCircumsphere(double radius) {
	tFace f;
	if (faces == NULL) return;
	f = faces;
	
	const double PI = 3.14159265359;
	double cir_delta = PI * 2;
	glVertex3d(radius, radius, 0);


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
	//double pt4 = 0;
	tTetra tetra;
	tsFace face;
	double volume;
	double area;
	double radius;

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
		
		volume = Volumei(&face, tetra->vertex[3]);
		area = areaTriangle(tetra);
		radius = area / (6 * volume);

		if (facet->normal[3] < 0.0 && volume) {
			tetra->face[0] = MakeFace(tetra->vertex[0], tetra->vertex[1], tetra->vertex[2], NULL);
			tetra->face[1] = MakeFace(tetra->vertex[3], tetra->vertex[1], tetra->vertex[0], NULL);
			tetra->face[2] = MakeFace(tetra->vertex[2], tetra->vertex[3], tetra->vertex[0], NULL);
			tetra->face[3] = MakeFace(tetra->vertex[1], tetra->vertex[2], tetra->vertex[3], NULL);
		}
		if (radius < alpha) {
			DELETE(tetras, tetra);
		}
		//to fill the teta: get vertices of facet and loop through each vertex
		//use FOREACHvertex_()
		FOREACHneighbor_(facet)
		{
			if (facet < neighbor) {
				tVertex vertices[4];
				vid = 0;

//				int a = qh_getdistance(facet, neighbor, alpha, alpha);
//				printf("%d\n", a);

				FOREACHvertex_(neighbor->vertices)
				{
					//get vertex
					vertices[vid++] = all_v[qh_pointid(vertex->point)];
				}
			}
		}//FOREACHneighbor_
	}
	

	

	//if (diameter > alpha) {

	//}


	//not used
	free(pt);
	free(all_v);
	pt = NULL;
	all_v = NULL;

	//free mem
	qh_freeqhull(!qh_ALL);
	qh_memfreeshort(&curlong, &totlong);
}

