#include "crust.h"
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

/*-------------------------------------------------------------------*/
void	Delaunay( void )
{
	printf("delaunay\n");
	// Use qhull to compute Delaunay triangulation of the points.
	// put results (vertices, edges, faces, and tetrahedra)
	// in the global structures above
	// see main.c for detail on how qhull does it


	//
	// create points in 4D (x,y,z,x^2+y^2+z^2)
	//

	//
	// compute convex hull in 4D by calling qhull
	// use flags: static char * options=(char *)"delaunay QJ Pp";
	//

	//
	//loop through all faces and call MakeNullTetra() to make a tetra 
	//remember that this is in 4D, so each face is a tetrahedron
	//
	//to fill the teta: get vertices of facet and loop through each vertex
	//use FOREACHvertex_()
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
	
	tTetra tetra;
	tsFace face;

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

		if (facet->normal[3] < 0.0 && Volumei(&face, tetra->vertex[3])) {
			MakeFace(tetra->vertex[0], tetra->vertex[1], tetra->vertex[2], NULL);
			MakeFace(tetra->vertex[3], tetra->vertex[1], tetra->vertex[0], NULL);
			MakeFace(tetra->vertex[2], tetra->vertex[3], tetra->vertex[0], NULL);
			MakeFace(tetra->vertex[1], tetra->vertex[2], tetra->vertex[3], NULL);
		}

		//to fill the teta: get vertices of facet and loop through each vertex
		//use FOREACHvertex_()
		FOREACHneighbor_(facet)
		{
			if (facet < neighbor) {
				tVertex vertices[4];
				vid = 0;
				FOREACHvertex_(neighbor->vertices)
				{
					//get vertex
					vertices[vid++] = all_v[qh_pointid(vertex->point)];
				}

			}
		}//FOREACHneighbor_
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



