#include "shape.h"

/* Global variable definitions */

//defined in chull.c
extern tVertex vertices;
extern tEdge edges;
extern tFace faces;
extern tTetra tetras;
extern bool debug;
extern bool check;


/*-------------------------------------------------------------------*/
void	Crust( void );
double getAngle(tVertex v1, tVertex v2);
double innerProduct(tVertex site, tVertex pole, tVertex v);

