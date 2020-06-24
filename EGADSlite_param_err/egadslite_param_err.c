/* ------------------------------------------------------------------------  */
/*   Simple program to illustrate EG_invEvaluate() sometimes gives parameter */
/*   values outsite the range of a Topology object. In this case, and EDGE   */
/*                                                                           */
/*   Notes:    * There is only 1 body in the EGADSlite model                 */
/*             * The coordinates in xyz[] are the coordinates of vertex 14   */
/*                 in the model. They are not exact but duplicate the issue  */
/* ------------------------------------------------------------------------  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "egads.h"

int main(int argc, char *argv[])
{
  ego      context, model, geom, *bodies, body, edge;
  int      ierr, bodyID = 1, edgeID = 18, peri;
  int      oclass, mtype, nbodies, *bsenses;
  double   range[4], params[4], results[3];
  double   xyz[3] = {10.647474, -20.650200, 6.789684};
  
  /* Initialize & read in Model Data */
  printf(" EG_open       = %d \n", EG_open(&context));
  printf(" EG_loadModel  = %d %s \n", EG_loadModel(context, 0, argv[1], &model), argv[1]);
  
  ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &bsenses);
  body = bodies[0];
  
  ierr = EG_objectBodyTopo(body, EDGE, edgeID, &edge);
  ierr = EG_invEvaluate(edge, xyz, params, results);
  ierr = EG_getRange(edge, range, &peri);
  
  // Print out results to screen
  printf(" \n");
  printf(" Edge %d \n", edgeID);
  printf("  range(tmin, tmax) = (%lf, %lf) \n", range[0], range[1]);
  printf("                xyz = (%lf, %lf, %lf) \n", xyz[0], xyz[1], xyz[2]);
  printf("                  t = %lf \n", params[0]);
  printf(" \n");  
}