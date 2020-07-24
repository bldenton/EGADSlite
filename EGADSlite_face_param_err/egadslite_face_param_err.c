/* ------------------------------------------------------------------------  */
/*   Simple program to illustrate EG_invEvaluate() sometimes gives parameter */
/*   values outsite the range of a Topology object. In this case, a FACE     */
/*                                                                           */
/*   Notes:    * There is only 1 body in the EGADSlite model                 */
/*             * The coordinates in xyz[] are the coordinates of location    */
/*                 in the model. They are not exact but duplicates the issue */
/* ------------------------------------------------------------------------  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "egads.h"

int main(int argc, char *argv[])
{
  ego      context, model, geom, *bodies, body, face;
  int      ierr, bodyID = 1, faceID = 4, peri;
  int      oclass, mtype, nbodies, *bsenses;
  double   range[4], params[4], results[3], eval[18];
  double   xyz[3] = {-0.020576, 0.215900, 0.020576};
  
  /* Initialize & read in Model Data */
  printf(" EG_open       = %d \n", EG_open(&context));
  printf(" EG_loadModel  = %d %s \n", EG_loadModel(context, 0, argv[1], &model), argv[1]);
  
  ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &bsenses);
  body = bodies[0];
  
  ierr = EG_objectBodyTopo(body, FACE, faceID, &face);
  ierr = EG_invEvaluate(face, xyz, params, results);
  ierr = EG_getRange(face, range, &peri);
  ierr = EG_evaluate(face, params, &eval);
  
  // Print out results to screen
  printf(" \n");
  printf(" FACE %d \n", faceID);
  printf("  range(umin, umax, vmin, vmax) = (%lf, %lf, %lf, %lf) \n", range[0], range[1], range[2], range[3]);
  printf("                            xyz = (%lf, %lf, %lf) \n", xyz[0], xyz[1], xyz[2]);
  printf("                         (u, v) = (%lf, %lf) \n", params[0], params[1]);
  printf("                        results = (%lf, %lf, %lf) \n", results[0], results[1], results[2]);
  printf("                           eval = (%lf, %lf, %lf) \n", eval[0], eval[1], eval[2]);
  printf(" \n");  
}