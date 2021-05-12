				EGADS Utility README


This subdirectory contains parts of EGADS that may be useful on their own, and 
code that augments EGADS functionality but is not considered integral.

EMP is the Explicit Multithread Package that is a small, simple machine
independent API that supports threaded coding. This is used in EGADS primarily 
during Face tessellation, where each thread builds a Face’s triangulation 
independently and concurrently within a Body. The entry points are part of the 
EGADS library but can clearly be used by itself. ThreadTest.c is code that 
displays the use of EMP and is also the unit tester.

retessFaces is code, that when built with EGADS, can retriangulate specific 
Faces (in a Body) using different tessellation parameters.

limitTessBody is example code that produces an EGADS’ Body tessellation but 
limits the Edge/Face vertex counts to those specified. The makefile 
‘limits.make’ (‘limits.mak’ for Windows nmake) builds a function that can be 
included along with the EGADS library and also a test executable. The tester 
loads an EGADS model from a file, determines the number of vertices desired per
entity, then invokes limitTessBody and finally reports the results.

egadsHOtess.c and the test program vHOtess takes an egads Tessellation object
and makes a new tessellation object by inserting the "high order" vertices --
the test program displays the results.

extractTess shows how to generate a multi-body tessellation that has Face and
Edge coherence.

uvmap takes a triangulation and generates a single uv map. This is from David
Marcum (MSU) and has its own license.
