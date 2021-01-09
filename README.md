# EGADS/EGADSlite Libraries for Petsc

Project consists of a reorganization of the EGADS/EGADSlite code developed by Dr. Haimes at MIT [ESP: Engineering Sketch Pad](https://acdl.mit.edu/ESP).
The current EGADS/EGADSlite version v1.18

The reorganization provided in this repository is intended to allow for:

  1. Concurrent and/or Independent use of the EGADS and EGADSlite libraries
     - Original coding resulted in conflicts due to functions having the same name & signature
	 - Resolution was achieved by:
	   - Copying & Renaming shared files (EGADSlite files are now identified with a _lite suffix)
	   - Renaming EGADSlite functions (Now start with EGlite_ prefix)
  2. Allow integration with [Petsc](https://www.mcs.anl.gov/petsc)
     - Integration allows Petsc to:
       - Read STEP, IGES, EGADS and EGADSlite files
	   - Generate 2D & 3D Meshes of the geometry
	     - Leans on Tetgen or cTetgen for 3D mesh generation
	   - Apply boundary conditions via geometric entity reference
	   - Snap-to-Geometry refinement of meshes

# Installation

Compiling of the EGADS/EGADSlite libraries are as follows

## Requirements

  - EGADSlite
    - Required Software
	  - C Compiler
	- Optional Software
	  - Python
  - EGADS
    - Require Software
	  - C Compiler
      - [OpenCascade v7.5.0](https://dev.opencascade.org/)
	    - Acquistions options outlined in `./OpenCascade-v7.5.0/OpenCascadev7.5.0.txt` of this repository
	    - OpenCascade compilation requires the following:
	      - GNU gcc 4.3+ or LLVM CLang 3.6+
		  - OpenGL 3.3+ or OpenGL ES 2.0+
		  - FreeType 2.4.11-2.7.1
		  - Tcl/Tk 8.6.3+ or ActiveTcl 8.6
		  - Doxygen 1.8.5+
		  - CMake
		  - Note: For optional additional 3rd Party libraries and tools see [here](https://old.opencascade.com/doc/occt-7.5.0/overview/html/index.html)
	- Optional Software
	  - Python
  - Petsc Integration
    - Petsc will automatically install and integrate the EGADSlite library when the `--download-egads=1` flag is used
	  - Note: It is recommended that the `--with-tetgen=true` and `--download-ctetgen=1` flags are also used for full meshing functionality

## OpenCascade

Compilation & Installation of OpenCascade is required **only** if compiling EGADS library. The EGADSlite library does not use OpenCascade.

  - Unzip tarball
    - See `./OpenCascade-v7.5.0/OpenCascadev7.5.0.txt` for options on getting the tarball or cloning repository
    - `tar -zxf opencascade-7.5.0.tgz` or clone `git clone https://github.com/bldenton/oce`
  - Create build director (Out of source builds are highly recommended)
    - `mkdir build`
  - Go to build directory & run cmake
    - `cd build`
	- `cmake ../opencascade-7.5.0 -DINSTALL_DIR=path/to/where/you/want/to/install`
  - Compile code with Make
    - `make`
	- Note: Supports parallel build. Can add -j#
  - Install Opencascade
    - `make install/strip`

OpenCascade is now installed. Take note of OpenCascade's installation directory. It will be needed for EGADS compilation.

## EGADS/EGADSlite Environment Setup

Environment Variables used during compilation and use of EGADS/EGADSlite libraries are generated as follows.

  - Generate EGADS/EGADSlite environment
    - EGADSlite library Only
	  - `./makeEnv`
	- EGADS library Only
	  - `./makeEnv /path/to/opencascade/installation/top/folder/`
	- EGADS & EGADSlite
	  - `./makeEnv /path/to/opencascade/installation/top/folder/`

This script creates the `./config`, `./include` and `./lib` folders that will store the relevant EGADS/EGADSlite files (./include) and libraries (./lib).
It will also generate files iESPenv.csh (tcsh) and iESPenv.sh (bash) scripts in the `./config` folder. The appropriate script for your system must
be sourced prior to compilation or use of the EGADS and/or EGADSlite libraries.

  `source ./config/iESPenv.sh`  or  `source ./config/iESPenv.csh`

## EGADSlite Only

Compilation of the EGADSlite library is achieved as follows:

  - Source appropriate ESP environment script
    - `source ./config/iESPenv.sh` or `source ./config/iESPenv.csh`
  - Run make
    - `make egadslite`

The required header files will be located in `./include` and the library will be located in `./lib`

## EGADS Only

Compilation of the EGADS library is achieved as follows:
**NOTE**: OpenCascade v7.5.0 installation is required. See above.

  - Source appropriate ESP environment script
    - `source ./config/iESPenv.sh` or `source ./config/iESPenv.csh`
  - Run make
    - `make egads`

The required header files will be located in `./include` and the library will be located in `./lib`

## EGADS & EGADSlite

Compilation of the EGADS & EGADSlite libraries are achieved as follows:
**NOTE**: OpenCascade v7.5.0 installation is required. See above.

  - Source appropriate ESP environment script
    - `source ./config/iESPenv.sh` or `source ./config/iESPenv.csh`
  - Run make
    - `make`

The required header files will be located in `./include` and the libraries will be located in `./lib`

## Petsc Integration

Petsc will automatically install and integrate the EGADSlite library when the `--download-egads=1` flag is used.
  - Note: It is recommended that the `--with-tetgen=true` and `--download-ctetgen=1` flags are also used for full meshing functionality.

Integration of the EGADS library is in work and will be online soon. For instructions on how to use Petsc with EGADSlite please review the Petsc documentation
at [Petsc Documentation](https://www.mcs.anl.gov/petsc/documentation/index.html)


# Credits

- Robert Haimes, MIT, [ESP: Engineering Sketch Pad](https://acdl.mit.edu/ESP) source code (EGADS/EGADSlite)
- [Dr. Matthew Knepley](https://cse.buffalo.edu/~knepley/), University @ Buffalo, Ph.D. Advisor & [Petsc](https://www.mcs.anl.gov/petsc) Integration assistance 
- [OpenCascade Project](https://dev.opencascade.org/)





		
=======
# EGADSlite
Project is intended to integrate the capabilities of EGADSlite into Petsc. 

When used with Petsc (option -with-egads = 1), Petsc will read a .egadslite CAD file and generate a 2D DM Plex with the CAD topology embedded.
The embedded information allows the Plex (mesh) to "snap to the geometry" the Plex is representing upon refinement. This technology allows the
Plex (mesh) to better approximate the domain and/or boundary conditions as it is refined.

The initially created 2D DM Plex can be used to generate a 3D Plex (mesh) via Petsc's implementation of cTetgen or Tetgen. In this case, the CAD
topology is transferred from the 2D Plex (mesh) to the 3D Plex (mesh) and gains the "snap to geometry" capability. In addition, boundary conditions
can be assigned via geometric association.

Additional Capabilies include:
     1) Inflate to Geometery :: When 3D Plexes (meshes) are generated via cTetgen or Tetgen with Steiner points, it must be "inflated" so the
                                Steiner points lie on the associated CAD surfaces.
     2) Print Model Topology :: Print Model Topology to screen
     
