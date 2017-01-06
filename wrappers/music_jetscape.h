// Copyright @ Bjoern Schenke, Sangyong Jeon, Charles Gale, and Chun Shen
#ifndef TEST_MUSIC_JETSCAPE_H_
#define TEST_MUSIC_JETSCAPE_H_

#include "music.h"
#include "fluid_dynamics.h"

class MPI_MUSIC: public FluidDynamics {
    // this is wrapper class for MUSIC so that it can be used as a external
    // library for the JETSCAPE integrated framework
 private:
    int mode;            // records running mode
    MUSIC *music_hydro_ptr;

 public:
     MPI_MUSIC();
     ~MPI_MUSIC();

     void initialize_hydro(Parameter parameter_list);

     void evolve_hydro();
     void get_hydro_info(real t, real x, real y, real z,
                         FluidCellInfo* fluid_cell_info_ptr) {};
     void get_hypersurface(real T_cut, SurfaceCellInfo* surface_list_ptr) {};
};

#endif  // TEST_MUSIC_JETSCAPE_H_
