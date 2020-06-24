

#include <iostream>
#include <fstream>
#include <math.h>


//! silo 2D functions
void silo_2DWrite_allCS(INT32 Field_ID);
void silo_2DWrite_prepare(void);
void silo_2DWrite_MeshCS(INT32 CS);
void silo_2DWrite_FieldCS(INT32 Field_ID, INT32 CS);
void silo_2DFlush_CleanUp(void);
void silo_2DWrite_Trajectory(INT32 Trajectory_ID, INT32 num_positions, double** positions, INT32 CS);
void silo_2DWrite_ParticleTrack(INT32 GROUP_ID, INT32 num_positions, double** positions, INT32 CS);
void silo_2DWrite_Line(INT32 Trajectory_ID, INT32 num_positions, double** positions, INT32 CS);



//! silo 3D functions
void silo_3DWrite_prepare(void);
void silo_3DWrite_Mesh(void);
void silo_3DWrite_Field(INT32 Field_ID);
void silo_3DFlush_CleanUp(void);
void silo_3DWrite_Trajectory(INT32 Trajectory_ID, INT32 num_positions, double** positions);
void silo_3DWrite_ParticleTrack(INT32 GROUP_ID, INT32 num_positions, double** positions);
void silo_3DWrite_Line(INT32 Trajectory_ID, INT32 num_positions, double** positions);



void silo_TestMesh(void);

void convert_particle_tracks_to_silo_mesh(void);