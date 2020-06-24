/***************************************************************************
 *   Copyright (C) 2007 by Joachim Mueller   *
 *   joachim@stardust   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifndef UTILS_H
#define UTILS_H


#include "defines.h"
#include "CBlk.h"

//! for Visualization
struct SBox_Info
{
    INT32 max_level;
    INT32 RB_[3];
    INT32 BlkNds_[3];

    FILE_REAL dt;

    FILE_REAL R_obstacle;
    FILE_REAL BoxOrigin[3];
    FILE_REAL BoxLength[3];



    INT32 Output2D;
    INT32 Output3D;

};

INT32 SFC_Indices_to_BlkNr(INT32 *indices, INT32 level);
void SFC_BlkNr_to_Indices(INT32 BlkNr, INT32 *indices, INT32 level);

INT32 LINEAR_Indices_to_BlkNr(INT32 ind[3]);
void LINEAR_BlkNr_to_Indices(INT32 blk_nr, INT32 ind[3]);


double vec_len(const double *v);
double vec_len(double *v);
float vec_len(float *v);

float vec_len2(float *v);
double vec_len2(double *v);


// template <typename T>
// T vec_scalar(T *v1, T *v2);

double vec_scalar(double *v1, double *v2);
// double vec_scalar(const double *v1, double *v2);
float  vec_scalar(float *v1, float *v2);

void vec_cross(float *res, float *v1, float *v2);
void vec_cross(double *res, double *v1, double *v2);
void vec_cross(double *res, float *v1, double *v2);
void vec_cross(long double *res, long double *v1, long double *v2);


void set_inflow_BField(D_REAL *B, PARTICLE_REAL *x);
void set_inflow_velocity(PARTICLE_REAL *v,  PARTICLE_REAL *x, INT32 species);
void set_inflow_density(PARTICLE_REAL &rho, PARTICLE_REAL *x, INT32 species);


void get_I_GN_equal_process(bool direc, CBlock* dest_Blk, CBlock* src_Blk, INT32 field_type);
void get_J_GN_equal_process(bool direc, CBlock* dest_Blk, CBlock* src_Blk, INT32 field_type);
void get_K_GN_equal_process(bool direc, CBlock* dest_Blk, CBlock* src_Blk, INT32 field_type);



void recv_MPI_Field(INT32 id_field, INT32 id_src_process, INT32 id_package, D_REAL* Field);
void send_MPI_Field(INT32 id_field, INT32 id_request,INT32 id_dest_process, INT32 tag, CBlock* src_Blk);


void cp_I_GN_send_MPI(bool direc, CBlock* dest_Blk, CBlock* src_Blk, INT32 id_package, INT32 field_type);
void cp_J_GN_send_MPI(bool direc, CBlock* dest_Blk, CBlock* src_Blk, INT32 id_package, INT32 field_type);
void cp_K_GN_send_MPI(bool direc, CBlock* dest_Blk, CBlock* src_Blk, INT32 id_package, INT32 field_type);

void cp_I_GN_recv_MPI(bool direc, CBlock* dest_Blk, CBlock* src_Blk, INT32 id_package, INT32 field_type);
void cp_J_GN_recv_MPI(bool direc, CBlock* dest_Blk, CBlock* src_Blk, INT32 id_package, INT32 field_type);
void cp_K_GN_recv_MPI(bool direc, CBlock* dest_Blk, CBlock* src_Blk, INT32 id_package, INT32 field_type);

void add_I_GN_same_process(CBlock* dest_Blk, CBlock* src_Blk, INT32 field_type);
void add_J_GN_same_process(CBlock* dest_Blk, CBlock* src_Blk, INT32 field_type);
void add_K_GN_same_process(CBlock* dest_Blk, CBlock* src_Blk, INT32 field_type);

void add_I_GN_send_MPI(CBlock* dest_Blk, CBlock* src_Blk, INT32 id_package, INT32 field_type);
void add_J_GN_send_MPI(CBlock* dest_Blk, CBlock* src_Blk, INT32 id_package, INT32 field_type);
void add_K_GN_send_MPI(CBlock* dest_Blk, CBlock* src_Blk, INT32 id_package, INT32 field_type);

void add_I_GN_receive_MPI(CBlock* dest_Blk, CBlock* src_Blk, INT32 id_package, INT32 field_type);
void add_J_GN_receive_MPI(CBlock* dest_Blk, CBlock* src_Blk, INT32 id_package, INT32 field_type);
void add_K_GN_receive_MPI(CBlock* dest_Blk, CBlock* src_Blk, INT32 id_package, INT32 field_type);

void im1_BV_from_parent(CBlock* child_block, INT32 field_type);
void ip1_BV_from_parent(CBlock* child_block, INT32 field_type);
void jm1_BV_from_parent(CBlock* child_block, INT32 field_type);
void jp1_BV_from_parent(CBlock* child_block, INT32 field_type);
void km1_BV_from_parent(CBlock* child_block, INT32 field_type);
void kp1_BV_from_parent(CBlock* child_block, INT32 field_type);


bool calc_cells_to_cut(INT32 X, INT32 *cell_indedx_L);

bool calc_BlockIndex_of_Pos(INT32**   index_blk,
		            D_REAL*  r_vec);


bool calc_BlockNodeXCell_of_Pos(INT32**  index_blk,
				INT32**  index_cell,
				D_REAL** x_cell,
				D_REAL*  r_vec);


void get_maxwellian_v(INT32 species, PARTICLE_REAL *vth);
void get_maxwellian_GSL(INT32 species, PARTICLE_REAL *vth);
void get_maxwellian_GSL(INT32 species, D_REAL ion_beta, PARTICLE_REAL *vth);


D_REAL calc_reaction_rate(PARTICLE_REAL* v, INT32 species, INT32 dest_species, INT32 neutral_species);




#endif
