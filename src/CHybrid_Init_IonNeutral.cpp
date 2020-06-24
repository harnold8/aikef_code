#include <math.h>
#include "utils.h"
#include "CHybrid.h"
#include "parameters.h"
#include "absolute_Globals.h"
#include <ctime>

#include "unistd.h"

#include <iostream>
#include <fstream>
#include <sstream>

//! time to wait in micro-seconds
#define TIME_TO_WAIT 200000


//! NOTE Apart from velocity dependence of chemical reactions, no modifications to this file are neccessary


#if defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)
//!--------------------------------------------------------
//!- init_Neutral_Profile:
//!--------------------------------------------------------
void CHybrid::init_Neutral_Profile(void)
{


	log_file <<  endl;
	log_file << " INITIALIZING NEUTRAL PROFILE ... " << endl;


	for(INT32 neutralSpec = 0;neutralSpec<num_Neutral_Species;neutralSpec++)	
	{	
		//! Enceladus
		if(neutral_profile_from_file)	
		read_extern_IonProfile_uniform_grid(neutralSpec,id_numberdensity_neutralSpecies1,id_velocity_neutralSpecies1);
					
		if(analytical_neutral_profile)
		set_analytical_neutral_profile(neutralSpec);
	}		

	log_file << " done. " << endl << endl;


}



//!------------------------------------------------------------------------
//!- set analytical neutral profile to neutral species field:
//!------------------------------------------------------------------------
void CHybrid::set_analytical_neutral_profile(INT32 neutral_species)
{


	log_file << "  Setting analytical neutral profile to field " << Field_Name[id_numberdensity_neutralSpecies1+neutral_species] << "...       ";
	
	
	//! in case density should not be added to existing density,
	//! field hast to be set to zero
	if(!neutral_profile_from_file && analytical_neutral_profile)
	set_zero_field(id_numberdensity_neutralSpecies1 + neutral_species);

	
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
	
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
			if(mpi_myRank==temp_Block->responsible_mpi_process)
			temp_Block->set_analytical_neutral_profile(neutral_species);
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	
	}
	
	log_file << "done." << endl;

}


#endif



#if defined(use_ion_production_as_field)
//!--------------------------------------------------------
//!- init_Neutral_Profile:
//!--------------------------------------------------------
void CHybrid::init_IonProduction_Profile(void)
{


	log_file <<  endl;
	log_file << " INITIALIZING ION PRODUCTION PROFILE ... " << endl;


	for(INT32 neutralSpec = 0;neutralSpec<num_ion_prod_fields;neutralSpec++)	
	{	
		//! Enceladus
		if(ion_prod_profile_from_file)	
		read_extern_IonProfile_uniform_grid(neutralSpec,id_density_ionProdSpecies1+neutralSpec,id_velocity_ionProdSpecies1+neutralSpec);

		if(calc_ionProd_from_neutral_field)
		calculate_Chapman();			

		if(set_analytical_ionProd_profile)
		set_ion_production_profile(neutralSpec);
	}		

	log_file << " done. " << endl << endl;


}



//!------------------------------------------------------------------------
//!- set analytical neutral profile to neutral species field:
//!------------------------------------------------------------------------
void CHybrid::set_ion_production_profile(INT32 neutral_species)
{


	log_file << "  Setting ion production profile to field " << Field_Name[id_density_ionProdSpecies1+neutral_species] << "   ...    ";
	
	
	//! in case density should not be added to existing density,
	//! field hast to be set to zero
	if(!ion_prod_profile_from_file || !calc_ionProd_from_neutral_field)
	set_zero_field(id_density_ionProdSpecies1 + neutral_species);

#if defined(use_neutral_species_as_field)
	if(!calc_ionProd_from_neutral_field)
	{	
		copy_Field(id_density_ionProdSpecies1 + neutral_species, id_numberdensity_neutralSpecies1+neutral_species);
		copy_Field(id_velocity_ionProdSpecies1 + neutral_species, id_velocity_neutralSpecies1+neutral_species);
	}	
#endif
	
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
	
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
			if(mpi_myRank==temp_Block->responsible_mpi_process)
			temp_Block->set_ion_production_profile(neutral_species);
		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	
	}
	
	log_file << "done." << endl;

}


#endif

//!--------------------------------------------------------
//!- read uniform_grid output as input for extern neutral/ion profile
//!  uniform_grid output must only consist of density and velocity
//! (in this order)
//!--------------------------------------------------------
void CHybrid::read_extern_IonProfile_uniform_grid(INT32 neutralSpec, INT32 id_densityfield, INT32 id_velocityfield)
{

	
	log_file << endl << endl;
	log_file << " Reading extern Ion Profile from file:" << endl <<"  "
	     << extern_IonProdField_name[neutralSpec] << endl
	     << " and storing it to field '"
	     << Field_Name[id_densityfield] << " ' ...  " << endl;

	INT32 tag = 0;
	INT32 is_finished = 0;
	INT32 length_in_byte = 1;
	if(serialize_reading_of_extern_field)
	{
		//! In order to avoid that all processes write at the same time,
		//! let p0 write first, than send message to p1 to trigger writing,
		//! then p2 .... up to pN

		//! - every process must receive except for first process
		//! - use Blocking receive
		if(serialize_writing_of_StateFile && mpi_myRank)
		mpi_myComm.Recv(&is_finished, length_in_byte, MPI_INT32, mpi_myRank-1, tag );
	}
	     
	     
	INT32 num_extern_box_nodes, num_extern_nodes[3];
	short num_extern_Nodes[3];
	FILE_REAL extern_Origin[3], extern_Length[3];
	D_REAL origin[3], length[3];
	char temp[42];
	int TL;
	FILE_REAL Radius, SI_Quantities[4];
	D_REAL *extern_rho, *extern_UiX, *extern_UiY, *extern_UiZ;


	char full_filename[200];
	sprintf(full_filename,"%s",extern_IonProdField_name[neutralSpec]);
	
	log_file << " full_filename " << full_filename << endl;

	ifstream extern_rho_file;
	extern_rho_file.open(full_filename, ifstream::binary);

	if(!extern_rho_file)
	{
		log_file << "no extern_rho file" << endl;
		finalize_MPI();
		exit(1);
	}


	//! start reading
	extern_rho_file.read(reinterpret_cast<char*> (temp), 42*sizeof(char));
	//! read number of mesh nodes on which extern_rho is defined
	extern_rho_file.read(reinterpret_cast<char*> (num_extern_Nodes), 3*sizeof(short));
	extern_rho_file.read(reinterpret_cast<char*> (&TL), sizeof(int));
	//! read Box Origin
	extern_rho_file.read(reinterpret_cast<char*> (extern_Origin), 3.*sizeof(FILE_REAL));
	//! read Box Length
	extern_rho_file.read(reinterpret_cast<char*> (extern_Length), 3.*sizeof(FILE_REAL));
	extern_rho_file.read(reinterpret_cast<char*> (&Radius), sizeof(FILE_REAL));
	extern_rho_file.read(reinterpret_cast<char*> (SI_Quantities), 4.*sizeof(FILE_REAL));

	for(int i=0;i<3;i++)
	{
		num_extern_nodes[i]= num_extern_Nodes[i];
		origin[i] = extern_Origin[i]*1560800./SI_x0/*R_Obstacle*/;
		length[i] = extern_Length[i]*1560800./SI_x0/*R_Obstacle*/;
	}
	

	log_file << "   Run-Name of extern profile : " << temp << endl;
	
	log_file << " 	Obstacle Radius of extern profile : " << Radius << endl;

	log_file << " 	Time level of extern profile : " << TL << endl;
	
	log_file << " 	x0 : " << SI_Quantities[0] << endl;
	log_file << " 	v0 : " << SI_Quantities[1] << endl;
	log_file << " 	B0 : " << SI_Quantities[2] << endl;
	log_file << " 	n0 : " << SI_Quantities[3] << endl;

	log_file << "  ->extern Nodes X: " << num_extern_nodes[0]<< endl;
	log_file << "  ->extern Nodes Y: " << num_extern_nodes[1]<< endl;
	log_file << "  ->extern Nodes Z: " << num_extern_nodes[2]<< endl;
	log_file << endl;

	log_file << "  ->extern LX=[ "<< -origin[0]<<" : "<<length[0]-origin[0]<<" ]"<< endl;
	log_file << "  ->extern LY=[ "<< -origin[1]<<" : "<<length[1]-origin[1]<<" ]"<< endl;
	log_file << "  ->extern LZ=[ "<< -origin[2]<<" : "<<length[2]-origin[2]<<" ]"<< endl;
	log_file << endl;

	log_file << "  ->intern LX=[ "<< -Box_Origin[0]<<" : "<<LX-Box_Origin[0]<<" ]"<< endl;
	log_file << "  ->intern LY=[ "<< -Box_Origin[1]<<" : "<<LY-Box_Origin[1]<<" ]"<< endl;
	log_file << "  ->intern LZ=[ "<< -Box_Origin[2]<<" : "<<LZ-Box_Origin[2]<<" ]"<< endl;
		log_file << endl;

		num_extern_box_nodes =   num_extern_nodes[0]
					*num_extern_nodes[1]
					*num_extern_nodes[2];
		

	//! in uniform_grid the positions of the nodes is also stored and have to be read
// 	FILE_REAL* mesh_nds_pos = new FILE_REAL[3*num_extern_box_nodes];
// 	extern_rho_file.read( reinterpret_cast<char*> (mesh_nds_pos),3*num_extern_box_nodes*sizeof(FILE_REAL));
// 	delete[] mesh_nds_pos;
	//! but this is not neccessary if Box Origin is given relative to Box length
	//! therefore, skip bytes
// 	extern_rho_file.seekg(3*num_extern_box_nodes*sizeof(FILE_REAL),ios::cur);
	
	INT32 uniform_grid_Field_Name_Size = 50;
	char temp2[50];
	int ftype;
	short num_comps;

	extern_rho = new D_REAL[4*num_extern_box_nodes];
	extern_UiX = extern_rho +num_extern_box_nodes;
	extern_UiY = extern_UiX +num_extern_box_nodes;
	extern_UiZ = extern_UiY +num_extern_box_nodes;
	//! MF: I only read ionprod and dont need velocity, so set to zero all and then get only the rate
	memset(extern_rho,0.,4*num_extern_box_nodes*sizeof(D_REAL));
	
	//! read rho
	extern_rho_file.read( reinterpret_cast<char*> (&num_comps),sizeof(short));
	extern_rho_file.read( reinterpret_cast<char*> (&ftype),sizeof(int));
	extern_rho_file.read( reinterpret_cast<char*> (temp2), uniform_grid_Field_Name_Size*sizeof(char));
	extern_rho_file.read( reinterpret_cast<char*> (extern_rho), num_comps*num_extern_box_nodes*sizeof(double));


	//! read U
// 	extern_rho_file.read( reinterpret_cast<char*> (&num_comps),sizeof(short));
// 	extern_rho_file.read( reinterpret_cast<char*> (&ftype),sizeof(int));
//         extern_rho_file.read( reinterpret_cast<char*> (temp2), uniform_grid_Field_Name_Size*sizeof(char));
//         extern_rho_file.read( reinterpret_cast<char*> (extern_UiX), num_comps*num_extern_box_nodes*sizeof(FILE_REAL));
  
	extern_rho_file.close();

	//! use indices a,b,c for extern mesh
	for(INT32 a_b_c=0; a_b_c<num_extern_box_nodes; a_b_c++)
	{
		//! normalization
		extern_rho[a_b_c] *= norm_externIonProdRate[neutralSpec];
		if(extern_rho[a_b_c]<0)log_file<<"error extern rho value negative"<<endl;
// 		extern_UiX[a_b_c] *= norm_externNeutralVelocity[neutralSpec];
// 		extern_UiY[a_b_c] *= norm_externNeutralVelocity[neutralSpec];
// 		extern_UiZ[a_b_c] *= norm_externNeutralVelocity[neutralSpec];
	}


	INT32 num_values_not_in_extern_box = 0;
	//!-----------------------------------------//
	//! Loop over Blks to set rho_extern	   //
	//!-----------------------------------------//
	
	//! NOTE: in order to avoid boundary problems, it may be better not set
	//! 	  rho_extern on Level0
	//!	  But be careful, if there is only one level...
	for(INT32 level=0; level<=MAX_LEVEL; level++)
	{
		CBlock* temp_Block = BlockList_of_Lev[level];
		while(temp_Block)
		{
			if(mpi_myRank == temp_Block->responsible_mpi_process)
			temp_Block->set_RhoUi_extern(id_densityfield,id_velocityfield,
						       num_extern_nodes,
						       origin,
						       length,
						       extern_rho,
						       num_values_not_in_extern_box);

		
			temp_Block = temp_Block->next_Blk_of_BlockList;
		}
	}
		
	if(num_values_not_in_extern_box>0)
		log_file << " WARNING in set_RhoUi_extern:" << endl
		     << " Requested Rho values not in extern simulation Box: " << num_values_not_in_extern_box << endl 
	             << " Setting extern_rho, v to zero." << endl;
	
	//! clean up
	delete[] extern_rho;

	if(serialize_reading_of_extern_field)
	{
		is_finished = 1;
		
		//! - every process must send except for last process
		//! - use Blocking Send
		if(serialize_writing_of_StateFile && mpi_myRank < mpi_num_processes-1)
		{	
			//! wait time in micro seconds
			usleep(TIME_TO_WAIT);
			mpi_myComm.Send(&is_finished, length_in_byte, MPI_INT32, mpi_myRank+1, tag);

		}
	}
	
	
	log_file << " done." << endl; 

}




void CHybrid::init_ion_neutral_variables(void)
{
	// NOTE global variable species_does_chemical_reactions
	
	//! check if species undergoes chemical reactions or collisions
	//! (used to speedup code)
	for(INT32 species=0; species<num_Particle_Species; species++)
	{
		D_REAL total_rate_of_species=0;
		
		for(INT32 neutral_species=0; neutral_species<num_Neutral_Species; neutral_species++)		
		{
			for(INT32 dest_species=0; dest_species<num_Particle_Species; dest_species++)
			total_rate_of_species += ReactionRate[species][dest_species][neutral_species];
		}
	
		total_rate_of_species += Recombination_for_Species[species];

		if(total_rate_of_species>0)
		species_does_chemical_reactions[species]=true;
		else
		species_does_chemical_reactions[species]=false;
		

		//! THE MOST UGLY PARAMETER EVER ;)
		for(INT32 neutral_species=0; neutral_species<NUM_NEUTRAL_SPECIES; neutral_species++)		
		obstacle_MP_weight_each_t0[species] += Global_IonProduction_Rate[species][neutral_species]*SI_t0/(SI_n0*SI_x0*SI_x0*SI_x0);
		
		//! How to insert a certain number of Ions each second:
		//! - Assume 1 particle in the Cell Centre of each cell with weight of 1.
		//! - Assume each cell to have a volume of 1 (=1*x0^3)
		//! - This does form the constant Background denity n0
		//!   -> a weight=1 particle represents n0*x0^3 particle

		//! Now inserting "obstacle_MP_weight_each_t0" means to insert
		//! obstacle_MP_weight_each_t0 *n0*x0^3 Particle.

		//! e.g. Mercury:
		//! B0 = 21nT
		//! n0 = 32.e6 1/m^3

		//! -> t0 ~ 0.5s
		//! -> x0 = 40.e3 m
		//! set obstacle_MP_weight_each_t0 = 1.

		//! ->   NP = 1. * 32.e6 * (40.e3)^3 = 2.048e21
		//! -> NP/s = 2.048e21 / 0.5s
		//! -> NP/s = 4.1e21/s

		//! In order to insert 1.e24 Ions/s
		//! obstacle_MP_weight_each_t0 = 244.
		//! is required
	
			
	}	
	
	if(use_velocity_dependent_reaction_rates)
		init_velocity_dependent_rates();
	
	
	for(INT32 neutral_species=0; neutral_species<num_Neutral_Species; neutral_species++)	
	{
		norm_IonProduction_Rate_fromNeutSpec[neutral_species]=0;
	
		for(INT32 species=0; species<num_Charged_Species; species++)
		norm_IonProduction_Rate_fromNeutSpec[neutral_species] += (PhotoionisationRate[species][neutral_species]+ElectronionisationRate[species][neutral_species])*SI_t0;
		
#ifdef use_neutral_species_as_field
		//! Beta of neutral gas
		Neutral_Betas[neutral_species]  = calcBeta*Neutral_Temperature[neutral_species];
#endif			
		
		Neutral_vth[neutral_species] = sqrt(Neutral_Temperature[neutral_species]*kB/(Neutral_Masses[neutral_species]*SI_m0))/SI_v0;

	
	}


}	


//!------------------------------------------------------------------------
//!- calc_electron_temperature:
//!------------------------------------------------------------------------
void CHybrid::calc_electron_temperature(INT32 Field_type)
{
    log_file << "  calc " << Field_Name[Field_type] << "...       ";
        
    set_zero_field(Field_type);
    
    //! Old approximation
//     for(int species=0; species<num_Inflow_Species; species++) 
//     {
//         add_multipliedField(Field_type,id_rhoSpecies1+species,Electron_Betas[species]/calcBeta); 
//     }
//     //! Teilen durch GesamtDichte gespeichert in id_rho_np1
//     Devide_fields(Field_type,Field_type,id_rho_np1,MCD_BField);
   
        
    //! Calc Electron Temperature by T_e= p_e/(n_e*k_b) 
    //! calcBeta = n_0*k_b/p_0 
    //! T_e = p_e^ast * p_0 / (n^ast * n_0 * k_b) = p_e^ast / (n^ast * calcBeta)
    add_multipliedField(Field_type,id_PEtotal,1./calcBeta);
    Devide_fields(Field_type,Field_type,id_rho_np1,MCD_BField);
    
    log_file << "done." << endl;
    
}

//!------------------------------------------------------------------------
//!- calc_RecombinationAlpha:
//!------------------------------------------------------------------------
void CHybrid::calc_RecombinationRate(void)
{
    
    bool does_any_species_recombine = false;
    for(INT32 species=0; species<num_Particle_Species; species++)
     if(Recombination_for_Species[species])
      does_any_species_recombine = true;
	
     if(does_any_species_recombine)
     {	     
	//! Calc electron Temperature 
	calc_electron_temperature(id_ElectronTemperature);
	
	calc_RecombinationAlphaField();
     }
     
}

//!--------------------------------------------------------------
//!- calc_RecombinationAlphaField:
//!--------------------------------------------------------------
void CHybrid::calc_RecombinationAlphaField(void)
{
    log_file << "  calc Recombination Rate ... ";
    
    //! set recombination field to zero
    //! to save time, this is done only for the recombination fields that are used
    //! and for all fields at the start of the simulation
    for(INT32 species=0; species<num_Particle_Species; species++)
     if(TL==1 || TL==TL_at_last_restore_state+1 || Recombination_for_Species[species])
      set_zero_field(id_recomb_Species1+species);
    
    for(INT32 level=0; level<= MAX_LEVEL; level++)
    {
        
        CBlock *temp_Block = BlockList_of_Lev[level];
        while(temp_Block)
        {
            
            if(mpi_myRank == temp_Block->responsible_mpi_process)
                temp_Block->calc_RecombinationAlphaField();
            
            temp_Block = temp_Block->next_Blk_of_BlockList;
        }
    }
    
    
    log_file << "done." << endl;
}


inline D_REAL sigma_temp(D_REAL mass, D_REAL Et, D_REAL ER, D_REAL v_ms, D_REAL exponent)
{
	D_REAL temp = ((0.5 * mass * m_p * v_ms*v_ms / (1.e3* e )) - Et )/(1.*ER);
	
	return pow(temp,exponent);
}

//! --------------------------------------------------------------------------
//! precalculate velocity-dependent part of reaction rates
//! because it is too time-consuming during run-time
//! --------------------------------------------------------------------------
void CHybrid::init_velocity_dependent_rates(void)
{
	memset(velocity_rate,0,30*NUM_PARTICLE_SPECIES*sizeof(D_REAL));
	//! velocity is supposed to be in km/s
	
	
	for(INT32 species=0;species<num_Particle_Species;species++)
	{
		
		//! velocity dependent part of reaction
		//! O+ + H2O -> H2O+ + O
		if(species==0)
		{
			for(INT32 v=0;v<30;v++)
			{	
				D_REAL vnorm = v*1e+03/SI_v0;
				
				const D_REAL A2 = 9e-08;
				const D_REAL B2 = 1.2e-08;	
				
				const D_REAL mu2 = 16.*18./(16.+18.)*m_p;	
				
				const D_REAL v_cms = vnorm*SI_v0*1.e+02;
				const D_REAL v_ms = vnorm*SI_v0;

				D_REAL sigma_1  = (A2 - B2*log10(mu2*v_ms*v_ms/e) )*(A2 - B2*log10(mu2*v_ms*v_ms/e) )
						- 1.5e-17*log10(v_cms);
				
				const D_REAL CoM_to_Lab = 1.*mu2/(16.*m_p);
				
				const D_REAL a1 = 10e9;
				const D_REAL a2 = 1.7;
				const D_REAL a3 = 0.001;
				const D_REAL a4 = 0.441;
				const D_REAL a5 = 5.86;
				const D_REAL a6 = 3.82;
				const D_REAL a7 = 3.2e-4;
				const D_REAL a8 = 10;
				const D_REAL ER = 25;
				const D_REAL Et = -0.001;
				
				const D_REAL mu4 = 16.*18./(16.+18.)*m_p;	
												
				D_REAL sigma_10;
				
				sigma_10 =  a1*sigma_temp(16,Et,ER,v_ms,a2)/
					( 1+sigma_temp(16,Et,a3,v_ms,a2+a4) + sigma_temp(16,Et,a5,v_ms,a2+a6) );
					+ a7* a1* sigma_temp(1,Et,a8/ER,v_ms,a2)/
					( 1+sigma_temp(1,Et,a8/a3,v_ms,a2+a4) + sigma_temp(1,Et,a8/a5,v_ms,a2+a6) );
				
				sigma_10 *= 1.e-16;	
					
 				D_REAL sigma = 1./3.*(sigma_10 + 2.*sigma_1);
			  
				velocity_rate[v][species] = v_cms*2.*sigma*CoM_to_Lab;
			  			
			}	
			species_does_chemical_reactions[species]=true;
		}
		
		
		//! velocity dependent part of reaction
		//! OH+ + H2O -> H2O+ + OH
		if(species==1)
		{
			for(INT32 v=0;v<30;v++)
			{	
				D_REAL vnorm = v*1e+03/SI_v0;
				
				const D_REAL A2 = 9e-08;
				const D_REAL B2 = 1.2e-08;	
				
				const D_REAL mu2 = 17.*18./(17.+18.)*m_p;	
				
				const D_REAL v_cms = vnorm*SI_v0*1.e+02;
				const D_REAL v_ms = vnorm*SI_v0;

				D_REAL sigma_1  = (A2 - B2*log10(mu2*v_ms*v_ms/e) )*(A2 - B2*log10(mu2*v_ms*v_ms/e) )
						- 1.5e-17*log10(v_cms);
				
				const D_REAL CoM_to_Lab = 1.*mu2/(17.*m_p);
				
				const D_REAL a1 = 10e9;
				const D_REAL a2 = 1.7;
				const D_REAL a3 = 0.001;
				const D_REAL a4 = 0.441;
				const D_REAL a5 = 5.86;
				const D_REAL a6 = 3.82;
				const D_REAL a7 = 3.2e-4;
				const D_REAL a8 = 10;
				const D_REAL ER = 25;
				const D_REAL Et = -0.001;
				
				const D_REAL mu4 = 17.*18./(17.+18.)*m_p;	
												
				D_REAL sigma_10;
				
				sigma_10 =  a1*sigma_temp(17,Et,ER,v_ms,a2)/
					( 1+sigma_temp(17,Et,a3,v_ms,a2+a4) + sigma_temp(17,Et,a5,v_ms,a2+a6) );
					+ a7* a1* sigma_temp(1,Et,a8/ER,v_ms,a2)/
					( 1+sigma_temp(1,Et,a8/a3,v_ms,a2+a4) + sigma_temp(1,Et,a8/a5,v_ms,a2+a6) );
				
				sigma_10 *= 1.e-16;	
					
 				D_REAL sigma = 1./3.*(sigma_10 + 2.*sigma_1);
			  
				velocity_rate[v][species] = v_cms*2.*sigma*CoM_to_Lab;
				
			 }
 			species_does_chemical_reactions[species]=true;
		}
		
		
		
		
		//! velocity dependent part of reaction
		//! H2O+ + H2O -> H2O+ + H2O
		if(species==2)
		{
			for(INT32 v=0;v<30;v++)
			{	
				D_REAL vnorm = v*1e+03/SI_v0;
				
				const D_REAL A2 = 4.3e-08;
				const D_REAL B2 = 9.5e-09;
				const D_REAL mu2 = 18*18/(18+18)*m_p;
				
				const D_REAL v_ms = vnorm*SI_v0;

				const D_REAL v_cms = v_ms*1.e+02;

				D_REAL sigma_1  = (A2 - B2*log10(mu2*v_ms*v_ms/e) );
				
				const D_REAL CoM_to_Lab = 1.*mu2/(18.*m_p);
			  
				velocity_rate[v][species] = v_cms*2.*sigma_1*sigma_1*CoM_to_Lab;
	
			}
			species_does_chemical_reactions[species]=true;
		}	
				
		//! velocity dependent part of reaction
		//! H+ + H2O -> H2O+ + H
		if(species==4)
		{
			for(INT32 v=0;v<30;v++)
			{
				D_REAL vnorm = v*1e+03/SI_v0;
				
				const D_REAL a1 = 5.85e5;
				const D_REAL a2 = 1.6;
				const D_REAL a3 = 0.1;
				const D_REAL a4 = 0.441;
				const D_REAL a5 = 5.86;
				const D_REAL a6 = 3.82;
				const D_REAL a7 = 3.2e-4;
				const D_REAL a8 = 10;
				const D_REAL ER = 25;
				const D_REAL Et = -0.001;
				
				const D_REAL mu4 = 1.*18/(1.+18.)*m_p;	
				
				const D_REAL v_ms = vnorm*SI_v0;

				const D_REAL v_cms = v_ms*1.e+02;

								
				D_REAL sigma_40;
				
				sigma_40 =  a1*sigma_temp(1,Et,ER,v_ms,a2)/
					( 1+sigma_temp(1,Et,a3,v_ms,a2+a4) + sigma_temp(1,Et,a5,v_ms,a2+a6) );
					+ a7* a1* sigma_temp(1,Et,a8/ER,v_ms,a2)/
					( 1+sigma_temp(1,Et,a8/a3,v_ms,a2+a4) + sigma_temp(1,Et,a8/a5,v_ms,a2+a6) );
				
				sigma_40 *= 1.e-16;	
					
				const D_REAL CoM_to_Lab = 1.*mu4/(1.*m_p);
			  
				velocity_rate[v][species] = v_cms*2.*sigma_40*CoM_to_Lab;
				
			}
			species_does_chemical_reactions[species]=true;
		}	
		
	}		
}	

