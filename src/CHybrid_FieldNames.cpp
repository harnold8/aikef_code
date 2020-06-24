


#include <math.h>
#include "CHybrid.h"
#include "parameters.h"
#include "absolute_Globals.h"

#include <iostream>
#include <fstream>



//!--------------------------------------------------------
//!- set_Field_Comps_Names_IDs: 
//!--------------------------------------------------------
void CHybrid::set_Field_Names_IDs(void)
{

      //!--------------------------------------------
      //!------ NAMES -------------------------------
      //!--------------------------------------------
      //! Field_Name,Field_Name_Size defined in absolute_Globals.h 


	//! Max FieldNameSize
	//! -> (adjust in Visu on change !!!)
	Field_Name_Size = 50;
	log_file<<"gecko0"<<endl;
	Field_Name = new char*[NUM_FIELDS];
	for(INT32 i=0; i<NUM_FIELDS; i++)
	{
		Field_Name[i] = new char[Field_Name_Size];
		memset(Field_Name[i],0,Field_Name_Size*sizeof(char));
	}
	log_file<<"gecko1"<<endl;
	//! is now set in parameters.cpp
// 	extern_Field_Path = new char*[num_externRhoVelocityFields];
// 	for(INT32 i=0; i<num_externRhoVelocityFields; i++)
// 	{
// 		extern_Field_Path[i] = new char[200];
// 		memset(extern_Field_Path[i], 0, 200*sizeof(char));
// 
// 		sprintf(extern_Field_Path[i], "/velocity_ionization_field.txt");
// 	}



	log_file<<"gecko1"<<endl;
	sprintf(Field_Name[       id_BTotal], "B");
	sprintf(Field_Name[        id_BEven], "BEven");
	sprintf(Field_Name[       id_EField], "E");
	sprintf(Field_Name[      id_UI_plus], "U");
	sprintf(Field_Name[     id_UI_minus], "Umin");
	sprintf(Field_Name[      id_rho_np1], "Rho");
	sprintf(Field_Name[          id_Eta], "Eta");
	sprintf(Field_Name[         id_rotB], "curlB");
	sprintf(Field_Name[         id_divB], "divB");
	sprintf(Field_Name[         id_divrotB], "divrotB");
	sprintf(Field_Name[         id_divU], "divU");
	sprintf(Field_Name[         id_gradPE], "gradPE");
	sprintf(Field_Name[	   id_gradB], "gradB");
// 	sprintf(Field_Name[id_Refine_Rating], "Ref_Rating");
	sprintf(Field_Name[        id_PhiDC], "PhiDC");
	sprintf(Field_Name[    id_PMagnetic], "B_pressure");
	sprintf(Field_Name[       id_PTotal], "total_pressure");
	sprintf(Field_Name[id_FieldTime    ], "Time_Field");
	sprintf(Field_Name[id_ParticleTime ], "Time_Particle");
        sprintf(Field_Name[id_gradPI0      ], "gradTI0");
        sprintf(Field_Name[id_Refine_Rating], "MPiC");
        sprintf(Field_Name[id_ElectronTemperature], "ElectronTemp");
        sprintf(Field_Name[id_rho_np1_recombined], "recombRho");
        log_file<<"gecko2"<<endl;
	sprintf(Field_Name[id_scratch_scalar ], "LAM");
	sprintf(Field_Name[id_scratch_vector ], "GAM");
        sprintf(Field_Name[         id_divE], "divE");
	log_file<<"gecko"<<endl;

	for(int species=0; species<num_Charged_Species; species++)
	{
		sprintf(Field_Name[id_UI_Species1    +species], "s%d_u"  ,species);
		sprintf(Field_Name[id_rhoSpecies1    +species], "s%d_rho",species);
                sprintf(Field_Name[id_gyro_Species1  +species], "s%d_gyroradius",species);
                sprintf(Field_Name[id_gyro_el_Species1  +species], "s%d_electron_gyroradius",species);
                sprintf(Field_Name[id_recomb_Species1  +species], "s%d_electron_gyroradius",species);
		sprintf(Field_Name[id_ForceSpecies1  +species], "s%d_force",species);
		sprintf(Field_Name[id_PISpecies1     +species], "s%d_vth2" ,species);
		sprintf(Field_Name[id_PESpecies1     +species], "s%d_Pe" ,species);

	}

#if defined(nonadiabatic_gradPE_TERM) || defined(use_neutral_species_as_field)
	for(int species=0; species<num_Neutral_Species; ++species)
	{
		sprintf(Field_Name[id_numberdensity_neutralSpecies1     + species], "ns%d_numberdensity", species);
		sprintf(Field_Name[id_velocity_neutralSpecies1          + species], "ns%d_velocity", species);
		sprintf(Field_Name[id_pressure_neutralSpecies1          + species], "ns%d_pressure", species);
		sprintf(Field_Name[id_new_electron_beta_neutralSpecies1 + species], "ns%d_beta_new", species);

	}
#endif

#if defined(use_dust_species_as_field)
	for(int species=0; species<num_Dust_Species; ++species)
	{
		sprintf(Field_Name[id_density_dustSpecies1     + species], "dust%d_density", species);
		sprintf(Field_Name[id_velocity_dustSpecies1          + species], "dust%d_velocity", species);
	}
#endif

#if defined(use_ion_production_as_field)
	for(int species=0; species<num_ion_prod_fields; ++species)
	{
		sprintf(Field_Name[id_density_ionProdSpecies1     + species], "ionProd%d_density", species);
		sprintf(Field_Name[id_velocity_ionProdSpecies1          + species], "ionProd%d_velocity", species);
	}
#endif

	sprintf(Field_Name[id_PEtotal ], "Pe");
        
	//! names for RHO & U  different from the default names
// 	sprintf(Field_Name[id_rhoSpecies1 +0], "sw_rho");
// 	sprintf(Field_Name[id_UI_Species1 +0], "sw_u");

// 	if(num_Charged_Species==2)
// 	{
// 		sprintf(Field_Name[id_rhoSpecies1 +1], "hi_rho");
// 		sprintf(Field_Name[id_UI_Species1 +1], "hi_u");
// 	}

	for(INT32 rho_extern=0; rho_extern<num_externRhoVelocityFields; rho_extern++)
	{
		sprintf(Field_Name[id_externRho1 +rho_extern], "RhoExtern%d", rho_extern);
		sprintf(Field_Name[id_extern_Ui1 +rho_extern], "VelExtern%d", rho_extern);

	}

	for(INT32 average_field=0; average_field<num_average_fields; average_field++) {
                sprintf(Field_Name[id_average_Field1 +average_field], "%s_%s",average_field_name_prefix, Field_Name[IDs_of_Fields_to_average[average_field]]);
        }



      //!--------------------------------------------
      //!------ VISU_NR -----------------------------
      //!--------------------------------------------
      //! VisuNr_FType defined in absolut_Globals.h

	//! define in which order fields shall appear in visualization
	VisuNr_FType = new INT32[NUM_FIELDS];

	VisuNr_FType[id_BTotal ] = 0;
	VisuNr_FType[id_EField ] = 1;

	VisuNr_FType[id_Eta    ] = 2;
	VisuNr_FType[id_rotB   ] = 3;
	VisuNr_FType[id_divB   ] = 4;
	VisuNr_FType[id_divE   ] = 5;
	VisuNr_FType[id_PMagnetic ] =  6;

	//! total fields
	VisuNr_FType[id_UI_plus ] = 7;
	VisuNr_FType[id_UI_minus] = 7;
	VisuNr_FType[id_rho_np1] = 8;
	VisuNr_FType[id_gradB] =  9;
	VisuNr_FType[id_PTotal       ] = 10;

	VisuNr_FType[id_FieldTime       ] = 11;
	VisuNr_FType[id_ParticleTime    ] = 12;
	VisuNr_FType[id_scratch_vector       ] = 13;
	VisuNr_FType[id_scratch_scalar        ] = 14;
        VisuNr_FType[id_rho_np1_recombined ] = 15;


	for(int species=0; species<num_Charged_Species; species++)
	{
		VisuNr_FType[id_UI_Species1 +species] =  16 + 5*species;
// 		log_file << " Species"<<species<<" velocity: "<< VisuNr_FType[id_UI_Species1 +species] << endl;

		VisuNr_FType[id_rhoSpecies1 +species] = VisuNr_FType[id_UI_Species1]    +1 +5*species;
// 		log_file << " Species"<<species<<" density: "<< VisuNr_FType[id_rhoSpecies1 +species] << endl;

		VisuNr_FType[id_ForceSpecies1  +species] = VisuNr_FType[id_UI_Species1] +2 +5*species;
// 		log_file << " Species"<<species<<" temperature: "<< VisuNr_FType[id_PISpecies1 +species] << endl;

		VisuNr_FType[id_PISpecies1  +species] = VisuNr_FType[id_UI_Species1]    +3 +5*species;
// 		log_file << " Species"<<species<<" temperature: "<< VisuNr_FType[id_PISpecies1 +species] << endl;

		VisuNr_FType[id_PESpecies1  +species] = VisuNr_FType[id_UI_Species1]    +4 +5*species;
// 		log_file << " Species"<<species<<" temperature: "<< VisuNr_FType[id_PISpecies1 +species] << endl;
		
	}

	for(INT32 rho_extern=0; rho_extern<num_externRhoVelocityFields; rho_extern++)
	{
		VisuNr_FType[id_externRho1 +rho_extern] = VisuNr_FType[id_PESpecies1 +num_Charged_Species-1]+1 +2*rho_extern;
// 		log_file << "ExSpecies"<<rho_extern<<" density:  "<<VisuNr_FType[id_externRho1 +rho_extern] << endl;

		VisuNr_FType[id_extern_Ui1 +rho_extern] = VisuNr_FType[id_PESpecies1 +num_Charged_Species-1]+2 +2*rho_extern;
// 		log_file << "ExSpecies"<<rho_extern<<" velocity: "<<VisuNr_FType[id_extern_Ui1 +rho_extern] << endl;
	}


        
}
