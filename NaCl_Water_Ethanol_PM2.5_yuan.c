/*****************************************
Copyright Reserved

writtern by Dr. Xiaole Chen at Southeast University in China
Email: xiaole_chennj@sina.com

Please keep this UDF within
Prof. Clement Kleinstreuer's Computational Multi-Physics Lab at North Carolina State University
Dr. Yu Feng's Computational Biofluidics and Biomechanics Lab at Oklahoma State University

*****************************************/



/*****************************************
IMPORTANT NOTE

Do use double presicion!

Single presicion would have error since small value
likes 298.15+1e-7 would still equal to 298.15

"Species List Sequence" in fluent should be:  Air   WaterVapor
*****************************************/
#include "udf.h"
#include "dpm.h"
#include "pdf_props.h"
#include "pdf_table.h"
#include "mem.h"
#include "materials.h"
#include "stdio.h"
#include "math.h"
#include "RH_T_P_201607.h"
//#include "multi_component_droplet.h"


long int NumberMark = 0;

/*******************************
EDIT THE FOLLOWING VARIABLES

diameter	1.00E-06			1.00E-09			5.00E-09			2.00E-08			1.00E-07			1.00E-05
water		3.8056890008E-16	3.8056890008E-25	4.7571112510E-23	3.0445512006E-21	3.8056890008E-19	3.8056890008E-13
alcohol		9.5142225019E-17	9.5142225019E-26	1.1892778127E-23	7.6113780015E-22	9.5142225019E-20	9.5142225019E-14
NaCl		9.5142225019E-17	9.5142225019E-26	1.1892778127E-23	7.6113780015E-22	9.5142225019E-20	9.5142225019E-14
fluorescein	2.3785556255E-18	2.3785556255E-27	2.9731945319E-25	1.9028445004E-23	2.3785556255E-21	2.3785556255E-15

********************************/
double ParticleInitDimater = 1e-9;

double TIMESTEP = 1.0e-7;

//double ParticleInitDimater = 1.0e-6;
/*************
Components' mass is calculated based on PM1 particle
**************/
double InitMassWater = 3.8056890008E-25;
double InitMassAlcohol = 9.5142225019E-26;
double InitMassNaCl = 9.5142225019E-26;
double InitMassFluorescein = 2.3785556255E-27;
double InitDropletTemp = 310.15;



int InitScalarNumForFF = 10;

/*
The list of variables defined in this UDF


Water.Density = WATERDENSITY;
Fluorescein.Density = FLUORESCEINDENSITY;
NaCl.Density = NACLDENSITY;
Alcohol.Density = ALCOHOLDENSITY;

P_USER_REAL(p, 0) = Fluorescein.Mass;
P_USER_REAL(p, 1) = NaCl.Mass;
P_USER_REAL(p, 2) = Water.Mass;
P_USER_REAL(p, 3) = Alcohol.Mass;
P_USER_REAL(p, 4) = NaCl.MoleFraction;
P_USER_REAL(p, 5) = Water.MoleFraction;
P_USER_REAL(p, 6) = Alcohol.MoleFraction;




*/


/*********************
Components Properties
*********************/

//#define OLEICACIDDENSITY 895.0


#define MolarMassNaCl 0.05844






typedef struct Component{
	double Density;
	double Mass;
	double MassLoss;
	double MassFraction;
	double MoleFraction;
	double MassFractionInCell;
	double MoleFractionInCell;
	double PartialPressureInCell;
	double MassFlux;
	double SchmidtNumber;
	double SpeciesDiffusivity;
	double SherwoodNumber;
	double MassFractionOnSurf;
	double SaturatedPressure;
	double KelvinEffectFactor;
	double ActivityCoefficient;
}EvaporativeSolute;

typedef struct Component1{
	double Density;
	double Mass;
	double MassFraction;
	double MoleFraction;

}NonEvaporativeSolute;

typedef struct Component2{
	double Density;
	double Mass;

}NonEvaporativeCore;

double particle_surface_tension(double GlycerolMass, double WaterMass)
{

	return 0.07199;

}

double mixture_density(double AlcoholMass, double AlcoholDensity,
	double WaterMass, double WaterDensity,
	double NaClMass, double NaClDensity,
	double FluoresceinMass, double FluoresceinDensity)
{

	/*************************
	molbiol.ru/eng/protocol/01_22.html

	NaCl+Water Solution density
	y=995.57+7.64*x
	y=density

	x = NaCl/Water*100.0;
	**************************/
	double WeightRatio = NaClMass / WaterMass*100.0;
	double Volume = 0.0;
	double NaCl1 = 0.0, NaCl2 = 0.0;
	if (WeightRatio > 35.9)
	{
		//NaCl begins to crystallize
		WeightRatio = 35.9;
		NaCl1 = WaterMass / 100.0*35.9;
		NaCl2 = NaClMass - NaCl1;
		Volume = FluoresceinMass / FluoresceinDensity + AlcoholMass / AlcoholDensity
			+ (NaCl1 + WaterMass) / (995.57 + 7.64*WeightRatio) + NaCl2 / NaClDensity;

	}
	else
	{
		Volume = FluoresceinMass / FluoresceinDensity + AlcoholMass / AlcoholDensity + (NaClMass + WaterMass) / (995.57 + 7.64*WeightRatio);


	}
	return ((AlcoholMass + WaterMass + FluoresceinMass + NaClMass) / Volume);

}





DEFINE_DPM_LAW(MultiComponentDroplet_NaCl_Water_Alcohol, p, ci)
{

	Thread *cc;
	cell_t cell;
	real RelVel[ND_ND], FluidVel[ND_ND], ParticleVel[ND_ND];
	double ParticleReynoldsNo;
	double ParticleVolume;
	double MFVaporinCell, MFAirinCell;
	double MoleFractionVaporinCell;
	double InitialVolumeRatio;
	double ParticleSherwoodTerm;
	double KnudsenNumber;
	double Cm, Ct;
	double ParticleTimeStep;
	double MoleInSolution;
	double CellGasDensity, CellTemperature;
	double ParticleDiameter;
	double ParticleSurfaceTension;
	double DTemperatureDt;
	double NusseltNumber;
	double CurrentTime;
	double TEMP;								// only used for debug

	double DryAirDensity;	
	double VaporDensity;		//value of 25 oC
	double VaporPressure;
	//0.023100513 for 25.0292 oC
	double TotalPressure;
	double DryAirMassFraction;
	double DryAirMoleFraction;
	EvaporativeSolute Water = { 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	EvaporativeSolute Alcohol = { 785.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	NonEvaporativeSolute NaCl = { 2165.0, 0.0, 0.0, 0.0 };
	NonEvaporativeCore Fluorescein = { 1601.0, 0.0 };
	//need to tell when Alcohol is all vaporated.



	cell = P_CELL(p);
	cc = P_CELL_THREAD(p);
	DryAirDensity = 101325 / 287.058 / C_T(cell,cc);
	VaporDensity = vapor_density(C_T(cell, cc));
	VaporPressure = VAPOR_PRESSURE(C_T(cell, cc));
	/*******************



	//	ParticleTimeStep = P_DT(p);
	//	CurrentTime = P_TIME(p);
	//	ParticleDiameter = P_DIAM(p);
	*************/

	NumberMark++;

	/**************************************
	Obtain the parameters in the cell
	**************************************/
	CellGasDensity = C_R(cell, cc);
//	Message("mark0 \n");
	Water.Density = WATERDENSITY;
	Fluorescein.Density = FLUORESCEINDENSITY;
	NaCl.Density = NACLDENSITY;
	Alcohol.Density = ALCOHOLDENSITY;
//	Message("mark1 \n");
	TotalPressure = C_P(cell, cc); // total pressure in the cell

	DryAirMassFraction = C_YI(cell, cc, 0);
	Water.MassFractionInCell = C_YI(cell, cc, 1);
//	Message("Air MassF %e, Water MassF %e, Alochol %e \n", C_YI(cell, cc, 0), C_YI(cell, cc, 1), (1 - C_YI(cell, cc, 1) - C_YI(cell, cc, 0)));
	Water.MoleFractionInCell = (Water.MassFractionInCell / MolarMassWater) 
		/ (Water.MassFractionInCell / MolarMassWater + DryAirMassFraction / MolarMassAir);
	DryAirMoleFraction = 1 - Water.MoleFractionInCell;
	Alcohol.MassFractionInCell = 0.0;
	Alcohol.MoleFractionInCell = 0.0;

	Water.PartialPressureInCell = TotalPressure*Water.MoleFractionInCell;
	Alcohol.PartialPressureInCell = 0.0;
	Water.SpeciesDiffusivity = 2.42e-5;
	Alcohol.SpeciesDiffusivity = 1.15e-5;

	if (P_USER_REAL(p, 1) == 0.0)
	{
//		Message("mark2 \n");
		/****************************
		initial condition:
		P_USER_REAL(p, 0) = Fluorescein.Mass;
		P_USER_REAL(p, 1) = NaCl.Mass;
		P_USER_REAL(p, 2) = Water.Mass;
		P_USER_REAL(p, 3) = Alcohol.Mass;
		P_USER_REAL(p, 4) = NaCl.MoleFraction;
		P_USER_REAL(p, 5) = Water.MoleFraction;
		P_USER_REAL(p, 6) = Alcohol.MoleFraction;

		P_USER_REAL(p, 8) = Water.MassLoss;
		P_USER_REAL(p, 9) = Alcohol.MassLoss;
		P_USER_REAL(p, 10) = Water.MassFractionOnSurf;
		P_USER_REAL(p, 11) = Alcohol.MassFractionOnSurf;
		P_USER_REAL(p, 12) = Water.MassFractionInCell;
		P_USER_REAL(p, 13) = Water.KelvinEffectFactor;

		****************************/
		P_DIAM(p) = ParticleInitDimater;
		//InitialVolumeRatio = pow((P_DIAM(p)/1.e-6), 3.0);
		P_USER_REAL(p, 0) = InitMassFluorescein;
		P_USER_REAL(p, 1) = InitMassNaCl;
		P_USER_REAL(p, 2) = InitMassWater;
		P_USER_REAL(p, 3) = InitMassAlcohol;


		InitialVolumeRatio = PI*pow(P_DIAM(p), 3.0) / 6.;

		P_RHO(p) = (P_USER_REAL(p, 0) + P_USER_REAL(p, 1) + P_USER_REAL(p, 2) + P_USER_REAL(p, 3))/InitialVolumeRatio;
		P_T(p) = InitDropletTemp;

		NumberMark = 0;

	}
	Fluorescein.Mass = P_USER_REAL(p, 0);
	NaCl.Mass = P_USER_REAL(p, 1);
	Water.Mass = P_USER_REAL(p, 2);
	Alcohol.Mass = P_USER_REAL(p, 3);

	NaCl.MassFraction = P_USER_REAL(p, 1) + P_USER_REAL(p, 2) + P_USER_REAL(p, 3);
	Water.MassFraction = P_USER_REAL(p, 2) / NaCl.MassFraction;
	Alcohol.MassFraction = P_USER_REAL(p, 3) / NaCl.MassFraction;
	NaCl.MassFraction = P_USER_REAL(p, 1) / NaCl.MassFraction;

	if(TP_USER_REAL(p, 12)<200)
	{
		
		
	//Message("1:%e\n",P_DIAM(p));
	}
	
	if (NaCl.MassFraction == 1.0)
	{
		NaCl.MoleFraction = 1.0;
		Water.MoleFraction = 0.0;
		Alcohol.MoleFraction = 0.0;
	}
	else
	{
		if ((NaCl.Mass / Water.Mass) >= 0.359)
		{
			MoleInSolution = Water.Mass*0.359 / MolarMassNaCl + Alcohol.Mass / MolarMassAlcohol + Water.Mass / MolarMassWater;
			NaCl.MoleFraction = Water.Mass*0.359 / MolarMassNaCl / MoleInSolution;

		}
		else
		{
			MoleInSolution = NaCl.Mass / MolarMassNaCl + Alcohol.Mass / MolarMassAlcohol + Water.Mass / MolarMassWater;
			NaCl.MoleFraction = P_USER_REAL(p, 1) / MolarMassNaCl / MoleInSolution;
		}
		Water.MoleFraction = P_USER_REAL(p, 2) / MolarMassWater / MoleInSolution;
		Alcohol.MoleFraction = P_USER_REAL(p, 3) / MolarMassAlcohol / MoleInSolution;
	}

//	Message("0. MASS NaCl= %e, Alcohol= %e, Water= %e \n", NaCl.Mass, Alcohol.Mass, Water.Mass);
//	Message("1. %e, %e, %e,%e \n", NaCl.MassFraction, Water.MassFraction, Alcohol.MassFraction, MoleInSolution);


	P_USER_REAL(p, 4) = NaCl.MoleFraction;
	P_USER_REAL(p, 5) = Water.MoleFraction;
	P_USER_REAL(p, 6) = Alcohol.MoleFraction;

	Water.SchmidtNumber = C_MU_L(cell, cc) / C_R(cell, cc) / Water.SpeciesDiffusivity;
	Alcohol.SchmidtNumber = C_MU_L(cell, cc) / C_R(cell, cc) / Alcohol.SpeciesDiffusivity;



	/***** Particle Reynolds Number *****/
	FluidVel[0] = C_U(cell, cc);
	FluidVel[1] = C_V(cell, cc);
	FluidVel[2] = C_W(cell, cc);
	ParticleVel[0] = P_VEL(p)[0];
	ParticleVel[1] = P_VEL(p)[1];
	ParticleVel[2] = P_VEL(p)[2];
	NV_VV(RelVel, = , ParticleVel, -, FluidVel);
	ParticleReynoldsNo = C_R(cell, cc)*P_DIAM(p)*NV_MAG(RelVel) / C_MU_L(cell, cc);

//	ParticleReynoldsNo = 10.0;  //only used this line for test
	/******** Sherwood Number *****************/
	Water.SherwoodNumber = ParticleReynoldsNo;
	if (Water.SherwoodNumber > 1)
	{
		Water.SherwoodNumber = pow(Water.SherwoodNumber, 0.0777);

	}
	Alcohol.SherwoodNumber = Water.SherwoodNumber;
	Water.SherwoodNumber = 1 + pow(1 + ParticleReynoldsNo*Water.SchmidtNumber, (1.0 / 3.0))
		*Water.SherwoodNumber;
	Alcohol.SherwoodNumber = 1 + pow(1 + ParticleReynoldsNo*Alcohol.SchmidtNumber, (1.0 / 3.0))
		*Alcohol.SherwoodNumber;



	/*** KnudsenNumberCorrection   Cm  Ct ***/
	KnudsenNumber = 2.0* MeanFreePath / P_DIAM(p);
	Cm = (1.0 + KnudsenNumber) / (1.0 + (4.0 / 3.0 + 0.377)*KnudsenNumber
		+ 4.0 / 3.0*KnudsenNumber*KnudsenNumber);
	Ct = Cm;



	/*** KelvinEffectFactor  ***/
	ParticleSurfaceTension = 0.07199;  //used the value of water   does not affact the kelvin effect factor much
	Water.KelvinEffectFactor = exp(4.0*ParticleSurfaceTension*MolarMassWater
		/ UniversalGasConstant / P_RHO(p) / P_DIAM(p) / P_T(p));

	Alcohol.KelvinEffectFactor = exp(4.0*ParticleSurfaceTension*MolarMassAlcohol
		/ UniversalGasConstant / P_RHO(p) / P_DIAM(p) / P_T(p));

	/***** Mass Fraction on Surface  ****/

	Water.ActivityCoefficient = 1.0;
	Alcohol.ActivityCoefficient = 1.0;

//	vapor_saturated_pressure(P_T(p));




	/*************
	Saturated Alochol pressure
	pow(10.0, (8.20417 - 1642.89 / (230.3 + P_T(p)-273.15)))*(101325.0/760.0)
	from
	http://ddbonline.ddbst.de/AntoineCalculation/AntoineCalculationCGI.exe?component=Ethanol
	
	*************/
	//Saturated Alochol pressure    P_USER_REAL(p, 7) 
	P_USER_REAL(p, 7) = pow(10.0, (8.20417 - 1642.89 / (230.3 + P_T(p) - 273.15)))*(101325.0 / 760.0);  

	Alcohol.MassFractionOnSurf = Alcohol.ActivityCoefficient*Alcohol.MoleFraction
		*P_USER_REAL(p, 7) * Alcohol.KelvinEffectFactor / C_R(cell, cc) / UniversalGasConstant
		*MolarMassAlcohol / P_T(p);

	Water.ActivityCoefficient = 1.0 - (Alcohol.MoleFraction + 1.85*NaCl.MoleFraction) / Water.MoleFraction;
	if (Water.ActivityCoefficient < 0)
	{
		Water.ActivityCoefficient = Water.MoleFraction / (Alcohol.MoleFraction + 1.85*NaCl.MoleFraction + Water.MoleFraction);
	}

	Water.MassFractionOnSurf = Water.ActivityCoefficient 
		* VAPOR_PRESSURE(P_T(p)) * Water.KelvinEffectFactor / C_R(cell, cc) / UniversalGasConstant
		*MolarMassWater / P_T(p);
//	Message("Water.MoleFraction %e, Mass Fraction %e,Water.KelvinEffectFactor %e,Cell Density%e, VAPOR_PRESSURE %e\n", Water.MoleFraction, VAPOR_PRESSURE(P_T(p)) / CellGasDensity / UniversalGasConstant*MolarMassWater / P_T(p), Water.KelvinEffectFactor, C_R(cell, cc), VAPOR_PRESSURE(P_T(p)));
//	Message(" Alcohol  Water  MassFractionOnSurf %e,%e \n", Alcohol.MassFractionOnSurf, Water.MassFractionOnSurf);

//	Message("%e,%e,%e,%e,%e,%e,%e,%e \n", Water.ActivityCoefficient, Water.MoleFraction, VAPOR_PRESSURE(P_T(p)), Water.KelvinEffectFactor, C_R(cell, cc), UniversalGasConstant, MolarMassWater, P_T(p));
//	Message("%e,%e,%e,%e,%e,%e,%e,%e \n", Alcohol.ActivityCoefficient, Alcohol.MoleFraction, P_USER_REAL(p, 7), Alcohol.KelvinEffectFactor, C_R(cell, cc), UniversalGasConstant, MolarMassWater, P_T(p));

	/*********Mass Flux**********/
	Water.MassFlux = CellGasDensity * Water.SherwoodNumber * Water.SpeciesDiffusivity * Cm / P_DIAM(p)
		*log((1.0 - Water.MassFractionInCell) / (1.0 - Water.MassFractionOnSurf));
	Alcohol.MassFlux = CellGasDensity * Alcohol.SherwoodNumber * Alcohol.SpeciesDiffusivity * Cm / P_DIAM(p)
		*log((1.0 - Alcohol.MassFractionInCell) / (1.0 - Alcohol.MassFractionOnSurf));

	Water.MassLoss = 4.0*PI*P_DIAM(p)*P_DIAM(p)*Water.MassFlux*P_DT(p);
	Alcohol.MassLoss = 4.0*PI*P_DIAM(p)*P_DIAM(p)*Alcohol.MassFlux*P_DT(p);


	/*** Mass ****/
	if ((Water.Mass - Water.MassLoss) > 0.0){ Water.Mass -= Water.MassLoss; }
	else{ Water.Mass = 0.0; }
	if ((Alcohol.Mass - Alcohol.MassLoss) > 0.0){ Alcohol.Mass -= Alcohol.MassLoss; }
	else{ Alcohol.Mass = 0.0; }

	/*** Density, Volume and Diameter ***/

	P_RHO(p) = mixture_density(Alcohol.Mass, ALCOHOLDENSITY, Water.Mass, WATERDENSITY, NaCl.Mass, NACLDENSITY,
		Fluorescein.Mass, FLUORESCEINDENSITY);

	ParticleVolume = (Alcohol.Mass + Water.Mass + NaCl.Mass + Fluorescein.Mass) / P_RHO(p);
	P_DIAM(p) = pow(6.0*ParticleVolume / PI, 1.0 / 3.0);


	P_USER_REAL(p, 2) = Water.Mass;
	P_USER_REAL(p, 3) = Alcohol.Mass;
	
	if(TP_USER_REAL(p, 12)<200)
	{
		
		
	Message("%e\n",P_DIAM(p));
	}
 TP_USER_REAL(p, 12)+=1;

	/*******   Temperature Calculation *********/

//	ParticleSherwoodTerm = P_USER_REAL(p, InitScalarNumForFF + 30);
//	if (ParticleSherwoodTerm > 1)ParticleSherwoodTerm = pow(ParticleSherwoodTerm, 0.0777);

//	NusseltNumber = pow((1 + P_USER_REAL(p, InitScalarNumForFF + 30)*PrandtlNumberOfAir), 1.0 / 3.0)*ParticleSherwoodTerm;

	DTemperatureDt = - Water.MassLoss*WaterEvaLatentHeat - Alcohol.MassLoss*AlocholEvaLatentHeat;



	DTemperatureDt /= (Water.Mass*WaterSpHeat + Alcohol.Mass*AlocholSpHeat+(NaCl.Mass+Fluorescein.Mass)*NaClSpHeat);

	//no Specific heat data for Fluorescein




	P_T(p) += DTemperatureDt*P_DT(p);

#if !RP_HOST 
#endif


}


