#include "CPistonGap.h"
#include "..\..\caspar_input\input.h"
#include "logger.h"
#pragma once


extern class CGapInput myGapInput;
extern struct sGapResult myGapResult;
extern class CMesh myMesh;
extern class CThermal myThermal;
extern class CFEMThermal myFEMThermal;
extern class CGapUtils myGapUtils;
extern class input myinput;


//Calculate istantaneous thermal flux generated from gap to piston and cylinder
void CPistonGap::PistonCylinderCalcGapThermalFlux(void)
{	
	double lambda;
	Range all = Range::all();


	//Heat fluxes clear
	QgapK = 0.0;
	QgapB = 0.0;

	//Oil conduction coefficient
	lambda = my_oil ->get_lambda();//oilpistongap.oillambda;

	//CONDUCTIVE HEAT FLUX TO SOLID PARTS FROM FLUID [W/m2]
	Array<double,1> To(N*M);
	//Cylinder
	To = T(Range(0,N*M-1));
	QgapB = lambda/dz2 * (To - TB_surf_gap);
	//Piston
	To = T(Range(N*M*(Q-1),N*M*Q-1));
	QgapK = lambda/dz2 * (To - TK_surf_gap);


};

//Calculate istantaneous thermal flux generated from gap to piston and cylinder
void CPistonGap::PistonCylinderCalcGapThermalFlux_new(void)
{	
	double lambda;
	Range all = Range::all();


	//Heat fluxes clear
	QgapK = 0.0;
	QgapB = 0.0;

	//Oil conduction coefficient
	lambda = my_oil ->get_lambda();//oilpistongap.oillambda;

	//CONDUCTIVE HEAT FLUX TO SOLID PARTS FROM FLUID [W/m2]
	Array<double,1> To(N*M);
	//Cylinder
	To = T(Range(0,N*M-1));
	QgapB = lambda/(dz2/2) * (To - TB_surf_gap) - lambda/dz2/Q * (TK_surf_gap - TB_surf_gap);
	//Piston
	To = T(Range(N*M*(Q-1),N*M*Q-1));
	QgapK = lambda/(dz2/2) * (To - TK_surf_gap) - lambda/dz2/Q * (TB_surf_gap - TK_surf_gap);


};

//Calculate thermal flux to solid surfaces considering boundary flux over shaft revolution
void CPistonGap::PistonCylinderCalcBodyThermalFlux(void)
{
	double phi_deg,phi_rad,speedK,timeold,TDC,TCase,
		AlphaCase,AlphaDC,dt,lvar,lA,lch,lB;
	int phioffset;
	Range all = Range::all();
	Array<double,1> Ebody;
	Array<double,1> Egap;
	Array<double,1> qDC;
	Array<double,1> qCase;
	

	//Parameters
	speedK = operatingpistongap.speedK;
	phi_rad = operatingpistongap.phi_rad;
	phi_deg = operatingpistongap.phi_deg;
	phioffset = (int) (floor(phi_rad/dphi)*speedK);
	timeold = myGapResult.time;


	//Geometry
	lvar = geometrypistongap.lvar;
	lA = geometrypistongap.lA;
	lB = geometrypistongap.lB;
	lch = geometrypistongap.lch;


	//Temperatures
	if(phi_deg < 180.0)
	{
		TDC = temperaturepistongap.THP;
	}
	else
	{ 
		TDC = temperaturepistongap.TLP;
	}
	TCase = temperaturepistongap.TCase;	
	//Covection Case side
	AlphaDC = oilpistongap.AlphaDC;								
	//Covection DC side
	AlphaCase = oilpistongap.AlphaCase;
	//Delta Time
	dt = timenew-timeold;


	//----------------------PISTON INTERPOLATION-------------------------//
	Egap.resize(QgapK.extent(0));					Egap = 0.0;
	Ebody.resize(myGapInput.xyzfK_th.extent(0));	Ebody = 0.0;
	qDC.resize(myGapInput.xyzfK_th.extent(0));		qDC = 0.0;
	qCase.resize(myGapInput.xyzfK_th.extent(0));	qCase = 0.0;
	//Gap istantaneous energy [J/m2] = [W/m2] * [s]
	Egap = QgapK * dt;
	//Energy flux rotation due to piston relative motion in gap
	if(phioffset>0)
	{
		Array<double,1> EgapCopy;
		EgapCopy.resize(N*M);
		EgapCopy = Egap;
		//Egap(Range(phioffset*M,N*M-1)) = EgapCopy(Range(0,(N-phioffset)*M-1));
		//Egap(Range(0,phioffset*M-1)) = EgapCopy(Range((N-phioffset)*M,N*M-1));
		Egap(Range(0,(N-phioffset)*M-1)) = EgapCopy(Range(phioffset*M,N*M-1));
		Egap(Range((N-phioffset)*M,N*M-1)) = EgapCopy(Range(0,phioffset*M-1));
	};
	//Energy flux to solid body interpolation
	Ebody = myGapUtils.InterpolateFields(FaceIdK_f2s_th,FaceDistK_f2s_th,Egap);
	//Boundary heat fluxes [W/m2]
	qDC = AlphaDC * (TDC - TK_surf);
	qCase = AlphaCase * (TCase - TK_surf);
	//Boundary energy flux [J/m2]
	qDC *= dt;	qCase *= dt;
	//Boundaries
	//Ebody = where( zfK_th<(lch+lA), qDC , Ebody );
	//Ebody = where( zfK_th>(lch+lA+lvar), qCase , Ebody );
	Ebody = where( zfK_th<(lch+lA), 0 , Ebody );
	Ebody = where( zfK_th>(lch+lA+lvar), 0 , Ebody );
	tbodyK_DC += where( zfK_th<(lch+lA), dt , 0 );
	tbodyK_case += where( zfK_th>(lch+lA+lvar), dt , 0 );
	//Final
	EbodyK += Ebody;


	//----------------------CYLINDER INTERPOLATION-------------------------//
	Egap.resize(QgapB.extent(0));					Egap = 0.0;
	Ebody.resize(myGapInput.xyzfB_th.extent(0));	Ebody = 0.0;
	qDC.resize(myGapInput.xyzfB_th.extent(0));		qDC = 0.0;
	qCase.resize(myGapInput.xyzfB_th.extent(0));	qCase = 0.0;
	//Gap istantaneous energy [J/m2] = [W/m2] * [s]
	Egap = QgapB * dt;
	//Energy flux to solid body interpolation
	Ebody = myGapUtils.InterpolateFields(FaceIdB_f2s_th,FaceDistB_f2s_th,Egap);
	//Boundary heat fluxes [W/m2]
	qDC = AlphaDC * (TDC - TB_surf);
	qCase = AlphaCase * (TCase - TB_surf);
	//Boundary energy flux [J/m2]
	qDC *= dt;	qCase *= dt;
	//Boundaries
	//Ebody = where( zfB_th<lB, qDC , Ebody );
	//Ebody = where( zfB_th>(lB+lvar), qCase , Ebody );
	Ebody = where( zfB_th<lB, 0 , Ebody );
	Ebody = where( zfB_th>(lB+lvar), 0 , Ebody );
	tbodyB_DC += where( zfB_th<lB, dt , 0 );
	tbodyB_case += where( zfB_th>(lB+lvar), dt , 0 );
	//Final
	EbodyB += Ebody;

}

void CPistonGap::PistonCylinderConstructWeightedMatrix(void)
{
	Array<double,1> hbody;
	Array<double,1> hgap;
	double phi_deg,phi_rad,speedK,timeold,phioffset,dt,lvar,lA,lB,lch;


	//Parameters
	speedK = operatingpistongap.speedK;
	phi_rad = operatingpistongap.phi_rad;
	phi_deg = operatingpistongap.phi_deg;
	phioffset = (int) (floor(phi_rad/dphi)*speedK);
	timeold = myGapResult.time;

	//Geometry
	lvar = geometrypistongap.lvar;
	lA = geometrypistongap.lA;
	lB = geometrypistongap.lB;
	lch = geometrypistongap.lch;

	//Delta Time
	dt = timenew-timeold;

	//--------------------------Cylinder to Piston----------------------------
	hgap.resize(h.extent(0));
	hgap = h;
	if(phioffset>0)
	{
		Array<double,1> hgapCopy;
		hgapCopy.resize(N*M);
		hgapCopy = hgap;
		//hgap(Range(phioffset*M,N*M-1)) = hgapCopy(Range(0,(N-phioffset)*M-1));
		//hgap(Range(0,phioffset*M-1)) = hgapCopy(Range((N-phioffset)*M,N*M-1));
		hgap(Range(0,(N-phioffset)*M-1)) = hgapCopy(Range(phioffset*M,N*M-1));
		hgap(Range((N-phioffset)*M,N*M-1)) = hgapCopy(Range(0,phioffset*M-1));
	};
	//Gap height to solid body interpolation
	hbody.resize(myGapInput.xyzfK_th.extent(0));	hbody = 10.0;
	hbody = myGapUtils.InterpolateFields(FaceIdK_f2s_th,FaceDistK_f2s_th,hgap);

	for(int i_face=0; i_face<hbody.size(); i_face++)
	{
		if(zfK_th(i_face)>(lch+lA) && zfK_th(i_face) < (lch+lA+lvar))
		{
			double denominator = 0;
				for(int i_nb=0; i_nb<FaceId_B2K_th.extent(1); i_nb++)
				denominator += 1 / FaceDist_B2K_th(i_face,i_nb);

			for(int i_nb=0; i_nb<FaceId_B2K_th.extent(1); i_nb++)
			{
				double weight = (1/FaceDist_B2K_th(i_face,i_nb)) / denominator;
				int FLAG = 0;
				for(int i=0; i<Mid_B2K[i_face].size(); i++)
				{
					if(FaceId_B2K_th(i_face,i_nb) == Mid_B2K[i_face][i])
					{
						FLAG = 1;
						Mwt_B2K[i_face][i] += 1 / hbody(i_face) * dt * weight;
						break;
					}
				}
				if (FLAG == 0)//new row
				{
					Mid_B2K[i_face].push_back(FaceId_B2K_th(i_face,i_nb));
					Mwt_B2K[i_face].push_back(1 / hbody(i_face) * dt * weight);
				}
			}
		}
		//cout<<"Mwt_B2K["<<i_face<<"].size(): "<<Mwt_B2K[i_face].size()<<endl;
	}

	//--------------------------Piston to Cylinder----------------------------
	hgap.resize(h.extent(0));
	hgap = h;
	//Gap height to solid body interpolation
	hbody.resize(myGapInput.xyzfB_th.extent(0));	hbody = 10.0;
	hbody = myGapUtils.InterpolateFields(FaceIdB_f2s_th,FaceDistB_f2s_th,hgap);

	for(int i_face=0; i_face<hbody.size(); i_face++)
	{
		if(zfB_th(i_face)>(lB) && zfB_th(i_face) < (lB+lvar))
		{
			double denominator = 0;
			for(int i_nb=0; i_nb<FaceId_K2B_th.extent(1); i_nb++)
				denominator += 1 / FaceDist_K2B_th(i_face,i_nb);
	
			for(int i_nb=0; i_nb<FaceId_K2B_th.extent(1); i_nb++)
			{
				double weight = (1/FaceDist_K2B_th(i_face,i_nb)) / denominator;
				int FLAG = 0;
				for(int i=0; i<Mid_K2B[i_face].size(); i++)
				{
					if(FaceId_K2B_th(i_face,i_nb) == Mid_K2B[i_face][i])
					{
						FLAG = 1;
						Mwt_K2B[i_face][i] += 1 / hbody(i_face) * dt * weight;
						break;
					}
				}
				if (FLAG == 0)
				{
					Mid_K2B[i_face].push_back(FaceId_K2B_th(i_face,i_nb));
					Mwt_K2B[i_face].push_back(1 / hbody(i_face) * dt * weight);
				}
			}
		}
		//cout<<"Mwt_K2B["<<i_face<<"].size(): "<<Mwt_K2B[i_face].size()<<endl;
	}

}

//Calculate thermal flux to solid surfaces considering boundary flux over shaft revolution
void CPistonGap::PistonCylinderCalcBodyThermalFlux_new(void)
{
	double phi_deg,phi_rad,speedK,timeold,TDC,TCase,
		AlphaCase,AlphaDC,dt,lvar,lA,lch,lB;
	int phioffset;
	Range all = Range::all();
	Array<double,1> Ebody;
	Array<double,1> Egap;
	Array<double,1> qDC;
	Array<double,1> qCase;
	

	//Parameters
	speedK = operatingpistongap.speedK;
	phi_rad = operatingpistongap.phi_rad;
	phi_deg = operatingpistongap.phi_deg;
	phioffset = (int) (floor(phi_rad/dphi)*speedK);
	timeold = myGapResult.time;


	//Geometry
	lvar = geometrypistongap.lvar;
	lA = geometrypistongap.lA;
	lB = geometrypistongap.lB;
	lch = geometrypistongap.lch;


	//Temperatures
	if(phi_deg < 180.0)
	{
		TDC = temperaturepistongap.THP;
	}
	else
	{ 
		TDC = temperaturepistongap.TLP;
	}
	TCase = temperaturepistongap.TCase;	
	//Covection Case side
	AlphaDC = oilpistongap.AlphaDC;								
	//Covection DC side
	AlphaCase = oilpistongap.AlphaCase;
	//Delta Time
	dt = timenew-timeold;


	//----------------------PISTON INTERPOLATION-------------------------//
	Egap.resize(QgapK.extent(0));					Egap = 0.0;
	Ebody.resize(myGapInput.xyzfK_th.extent(0));	Ebody = 0.0;
	qDC.resize(myGapInput.xyzfK_th.extent(0));		qDC = 0.0;
	qCase.resize(myGapInput.xyzfK_th.extent(0));	qCase = 0.0;
	//Gap istantaneous energy [J/m2] = [W/m2] * [s]
	Egap = QgapK * dt;
	//Energy flux rotation due to piston relative motion in gap
	if(phioffset>0)
	{
		Array<double,1> EgapCopy;
		EgapCopy.resize(N*M);
		EgapCopy = Egap;
		Egap(Range(phioffset*M,N*M-1)) = EgapCopy(Range(0,(N-phioffset)*M-1));
		Egap(Range(0,phioffset*M-1)) = EgapCopy(Range((N-phioffset)*M,N*M-1));
	};
	//Energy flux to solid body interpolation
	Ebody = myGapUtils.InterpolateFields(FaceIdK_f2s_th,FaceDistK_f2s_th,Egap);
	//Boundary heat fluxes [W/m2]
	qDC = AlphaDC * (TDC - TK_surf);
	qCase = AlphaCase * (TCase - TK_surf);
	//Boundary energy flux [J/m2]
	qDC *= dt;	qCase *= dt;
	//Boundaries
	Ebody = where( zfK_th<(lch+lA), qDC , Ebody );
	Ebody = where( zfK_th>(lch+lA+lvar), qCase , Ebody );
	//Final
	EbodyK += Ebody;


	//----------------------CYLINDER INTERPOLATION-------------------------//
	Egap.resize(QgapB.extent(0));					Egap = 0.0;
	Ebody.resize(myGapInput.xyzfB_th.extent(0));	Ebody = 0.0;
	qDC.resize(myGapInput.xyzfB_th.extent(0));		qDC = 0.0;
	qCase.resize(myGapInput.xyzfB_th.extent(0));	qCase = 0.0;
	//Gap istantaneous energy [J/m2] = [W/m2] * [s]
	Egap = QgapB * dt;
	//Energy flux to solid body interpolation
	Ebody = myGapUtils.InterpolateFields(FaceIdB_f2s_th,FaceDistB_f2s_th,Egap);
	//Boundary heat fluxes [W/m2]
	qDC = AlphaDC * (TDC - TB_surf);
	qCase = AlphaCase * (TCase - TB_surf);
	//Boundary energy flux [J/m2]
	qDC *= dt;	qCase *= dt;
	//Boundaries
	Ebody = where( zfB_th<lB, qDC , Ebody );
	Ebody = where( zfB_th>(lB+lvar), qCase , Ebody );
	//Final
	EbodyB += Ebody;

}

//Calculate piston and cylinder body temperatures and thermal deformation
void CPistonGap::PistonCylinderSolveBodyThermal(void)
{

	/*Log << "EbodyK = " << "\n";
	for( int i = 0; i < EbodyK.size();i++){
		Log << EbodyK(i) << "\n";
	}

		Log << "EbodyB = " << "\n";
	for( int i = 0; i < EbodyB.size();i++){
		Log << EbodyB(i) << "\n";
	}*/

	//Dump IM's dwm
	/*myGapInput.xyzfK_p.free(); myGapInput.xyznK_p.free(); myGapInput.xyzfB_p.free();
	myGapInput.xyznB_p.free();*/ myGapInput.IM_piston.free(); myGapInput.IM_cylinder.free();

	//shaft speed [rev/s]
	double speed = operatingpistongap.speed/60.0;

	//---------------Calculate temperature piston------------//
	qbi_piston[0] = 0.0;
	//[J/(m2 rev) * rev/s] = [W/m2] 
	EbodyK *= speed;
	//limit flux
	//EbodyK = where(EbodyK>=1.0e5,1.0e5,EbodyK);
	EbodyK = where(EbodyK<=-1.0e3,-1.0e3,EbodyK);
	EbodyK = where(EbodyK>=5.0e4,5.0e4,EbodyK);
	//EbodyK = where(EbodyK<=-5.0e4,-5.0e4,EbodyK);
	//damp flux
	EbodyK = EbodyK_old + AlphaTh * (EbodyK - EbodyK_old);
	//assign flux to flux vector - from faces heat flux to surface nodes
	qbi_piston[0] = EbodyK;
	/*Log << "Boundary EbodyK " << "\n";
	for(int b = 0; b < EbodyK.size(); b++){
		Log << EbodyK(b) << "\n";
	}*/
	for(int i = 0;i<myinput.data.thermal.piston.neumann_bc.size();i++){
		Array<double,1> temp ;
		temp.resize(myMesh.nCells,myinput.data.thermal.piston.neumann_bc[i].q[0]);
		temp[0]= myinput.data.thermal.piston.neumann_bc[i].q[0];
		qbi_piston.push_back(temp);
	};
	qbi_piston[0] = where(qbi_piston[0]==0.0,1.0e-9,qbi_piston[0]);
	//read mesh and solve thermal
	myThermal.readMeshThermal("piston");
	for(int a = 0; a < qbi_piston.size(); a++){
		/*Log << "Boundary qbi_piston " << a << "\n";
		for(int b = 0; b < qbi_piston[a].size(); b++){
			Log << qbi_piston[a](b) << "\n";
		}*/
	}
	myThermal.ThermalSolve("piston",qbi_piston,TK_body,TK_surf);
	//assign flux to old array
	EbodyK_old = EbodyK;
	//reset flux
	EbodyK = 0.0;

	//---------------Solve piston thermal deformation------------//
	if(ThermalDeformation)
	{
		myFEMThermal.FEMThermalSolve("piston",TK_body,defK_th,myMesh.piston_name);
	}
	//free variables
	myMesh.xyz.free(); myMesh.conn.free(); myMesh.cxyz_th.free(); myMesh.E.free();
	myMesh.v.free(); myMesh.alpha.free(); myThermal.nodeid_qb.clear(); myThermal.phi.clear();



	//---------------Calculate temperature cylinder------------//
	qbi_cylinder[0] = 0.0;
	//[J/(m2 rev) * rev/s] = [W/m2] 
	EbodyB *= speed;
	//limit flux
	//EbodyB = where(EbodyB>=1.0e5,1.0e5,EbodyB);
	EbodyB = where(EbodyB<=-1.0e4,-1.0e4,EbodyB);
	EbodyB = where(EbodyB>=5.0e4,5.0e4,EbodyB);
	//EbodyB = where(EbodyB<=-5.0e4,-5.0e4,EbodyB);
	//damp flux
	EbodyB = EbodyB_old + AlphaTh * (EbodyB - EbodyB_old);
	//assign flux to flux vector
	qbi_cylinder[0] = EbodyB;
	/*Log << "Boundary EbodyB " << "\n";
	for(int b = 0; b < EbodyB.size(); b++){
		Log << EbodyB(b) << "\n";
	}*/
	qbi_cylinder[0] = where(qbi_cylinder[0]==0.0,1.0e-9,qbi_cylinder[0]);
	//read mesh
	myThermal.readMeshThermal("cylinder");
	//calculate other bores heat flux
	if(myinput.data.options_piston.general.EHDTestRig || myinput.data.options_piston.general.TriboTestRig)
		qbi_cylinder.resize(1);
	else
		qbi_cylinder.resize(myinput.data.operating_conditions.npistons);
	PistonCylinderBlockFlux(qbi_cylinder);
	//add other defined flux boundaries
	for(int i = 0;i<myinput.data.thermal.block.neumann_bc.size();i++){
		Array<double,1> temp; 
		temp.resize(myMesh.nNodes);//should be more than plenty...
		temp = myinput.data.thermal.block.neumann_bc[i].q[0];
		//Log << "defining qb " << i << "\n";
		qbi_cylinder.push_back(temp);
	};
	//solve thermal
	for(int a = 0; a < qbi_cylinder.size(); a++){
		//Log << "Boundary qbi_cylinder " << a << "\n";
		for(int b = 0; b < qbi_cylinder[a].size(); b++){
			//Log << qbi_cylinder[a](b) << "\n";
		}
	}
	myThermal.ThermalSolve("cylinder",qbi_cylinder,TB_body,TB_surf);
	//assign flux to old array
	EbodyB_old = EbodyB;
	//reset flux
	EbodyB = 0.0;
	
	//---------------Solve cylinder thermal deformation------------//
	if(ThermalDeformation)
	{
		myFEMThermal.FEMThermalSolve("cylinder",TB_body,defB_th,myMesh.cylinder_name);
	}
	//free variables
	myMesh.xyz.free(); myMesh.conn.free(); myMesh.cxyz_th.free(); myMesh.E.free();
	myMesh.v.free(); myMesh.alpha.free(); myThermal.nodeid_qb.clear(); myThermal.phi.clear(); 

	//Load IM's dwm
	if(myinput.data.options_piston.general.PressureDeformation)
	{
		//read surface coordinates
		//Log << "Reading pressure mesh coordinates file... " << "\t";
		//myGapInput.readBodySurfacexyzPressure();
		//Log << "done!" << "\n";
		//Log << " " << "\n";
		//read matrices
		Log << "Reading influence matrices... " << "\t";
		myGapInput.readInfluenceMatricesPistonCylinder(myinput.data.options_piston.general.IM_piston_path,myinput.data.options_piston.general.IM_bushing_path);
		Log << "done!" << "\n";
		Log << " " << "\n";
	};

};

//Calculate piston and cylinder body temperatures and thermal deformation
void CPistonGap::PistonCylinderSolveBodyThermal_new(void)
{

	myGapInput.IM_piston.free(); myGapInput.IM_cylinder.free();

	//shaft speed [rev/s]
	double speed = operatingpistongap.speed/60.0;
	double T_DC = 0.5 * (temperaturepistongap.THP + temperaturepistongap.TLP);
	double T_case = temperaturepistongap.TCase;
	vector<double> Tb;
	Tb.push_back(T_DC);
	Tb.push_back(T_case);

	//---------------Calculate temperature piston------------//
	qbi_piston[0] = 0.0;
	//[J/(m2 rev) * rev/s] = [W/m2] 
	EbodyK *= speed;
	tbodyK_DC *= speed*oilpistongap.AlphaDC;
	tbodyK_case *= speed*oilpistongap.AlphaCase;
	//limit flux
	//EbodyK = where(EbodyK>=1.0e5,1.0e5,EbodyK);
	//EbodyK = where(EbodyK<=-1.0e3,-1.0e3,EbodyK);
	//EbodyK = where(EbodyK>=5.0e4,5.0e4,EbodyK);
	//EbodyK = where(EbodyK<=-5.0e4,-5.0e4,EbodyK);
	//damp flux
	EbodyK = EbodyK_old + AlphaTh * (EbodyK - EbodyK_old);
	//assign flux to flux vector - from faces heat flux to surface nodes
	qbi_piston[0] = EbodyK;
	vector<Array<double,1>> wtd_h_K;
	wtd_h_K.push_back(tbodyK_DC);
	wtd_h_K.push_back(tbodyK_case);
	/*Log << "Boundary EbodyK " << "\n";
	for(int b = 0; b < EbodyK.size(); b++){
		Log << EbodyK(b) << "\n";
	}*/
	for(int i = 0;i<myinput.data.thermal.piston.neumann_bc.size();i++){
		Array<double,1> temp ;
		temp.resize(myMesh.nCells,myinput.data.thermal.piston.neumann_bc[i].q[0]);
		temp[0]= myinput.data.thermal.piston.neumann_bc[i].q[0];
		qbi_piston.push_back(temp);
	};
	qbi_piston[0] = where(qbi_piston[0]==0.0,1.0e-9,qbi_piston[0]);
	//read mesh and solve thermal
	myThermal.readMeshThermal("piston");
	myThermal.ThermalSolve_Piston(qbi_piston,wtd_h_K,Tb,TK_body,TK_surf);
	//assign flux to old array
	EbodyK_old = EbodyK;
	//reset flux
	EbodyK = 0.0;
	tbodyK_DC = 0.0;
	tbodyK_case = 0.0;
	//free variables
	myMesh.xyz.free(); myMesh.conn.free(); myMesh.cxyz_th.free(); myMesh.E.free();
	myMesh.v.free(); myMesh.alpha.free(); myThermal.nodeid_qb.clear(); myThermal.phi.clear();
	wtd_h_K.clear();
	cout<<"Finish piston-piston matrix"<<endl;
	//---------------Calculate temperature cylinder------------//
	qbi_cylinder[0] = 0.0;
	//[J/(m2 rev) * rev/s] = [W/m2] 
	EbodyB *= speed;
	tbodyB_DC *= speed*oilpistongap.AlphaDC;
	tbodyB_case *= speed*oilpistongap.AlphaCase;;
	//limit flux
	//EbodyB = where(EbodyB>=1.0e5,1.0e5,EbodyB);
	//EbodyB = where(EbodyB<=-1.0e4,-1.0e4,EbodyB);
	//EbodyB = where(EbodyB>=5.0e4,5.0e4,EbodyB);
	//EbodyB = where(EbodyB<=-5.0e4,-5.0e4,EbodyB);
	//damp flux
	EbodyB = EbodyB_old + AlphaTh * (EbodyB - EbodyB_old);
	//assign flux to flux vector
	qbi_cylinder[0] = EbodyB;
	vector<Array<double,1>> wtd_h_B;
	wtd_h_B.push_back(tbodyB_DC);
	wtd_h_B.push_back(tbodyB_case);
	/*Log << "Boundary EbodyB " << "\n";
	for(int b = 0; b < EbodyB.size(); b++){
		Log << EbodyB(b) << "\n";
	}*/
	qbi_cylinder[0] = where(qbi_cylinder[0]==0.0,1.0e-9,qbi_cylinder[0]);
	//read mesh
	myThermal.readMeshThermal("cylinder");
	//calculate other bores heat flux
	if(myinput.data.options_piston.general.EHDTestRig || myinput.data.options_piston.general.TriboTestRig)
		qbi_cylinder.resize(1);
	else
		qbi_cylinder.resize(myinput.data.operating_conditions.npistons);
	PistonCylinderBlockFlux(qbi_cylinder);
	//add other defined flux boundaries
	for(int i = 0;i<myinput.data.thermal.block.neumann_bc.size();i++){
		Array<double,1> temp; 
		temp.resize(myMesh.nNodes);//should be more than plenty...
		temp = myinput.data.thermal.block.neumann_bc[i].q[0];
		//Log << "defining qb " << i << "\n";
		qbi_cylinder.push_back(temp);
	};
	//solve thermal
	for(int a = 0; a < qbi_cylinder.size(); a++){
		//Log << "Boundary qbi_cylinder " << a << "\n";
		for(int b = 0; b < qbi_cylinder[a].size(); b++){
			//Log << qbi_cylinder[a](b) << "\n";
		}
	}
	myThermal.ThermalSolve_Cylinder(qbi_cylinder,wtd_h_B,Tb,TB_body,TB_surf);
	//assign flux to old array
	EbodyB_old = EbodyB;
	//reset flux
	EbodyB = 0.0;
	tbodyB_DC = 0.0;
	tbodyB_case = 0.0;
	//free variables
	myMesh.xyz.free(); myMesh.conn.free(); myMesh.cxyz_th.free(); myMesh.E.free();
	myMesh.v.free(); myMesh.alpha.free(); myThermal.nodeid_qb.clear(); myThermal.phi.clear();
	wtd_h_B.clear();
	cout<<"Finish cylinder-cylinder matrix"<<endl;

	//---------------Solve interbody conduction------------//
	double lambda_oil = my_oil -> get_lambda();//oilpistongap.oillambda;
	double C = operatingpistongap.speed * lambda_oil / 60 / 9;
	myThermal.ThermalSolve_interbody_conduction(Mid_K2B, Mwt_K2B, Mid_B2K, Mwt_B2K, C);
	cout<<"Finish piston-cylinder and cylinder piston matrix"<<endl;

	//---------------Solve both body temperature distribution------------//
	myThermal.ThermalSolve_new();
	cout<<"Finish both body temperature distribution"<<endl;

	//read mesh 
	myThermal.readMeshThermal("piston");
	//read surface coordinates piston body
	myGapInput.readBodySurfacexyzThermal("piston",myinput.data.thermal.piston.meshfile);
	myThermal.ThermalSolve_post("piston",qbi_piston,TK_body,TK_surf);
	if(ThermalDeformation)
	{
		myFEMThermal.FEMThermalSolve("piston",TK_body,defK_th,myMesh.piston_name);
	};
	myMesh.xyz.free(); myMesh.conn.free(); myMesh.E.free(); myMesh.v.free();
	myMesh.alpha.free(); myThermal.nodeid_qb.clear(); myThermal.phi.clear(); 

	//read mesh 
	myThermal.readMeshThermal("cylinder");
	//read surface coordinates cylinder body
	myGapInput.readBodySurfacexyzThermal("cylinder",myinput.data.thermal.block.meshfile);
	myThermal.ThermalSolve_post("cylinder",qbi_cylinder,TB_body,TB_surf);
	if(ThermalDeformation)
	{
		myFEMThermal.FEMThermalSolve("cylinder",TB_body,defB_th,myMesh.cylinder_name);
	}
	myMesh.xyz.free(); myMesh.conn.free(); myMesh.E.free(); myMesh.v.free();
	myMesh.alpha.free(); myThermal.nodeid_qb.clear(); myThermal.phi.clear(); 












	/*//---------------Solve piston thermal deformation------------//
	if(ThermalDeformation)
	{
		myFEMThermal.FEMThermalSolve("piston",TK_body,defK_th,myMesh.piston_name);
	}
	//free variables
	myMesh.xyz.free(); myMesh.conn.free(); myMesh.cxyz_th.free(); myMesh.E.free();
	myMesh.v.free(); myMesh.alpha.free(); myThermal.nodeid_qb.clear(); myThermal.phi.clear();



	
	
	//---------------Solve cylinder thermal deformation------------//
	if(ThermalDeformation)
	{
		myFEMThermal.FEMThermalSolve("cylinder",TB_body,defB_th,myMesh.cylinder_name);
	}
	//free variables
	myMesh.xyz.free(); myMesh.conn.free(); myMesh.cxyz_th.free(); myMesh.E.free();
	myMesh.v.free(); myMesh.alpha.free(); myThermal.nodeid_qb.clear(); myThermal.phi.clear(); */

	//Load IM's dwm
	if(myinput.data.options_piston.general.PressureDeformation)
	{
		//read surface coordinates
		//Log << "Reading pressure mesh coordinates file... " << "\t";
		//myGapInput.readBodySurfacexyzPressure();
		//Log << "done!" << "\n";
		//Log << " " << "\n";
		//read matrices
		Log << "Reading influence matrices... " << "\t";
		myGapInput.readInfluenceMatricesPistonCylinder(myinput.data.options_piston.general.IM_piston_path,myinput.data.options_piston.general.IM_bushing_path);
		Log << "done!" << "\n";
		Log << " " << "\n";
	};

};

//Interpolate bore flux to other cylinder bores
void CPistonGap::PistonCylinderBlockFlux(vector<Array<double,1>> &qbi)
{

	Range all = Range::all();

	//Assign cylinder gap coordinates for reference bore
	int nFaces_gap = (int) myThermal.faceid_qb[0].size();
	Array<double,2> xyz_0;	xyz_0.resize(nFaces_gap,3);	xyz_0 = 0.0;
	for(int i=0;i<nFaces_gap;i++)
	{
		xyz_0(i,0) = myThermal.xyzf(myThermal.faceid_qb[0][i],0);
		xyz_0(i,1) = myThermal.xyzf(myThermal.faceid_qb[0][i],1);
		xyz_0(i,2) = myThermal.xyzf(myThermal.faceid_qb[0][i],2);
	}
	//Array of coordinates rotation of reference bore to other bores
	Array<double,2> xyz_r;	xyz_r.resize(nFaces_gap,3);	xyz_r = 0.0;
	//Array of coordinates nodes actual i-th bore
	Array<double,2> xyz_i;
	//Angle between cylinder bores
	double dtheta = 2.0*PI/(double) qbi.size();
	double theta = -1.0 * dtheta;
	//Interpolate reference heat flux to other bores coordinates
	int nqb = (int) qbi.size();
	//Log << "nqb=" << nqb << "\n";
	for(int i=0;i<nqb-1;i++)
	{
		//Rotate reference bore according to actual i-th bore position
		xyz_r(all,0) = xyz_0(all,0)*cos(theta) - 1.0 * xyz_0(all,1)*sin(theta);
		xyz_r(all,1) = xyz_0(all,0)*sin(theta) + xyz_0(all,1)*cos(theta);
		xyz_r(all,2) = xyz_0(all,2);

		//Increment theta
		theta -= dtheta;
		//Assign actual i-th cylinder gap coordinates
		nFaces_gap = (int) myThermal.faceid_qb[i+1].size();
		xyz_i.resize(nFaces_gap,3);	xyz_i = 0.0;
		for(int j=0;j<nFaces_gap;j++)
		{
			xyz_i(j,0) = myThermal.xyzf(myThermal.faceid_qb[i+1][j],0);
			xyz_i(j,1) = myThermal.xyzf(myThermal.faceid_qb[i+1][j],1);
			xyz_i(j,2) = myThermal.xyzf(myThermal.faceid_qb[i+1][j],2);
		}
		
		//Search neighbours nodes from reference coordinates to actual i-th cylinder
		Array<double,2> FaceDist; Array<int,2> FaceId;
		myGapUtils.SearchNeighbours(xyz_r,xyz_i,FaceId,FaceDist,10);
		//Interpolate reference heat flux to i-th heat flux
		qbi[i+1] = myGapUtils.InterpolateFields(FaceId,FaceDist,qbi[0]);			
	}
};

void CPistonGap::Heatfluxrecovering(int option)
{
	string heatflux_file;
	string line;
	vector <double> readvalue;
	if (option == 1)
		heatflux_file = "./piston_flux.txt";
	else
		heatflux_file = "./block_flux.txt";
	ifstream fheatflux (heatflux_file);
	//cout<<"check1"<<endl;
	if (fheatflux.is_open())
	{
		//cout<<"check2"<<endl;
		while(fheatflux)
		{
			//cout<<"check3"<<endl;
			getline(fheatflux,line);
			istringstream iss(line);
			double temp;
			while (iss >> temp)
				readvalue.push_back(temp);
			iss.clear();
		}
		fheatflux.clear();
		fheatflux.close();

		if (option == 1)
		{
			for (int i=0;i<EbodyK.size();i++)
			{
				EbodyK(i) = readvalue[i];
				//cout<<readvalue[i]<<endl;
			}
			int i_revrow = 0;
			for (int i=0;i<qbi_piston.size();i++)
			{
				for (int j=0;j<qbi_piston[i].extent(0);j++)
				{
					qbi_piston[i](j) = readvalue[i_revrow];
					i_revrow++;
				}
			}
			EbodyK_old=EbodyK;
		}
		else
		{
			for (int i=0;i<EbodyB.size();i++)
			{
				EbodyB(i) = readvalue[i];
				//cout<<readvalue[i]<<endl;
			}
			int i_revrow = 0;
			for (int i=0;i<qbi_cylinder.size();i++)
			{
				for (int j=0;j<qbi_cylinder[i].extent(0);j++)
				{
					qbi_cylinder[i](j) = readvalue[i_revrow];
					i_revrow++;
				}
			}
			EbodyB_old=EbodyB;
		}
	}
	readvalue.clear();
}