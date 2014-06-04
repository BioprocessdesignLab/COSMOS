
void formula(complex<double> *x,complex<double> *fx,complex<double> **v,complex<double> **sv,int nd1,int ni,int cp,int cm)
{
	int i,j;
	complex<double> s;
 //*****************************************************************************
	
	double
		THA = 0,
		LKR = 0,
		AdoMet = 20,
		Cys = 15,
		Phosphate = 10000,
		Val = 100,
		AK1_kforward_app_exp = 5.65,
		AK1_kreverse_app_exp = 1.57,
		AK1_Lys_Ki_app_exp = 550,
		AK1_AdoMet_Ka_app_exp = 3.5,
		AK1_h_exp = 2,
		AK2_kforward_app_exp = 3.15,
		AK2_kreverse_app_exp = 0.88,
		AK2_Lys_Ki_app_exp = 22,
		AK2_h_exp = 1.1,
		AKI_kforward_app_exp = 0.36,
		AKI_kreverse_app_exp = 0.10,
		AKI_Thr_Ki_app_exp = 124,
		AKI_h_exp = 2,
		AKII_kforward_app_exp = 1.35,
		AKII_kreverse_app_exp = 0.38,
		AKII_Thr_Ki_app_exp = 109,
		AKII_h_exp = 2,
		ASADH_kforward_app_exp = 0.9,
		ASADH_kreverse_app_exp = 0.23,
		DHDPS1_k_app_exp = 1,
		DHDPS1_Lys_Ki_app_exp = 10,
		DHDPS1_h_exp = 2,
		DHDPS2_k_app_exp = 1,
		DHDPS2_Lys_Ki_app_exp = 33,
		DHDPS2_h_exp = 2,
		HSDHI_kforward_app_exp = 0.84,
		HSDHI_Thr_relative_residual_activity_app_exp = 0.15,
		HSDHI_Thr_relative_inhibition_app_exp = 0.85,
		HSDHI_Thr_Ki_app_exp = 400,
		HSDHII_kforward_app_exp = 0.64,
		HSDHII_Thr_relative_residual_activity_app_exp = 0.25,
		HSDHII_Thr_relative_inhibition_app_exp = 0.75,
		HSDHII_Thr_Ki_app_exp = 8500,
		HSK_kcat_app_exp = 2.8,
		HSK_Hser_app_exp = 14,
		TS1_kcatmin_exp = 0.42,
		TS1_AdoMet_kcatmax_exp = 3.5,
		TS1_AdoMEt_Km_no_AdoMet_exp = 250,
		TS1_AdoMet_Ka1_exp = 73,
		TS1_AdoMet_Ka2_exp = 0.5,
		TS1_AdoMet_Ka3_exp = 1.09,
		TS1_AdoMet_Ka4_exp = 142,
		TS1_Phosphate_Ki_exp = 1000,
		TS1_h_exp = 2,
		CGS_kcat_exp = 30,
		CGS_Cys_Km_exp = 460,
		CGS_Phser_Km_exp = 2500,
		CGS_Phosphate_Ki_exp = 2000,
		TD_k_app_exp = 0.0124,
		TD_Ile_Ki_no_Val_app_exp = 30,
		TD_Val_Ka1_app_exp = 73,
		TD_Val_Ka2_app_exp = 615,
		TD_h_app_exp = 3,
		Lys_tRNAS_Lys_Km = 25,
		Thr_tRNAS_Vmax = 0.43,
		Thr_tRNAS_Thr_Km = 100,
		Ile_tRNAS_Vmax = 0.43,
		Ile_tRNAS_Ile_Km = 20,
		THA_kcat_exp = 1.7,
		THA_Thr_Km_exp = 7100,
		LKR_kcat_exp = 3.1,
		LKR_Lys_Km_exp = 13000
		;

	complex<double>
		Vak1 = x[8]*(AK1_kforward_app_exp - AK1_kreverse_app_exp*x[1])/(1.0+pow((x[3]/(AK1_Lys_Ki_app_exp/(1.0+AdoMet/ AK1_AdoMet_Ka_app_exp))), AK1_h_exp)),
		Vak2 = x[9]*(AK2_kforward_app_exp - AK2_kreverse_app_exp*x[1])/(1.0+pow((x[3]/ AK2_Lys_Ki_app_exp), AK2_h_exp)),
		VakI = x[10]*(AKI_kforward_app_exp - AKI_kreverse_app_exp*x[1])*1.0/(1.0+pow((x[6]/ AKI_Thr_Ki_app_exp), AKI_h_exp)),
		VakII = x[11]*(AKII_kforward_app_exp - AKII_kreverse_app_exp*x[1])/(1.0+pow((x[6]/ AKII_Thr_Ki_app_exp), AKII_h_exp)),
		Vasadh = x[12]*(ASADH_kforward_app_exp*x[1]- ASADH_kreverse_app_exp*x[2]),
		Vdhdps1 = x[13]* DHDPS1_k_app_exp *x[2]*(1.0/(1.0+pow((x[3]/DHDPS1_Lys_Ki_app_exp), DHDPS1_h_exp))),
		Vdhdps2 = x[14]* DHDPS2_k_app_exp *x[2]*(1.0/(1.0+pow((x[3]/ DHDPS2_Lys_Ki_app_exp), DHDPS2_h_exp))),
		VhsdhI = x[10]* HSDHI_kforward_app_exp *x[2]*( HSDHI_Thr_relative_residual_activity_app_exp + HSDHI_Thr_relative_inhibition_app_exp /(1.0+x[6]/HSDHI_Thr_Ki_app_exp)),
		VhsdhII = x[11]* HSDHII_kforward_app_exp *x[2]*( HSDHII_Thr_relative_residual_activity_app_exp+ HSDHII_Thr_relative_inhibition_app_exp /(1.0+x[6]/ HSDHII_Thr_Ki_app_exp)),
		Vhsk = x[16]* HSK_kcat_app_exp *x[4]/( HSK_Hser_app_exp +x[4]),
		Vts1 = x[18]*(TS1_kcatmin_exp + TS1_AdoMet_kcatmax_exp *pow(AdoMet, TS1_h_exp) /TS1_AdoMet_Ka1_exp)/(1.0+pow(AdoMet, TS1_h_exp) / TS1_AdoMet_Ka1_exp)*x[5]/((1.0+Phosphate/TS1_Phosphate_Ki_exp)*(TS1_AdoMEt_Km_no_AdoMet_exp *(1.0+AdoMet/ TS1_AdoMet_Ka2_exp)/(1.0+AdoMet/ TS1_AdoMet_Ka3_exp))/(1.0+pow(AdoMet, TS1_h_exp) / TS1_AdoMet_Ka4_exp)+x[5]),
		Vcgs = x[17]*( CGS_kcat_exp /(1.0+ CGS_Cys_Km_exp /Cys))*x[5]/(( CGS_Phser_Km_exp /(1.0+ CGS_Cys_Km_exp /Cys))*(1.0+Phosphate/ CGS_Phosphate_Ki_exp)+ x[5]),
		Vtd = x[19]*TD_k_app_exp*x[6]/(1.0+pow((x[7]/( TD_Ile_Ki_no_Val_app_exp + TD_Val_Ka1_app_exp *Val/( TD_Val_Ka2_app_exp +Val))), TD_h_app_exp)),
		VLys_tRNAS = x[15] *x[3]/(Lys_tRNAS_Lys_Km+x[3]),
		VThr_tRNAS = x[15] *x[6]/(Thr_tRNAS_Thr_Km+x[6]),
		VIle_tRNAS = x[15] *x[7]/(Ile_tRNAS_Ile_Km+x[7]),
		Vtha = THA*THA_kcat_exp*x[6]/(THA_Thr_Km_exp+x[6]),
		Vlkr = LKR*LKR_kcat_exp*x[3]/(LKR_Lys_Km_exp+x[3])
		;

	sv[1][1]=Vak1;		sv[1][2]=Vak2;		sv[1][3]=VakI;  sv[1][4]=VakII; sv[1][5]=Vasadh;	sv[1][6]=0.0;		sv[1][7]=0.0;	sv[1][8]=0.0;
	sv[2][1]=Vasadh;	sv[2][2]=0.0; 		sv[2][3]=0.0;   sv[2][4]=0.0; 	sv[2][5]=Vdhdps1;	sv[2][6]=Vdhdps2;	sv[2][7]=VhsdhI;sv[2][8]=VhsdhII;
	sv[3][1]=Vdhdps1;  	sv[3][2]=Vdhdps2;  	sv[3][3]=0.0;   sv[3][4]=0.0;	sv[3][5]=VLys_tRNAS;sv[3][6]=Vlkr;		sv[3][7]=0.0;	sv[3][8]=0.0;
	sv[4][1]=VhsdhI;	sv[4][2]=VhsdhII;	sv[4][3]=0.0;	sv[4][4]=0.0;	sv[4][5]=Vhsk;		sv[4][6]=0.0;		sv[4][7]=0.0;	sv[4][8]=0.0;
	sv[5][1]=Vhsk;  	sv[5][2]=0.0;		sv[5][3]=0.0;	sv[5][4]=0.0;	sv[5][5]=Vcgs;		sv[5][6]=Vts1;		sv[5][7]=0.0;	sv[5][8]=0.0;
	sv[6][1]=Vts1;  	sv[6][2]=0.0;		sv[6][3]=0.0;	sv[6][4]=0.0;	sv[6][5]=Vtd;		sv[6][6]=VThr_tRNAS;sv[6][7]=Vtha;	sv[6][8]=0.0;
	sv[7][1]=Vtd;		sv[7][2]=0.0;		sv[7][3]=0.0;	sv[7][4]=0.0;	sv[7][5]=VIle_tRNAS;sv[7][6]=0.0;		sv[7][7]=0.0;	sv[7][8]=0.0;

	//****************************************************************

	for (i=1;i<nd1;i++){
		s=0.0;
		for (j=1;j<cp+1;j++){
			s+=sv[i][j];
		}
		v[i][0]=s;
		s=0.0;
		for (j=cp+1;j<cp+cm+1;j++){
			s+=sv[i][j];
		}
		v[i][1]=s;
		fx[i]=v[i][0]-v[i][1];
	}
}
