#include <stdio.h>
#include <math.h>

    double table1[12][8] = {\
/*	| r60	 |   Tc  |  zc	   | pc   |	 	k1		|		k2	      |	    k3	    |	    k4	    |*/
	0.325022, 298.11, 0.27998, 6.250, 2.54616855327, -0.058244177754, 0.803398090807, -0.745720314137,\
	0.355994, 305.33, 0.28220, 6.870, 1.89113042610, -0.370305782347, -0.544867288720, 0.337876634952,\
	0.429277, 333.67, 0.28060, 5.615, 2.20970078464, -0.294253708172, -0.405754420098, 0.319443433421,\
	0.470381, 352.46, 0.27930, 5.110, 2.25341981320, -0.266542138024, -0.372756711655, 0.384734185665,\
	0.507025, 369.78, 0.27626, 5.000, 1.96568366933, -0.327662435541, -0.417979702538, 0.303271602831,\
	0.562827, 407.85, 0.28326, 3.860, 2.04748034410, -0.289734363425, -0.330345036434, 0.291757103132,\
	0.584127, 425.16, 0.27536, 3.920, 2.03734743118, -0.299059145695, -0.418883095671, 0.380367738748,\
	0.624285, 460.44, 0.27026, 3.247, 2.06541640707, -0.238366208840, -0.161440492247, 0.258681568613,\
	0.631054, 469.65, 0.27235, 3.200, 2.11263474494, -0.261269413560, -0.291923445075, 0.308344290017,\
	0.657167, 498.05, 0.26706, 2.727, 2.02382197871, -0.423550090067, -1.152810982570, 0.950139001678,\
	0.664064, 507.35, 0.26762, 2.704, 2.17134547773, -0.232997313405, -0.267019794036, 0.378629524102,\
	0.688039, 540.15, 0.26312, 2.315, 2.19773533433, -0.275056764147, -0.447144095029, 0.493770995799};


double check_gamma_range(double num, double lower, double higher)
{
	if (num > higher)
	{
		printf ("gamma_x_60 = %f, is greater than %f. No answer exists.\n", num, higher);
		return -1;
	}
	if (num < lower)
	{
		printf ("gamma_x_60 = %f, is less than %f.\n", num, lower);
		return -1;
	}
	return num;
}

double calc_ctl(double gamma60, double TF)
{
	double gamma60_1, gamma60_2;
	double Tc, Tc_1, Tc_2;
	double Pc_1, Pc_2;
	double delta, Tr_60, h2, tau, tau_x, C_TL, Tr_x, X;
	double Zc_1, Zc_2;
	double sat_60_1, sat_60_2, sat_x_1, sat_x_2;
	double k1_1, k1_2, k2_1, k2_2, k3_1, k3_2, k4_1, k4_2;
	
	for (int i = 11; i >= 0; i--)
	{
		if (table1[i][0] >= gamma60)//greater than(>=) gamma_x 
		{
			gamma60_2 = table1[i][0];
			Tc_2 = table1[i][1];
			Zc_2 = table1[i][2];
			Pc_2 = table1[i][3];
			k1_2 = table1[i][4];
			k2_2 = table1[i][5];
			k3_2 = table1[i][6];
			k4_2 = table1[i][7];
		}
		else						//less than (<) gamma_x
		{
			gamma60_1 = table1[i][0];
			Tc_1 = table1[i][1];
			Zc_1 = table1[i][2];
			Pc_1 = table1[i][3];
			k1_1 = table1[i][4];
			k2_1 = table1[i][5];
			k3_1 = table1[i][6];
			k4_1 = table1[i][7];
			break;
		}		
	}
	
	delta = (gamma60 - gamma60_1)/(gamma60_2 - gamma60_1);
    Tc = Tc_1 + delta*(Tc_2 - Tc_1);
//	printf("gamma60 = %f, gamma60_1 = %f, gamma60_2 = %f\n", gamma60, gamma60_1, gamma60_2);
//	printf("Tc = %f, Tc_1 = %f, Tc_2 = %f, delta = %f.\n", Tc, Tc_1, Tc_2, delta);

	Tr_x = (TF+459.67)/1.8/Tc;
    if (Tr_x > 1.0)
    {
	printf("Tr_x = %f, is greater than 1.0, no solution.\n", Tr_x);
	return -1;
    }
	
	Tr_60 = 519.67/1.8/Tc;
	h2 = Zc_1*Pc_1/Zc_2/Pc_2;
	tau = 1 - Tr_60;
	sat_60_1 = Pc_1*(1+((k1_1*pow(tau, 0.35)+k3_1*pow(tau, 2)+k4_1*pow(tau, 3))/(1+k2_1*pow(tau, 0.65))));
	sat_60_2 = Pc_2*(1+((k1_2*pow(tau, 0.35)+k3_2*pow(tau, 2)+k4_2*pow(tau, 3))/(1+k2_2*pow(tau, 0.65))));
	
	X = sat_60_1/(1+delta*(sat_60_1/h2/sat_60_2-1));
	
	tau_x = 1-Tr_x;
	sat_x_1 = Pc_1*(1+((k1_1*pow(tau_x, 0.35)+k3_1*pow(tau_x, 2)+k4_1*pow(tau_x, 3))/(1+k2_1*pow(tau_x, 0.65))));
	sat_x_2 = Pc_2*(1+((k1_2*pow(tau_x, 0.35)+k3_2*pow(tau_x, 2)+k4_2*pow(tau_x, 3))/(1+k2_2*pow(tau_x, 0.65))));
	
	C_TL = sat_x_1/X/(1+delta*(sat_x_1/h2/sat_x_2-1));
	
	return C_TL;
}

double calc_gamma60(double gamma_array[2][5], double table1[12][8], double TF)
{
	double delta;
	
	//**************T23/7**************
	if (gamma_array[0][0] != -1)	//gamma60_low exists.
	{
		delta = (gamma_array[1][3]-gamma_array[1][0])/(gamma_array[1][1]-gamma_array[1][0]);
		if (delta < 0.001)
			delta = 0.001;
		if (delta > 0.999)
			delta = 0.999;
		gamma_array[0][2] = gamma_array[0][0]+delta*(gamma_array[0][1]-gamma_array[0][0]);		
	}
	if (gamma_array[1][0] == -1)	//gamma_x_low does not exist.
		gamma_array[0][2] = (gamma_array[0][0]+gamma_array[0][1])/2;
	double C_TL = calc_ctl(gamma_array[0][2], TF);	
	gamma_array[1][2] = C_TL*gamma_array[0][2];
	
	//**************T23/8**************
	
	if ((gamma_array[1][3] >= gamma_array[1][1] && gamma_array[1][3] <= gamma_array[1][2] && fabs(gamma_array[0][0]-gamma_array[0][2]) <= 0.00000001)||\
		(gamma_array[1][3] >= gamma_array[1][2] && gamma_array[1][3] <= gamma_array[1][1] && fabs(gamma_array[0][1]-gamma_array[0][2]) <= 0.00000001))
	{
		gamma_array[0][3] = gamma_array[0][2];
		//	T54E使用的时候，不需要精确
		//gamma_array[0][3] = (int)(gamma_array[0][3]*10000.0+0.5)/10000.0;
	
		return 1;
	}

	//**************T23/9**************
	double alpha = gamma_array[0][1]-gamma_array[0][0];
	double beta = gamma_array[1][1]*gamma_array[1][1]-gamma_array[1][0]*gamma_array[1][0];
	double phi = (gamma_array[1][1]-gamma_array[1][0])/(gamma_array[1][2]-gamma_array[1][0]);
	double A = (alpha-phi*(gamma_array[0][2]-(gamma_array[0][0])))\
				/(beta-phi*(gamma_array[1][2]*gamma_array[1][2]-gamma_array[1][0]*gamma_array[1][0]));
	double B = (alpha-A*beta)/(gamma_array[1][1]-gamma_array[1][0]);
	double C = gamma_array[0][0]-B*gamma_array[1][0]-A*pow(gamma_array[1][0], 2);
	
	gamma_array[0][4] = A*pow(gamma_array[1][3], 2)+B*gamma_array[1][3]+C;
	
	if (gamma_array[0][4] < gamma_array[0][0])
		gamma_array[0][4] = gamma_array[0][0]+\
		(gamma_array[0][2]-gamma_array[0][0])*(gamma_array[1][3]-gamma_array[1][0])\
		/(gamma_array[1][2]-gamma_array[1][0]);
	if (gamma_array[0][4] > gamma_array[0][1])
		gamma_array[0][4] = gamma_array[0][2]+\
		(gamma_array[0][1]-gamma_array[0][2])*(gamma_array[1][3]-gamma_array[1][2])\
		/(gamma_array[1][1]-gamma_array[1][2]);
	
	//**************T23/10**************
	C_TL = calc_ctl(gamma_array[0][4], TF);
	gamma_array[1][4] = C_TL*gamma_array[0][4];
	
	if (fabs(gamma_array[1][4]-gamma_array[1][3]) <= 0.00000001)
	{
		gamma_array[0][3] = gamma_array[0][4];
		//	T54E使用的时候，不需要精确
		//gamma_array[0][3] = (int)(gamma_array[0][3]*10000.0+0.5)/10000.0;
		return 1;
	}
	
	//**************T23/11**************

	if (gamma_array[1][4] > gamma_array[1][3])
	{
		gamma_array[1][1] = gamma_array[1][4];
		gamma_array[0][1] = gamma_array[0][4];
	}
	if (gamma_array[1][2] < gamma_array[1][3])
	{
		gamma_array[1][0] = gamma_array[1][2];
		gamma_array[0][0] = gamma_array[0][2];
	}
	if (gamma_array[1][4] < gamma_array[1][3])
	{
		gamma_array[1][0] = gamma_array[1][4];
		gamma_array[0][0] = gamma_array[0][4];
	}
	if (gamma_array[1][2] > gamma_array[1][3])
	{
		gamma_array[1][1] = gamma_array[1][2];
		gamma_array[0][1] = gamma_array[0][2];
	}	
}


double table_23e(double gamma_x, double TF)
{
	double Tx, Tr_x, Tr_60, tau, tau_x, sat_60, sat_x, gamma_x_ref;
	int Tx_upper = 36615, Tx_lower = 22715, time = 0;
	double density_upper = 0.7400, density_lower = 0.2100;
	double gamma60_1 = 0, gamma60_2 = 0;
	double Tr_x_1, Tr_x_2;
	double gamma60 = 0, gamma_x_1 = 0, gamma_x_2 = 0;
	double Tc_2, Tc_1;
	/*
	gamma60: 	gamma60_low, gamma60_high, gamma60_mid, gamma60, gamma60_trial
	gamma_x:	gamma_x_low, gamma_x_high, gamma_x_mid, gamma_x, gamma_x_trial
	*/
	double gamma_array[2][5] = {0};
    
/*		Tr_60		sat_60		gamma_x_ref		Tr_x	sat_x*/	
	double table2[12][5] = {0};
	
	/*	T54E使用的时候，不需要精确
	if (TF < 0)
	{
		TF = (int)(TF*10.0-0.5)/10.0;
	}
	else
	{
		TF = (int)(TF*10.0+0.5)/10.0;
	}
	gamma_x = (int)(gamma_x*10000.0+0.5)/10000.0;*/
	gamma_array[1][3] = gamma_x;
	 
	Tx = (TF+459.67)/1.8;
	if	((int)(Tx*100.0+0.5) > Tx_upper )
    {
		printf("Tx = %f, is greater than 366.15, no solution\n", Tx);
		return -1;
	}
	if	((int)(Tx*100.0+0.5) < Tx_lower )
    {
		printf("Tx = %f, is less than 227.15, no solution.\n", Tx);
		return -1;
    }
	
    if (gamma_x > density_upper)
    {
		printf("gamma_x = %f, is greater than %f. no solution.\n", gamma_x, density_upper);
		return -1;
    }
	if (gamma_x < density_lower)
    {
		printf("gamma_x = %f, is less than %f. no solution.\n", gamma_x, density_lower);
		return -1;
    }
		
	for (int i = 0; i < 12; i++)
	{
		Tr_x = table2[i][3] = Tx/table1[i][1];
		Tr_60 = table2[i][0] = 519.67/1.8/table1[i][1];
		
		
		if (Tr_x <= 1.0)
		{	
			tau = 1 - Tr_60;
			sat_60 = table2[i][1]  = table1[i][3]*(1+((table1[i][4]*pow(tau, 0.35)\
					+table1[i][6]*pow(tau, 2)+table1[i][7]*pow(tau, 3))/(1+table1[i][5]*pow(tau, 0.65))));
			tau_x = 1-Tr_x;
			sat_x = table2[i][4] = table1[i][3]*(1+((table1[i][4]*pow(tau_x, 0.35)\
					+table1[i][6]*pow(tau_x, 2)+table1[i][7]*pow(tau_x, 3))/(1+table1[i][5]*pow(tau_x, 0.65))));
			gamma_x_ref = table2[i][2] = table1[i][0]*sat_x/sat_60;			
		}
		else
		{
			printf ("Tr_x = %f, is larger than 1.\n", Tr_x);
			printf ("The reference fluid will not be a liquid at this observed temperature.\n");
			gamma_x_ref = table2[i][2] = -1*table1[i][0];
		}
	}
	
	for (int i = 11; i >= 0; i--)
	{
		if (table2[i][2] >= gamma_x)//greater than(>=) gamma_x 
		{
			gamma_x_2 = table2[i][2];
			gamma60_2 = table1[i][0];
			Tr_x_2 = table2[i][3];
			Tc_2 = table1[i][1];
		}
		else						//less than (<) gamma_x
		{
			gamma_x_1 = table2[i][2];
			gamma60_1 = table1[i][0];
			Tr_x_1 = table2[i][3];
			Tc_1 = table1[i][1];
			break;
		}		
	}
	if ((int)(gamma_x_1*10000) == 0)//gamma_x is lower than "EE 68/32"
	{	
		gamma_x_1 = table2[0][2];	//"EE (68/32)" 	for fluid '1'
		gamma_x_2 = table2[1][2];	//"ethane" 		for fluid '2'
		gamma60_1 = table1[0][0];
		gamma60_2 = table1[1][0];
	}
	if ((int)(gamma_x_2*10000) == 0)//gamma_x is greater than "n-heptane"
	{	
		gamma_x_2 = gamma_x_1;		//"n-heptane" for fluid '1' and fluid '2'
		gamma60_2 = table1[11][0];
	}
	//**************T23/6**************
	if (gamma_x > gamma_x_2)
	{
		printf ("gamma_x = %f, is greater than gamma_x_2(%f), no answer exists.\n", gamma_x, gamma_x_2);
		return -1;
	}
	gamma_array[0][0] = gamma60_1;
	gamma_array[0][1] = gamma60_2;
	gamma_array[1][1] = gamma_x_2;
	gamma_array[1][0] = gamma_x_1;
	
	if (gamma_x < gamma_x_1)
	{
		printf ("gamma_x = %f, is less than gamma_x_1(%f), no answer exists.\n", gamma_x, gamma_x_1);
		return -1;
	}
	
	if (Tr_x_1 > 1.0)
	{
		gamma_array[0][0] = (Tx-Tc_1)/(Tc_2-Tc_1)*(gamma60_2-gamma60_1)+gamma60_1;
		double C_TL = calc_ctl(gamma_array[0][0], TF);						//计算C_TL的值		
		gamma_array[1][0] = C_TL*gamma_array[0][0];
	}
	else
		gamma_array[0][0] = gamma60_1;
	if (gamma_array[0][0] < 0.3500)
	{
		gamma_array[0][0] = 0.3500;
		double C_TL = calc_ctl(gamma_array[0][0], TF);						//计算C_TL的值
		
		gamma_array[1][0] = C_TL*gamma_array[0][0];
		if (gamma_x < gamma_array[1][0])
		{
			printf ("gamma_x = %f, is lower than gamma_x_low(%f), no answer exists.\n", gamma_x, gamma_array[1][0]);
			return -1;
		}
	}
	for (time = 0; time < 10; time++)
	{
		if (calc_gamma60(gamma_array, table1, TF) == 1)
			break;
	}
	if (time == 10)
	{
		printf ("10 iterations are reached, no solution\n");
		return -1;
	}
	else
		return gamma_array[0][3];
}

/*
int main(void)
{
	//table_23e(0.67432, -23.33);
	//table_23e(0.24573, 189.98);
	//table_23e(0.50004, 190.04);
	//table_23e(0.22238, 87.25);
	//table_23e(0.34006, 64.63);
	//table_23e(0.72858, -27.53);
	//table_23e(0.45572, -24.67);
	//table_23e(0.25776, 179.28);
	//table_23e(0.39548, 59.78);
	//table_23e(0.21056, 87.46);
	//table_23e(0.45003, 199.43);
	//table_23e(0.60133, 177.17);
	//table_23e(0.73592, -44.13);
	//table_23e(0.21, 189.4);
	//table_23e(0.20993, 187.94);
	//table_23e(0.74005, -28.48);
	//table_23e(0.74, -50);
	//table_23e(0.5, -50.8);
	table_23e(0.5, -50.9);
	
	return 0;
}*/