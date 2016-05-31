#include <stdio.h>
#include <math.h>
#include "LPG.h"

#define TEMP_LOWER 22715
#define TEMP_UPPER 36615
#define DEN_LOWER 3317
#define	DEN_UPPER 6836

double table_60e(double rho_20, double TF)
{
//	printf ("T54/1: ");
	rho_20 = (int)(rho_20*10.0+0.5)/10.0;
	if (TF > 0)
		TF = (int)(TF*20.0+0.5)/20.0;
	else 
		TF = (int)(TF*20.0-0.5)/20.0;
//	printf ("T54/1: rho_20 = %f, TF = %f.\n", rho_20, TF);
	
	double Tx = TF+273.15;
//	printf ("Tx = %f\n", Tx);
	if ((int)(Tx*100.0+0.5) < TEMP_LOWER || (int)(Tx*100.0+0.5) > TEMP_UPPER)
	{
		printf ("Tx is out of range.\n");
		return -1;
	}
	if ((int)(rho_20*10.0+0.5) < DEN_LOWER || (int)(rho_20*10.0+0.5) > DEN_UPPER)
	{
		printf ("rho_20 is out of range.\n");
		return -1;
	}
//	printf ("T54/2: Tx = %f.\n", Tx);
//	printf ("Tx and rho_20 are within range, continue.\n");
	
	double gamma_TB = rho_20/999.016;
//	printf ("T54/4: gamma_TB = %f.\n", gamma_TB);
	
	double T_BK = 293.15*1.8-459.67;
	//LPG-23E中的double gamma_x、double TF、gamma_array[0][3]不需要精确
	double gamma_60 = table_23e(gamma_TB, T_BK);
//	printf ("T54/5: gamma_60 = %f.\n", gamma_60);
	
	if (gamma_60 < 0.34995 || gamma_60 >= 0.68805)
	{
//		printf ("T54/6: gamma_60 = %f, is out of range[0.34995, 0.68805), the standard don't apply.\n", gamma_60);
		return -1;
	}
//	printf ("T54/6:  amma_60 = %.12f, is within the range[0.34995, 0.68805), continue.\n", gamma_60);
	
	double C_TL1 = calc_ctl(gamma_60, Tx*1.8-459.67);
	if (C_TL1 == -1)
	{
//		printf ("T54/7:  C_TL1 == -1, quit.\n");
		return -1;
	}
//	printf ("T54/7:  C_TL1 = %f.\n", C_TL1);
	
	double C_TL2 = calc_ctl(gamma_60, T_BK);
	if (C_TL2 == -1)
	{
//		printf ("T54/8:  C_TL2 == -1, quit.\n");
		return -1;
	}
//	printf ("T54/8:  C_TL2 = %f.\n", C_TL2);
	
	double C_TL = C_TL1/C_TL2;
//	printf ("T54/9:  C_TL = %f.\n", C_TL);
	
	if (C_TL <= 0)
	{
//		printf ("T54/10: C_TL = %f, is less than or equal 0. Quit.\n", C_TL);
		return -1;
	}
//	printf ("T54/10: C_TL = %f, is within the range. Continue.\n", C_TL);
		
	C_TL = (int)(C_TL*100000.0+0.5)/100000.0;
//	printf ("T54/11: C_TL = %f, has been rounded. \n", C_TL);
	return C_TL;
}

int main(void)
{
	table_60e(332.69, -5.020);
	table_60e(399.55, -3.920);
	table_60e(451.09, 30.774);
	table_60e(489.92, 84.975);
	table_60e(539.49, 68.360);
	table_60e(569.42, -16.090);
	table_60e(599.37, 43.360);
	table_60e(624.42, 76.650);
	table_60e(639.41, -24.460);
	table_60e(659.38, 80.580);
	table_60e(669.38, 82.790);
	table_60e(399.83, 90.570);
	table_60e(449.56, -46.030);
	table_60e(331.59, 64.440);
	table_60e(449.56, 93.030);
	table_60e(683.65, -17.780);
	table_60e(683.64, 93.020);
	table_60e(331.67, -46.020);
	
	
	FILE * fp;
	int density_low = 335, density_high = 680;
	//int density_low = 400, density_high = 435;
	int length1 = (density_high-density_low)/5+2;
	int temp_low = -400, temp_high = 500;
	//int temp_low = 100, temp_high = 295;
	int length2 = (temp_high-temp_low)/5+2;
	fp = fopen("table_60e.txt","w");
	fprintf(fp, "//该系数为 实际值*100000\n//注意：若取值<0，则无效。\n");
	fprintf(fp, "//密度范围：335～680，单位：kg/m³\n//温度范围：-40.0～+50.0，单位：℃\n");
	
	for (int i = 0; i < length2; i++)
	{
		for (int j = 0; j < length1; j++)
		{
			if (i == 0 ) //第一行显示密度值
				{
				if ( j == 0)//第一行 第一列显示空
					fprintf(fp, "%s", "//");
				else
					fprintf(fp, "%-9d", density_low+j*5-5);
				}
			//if (j == 0 && i != 0)	//第一列显示温度值
			//	fprintf(fp, "%-6d", temp_low+i*5-5);
			if (i*j != 0)
			{
				fprintf(fp, "%7.0f, ", (table_60e(density_low+j*5-5, (temp_low+i*5-5)/10.0)*100000));
				//fprintf(fp, "%10f", table_60e(density_low+j*5-5, (temp_low+i*5-5)/10.0));
			}
		}
		if (i)
		{
			fprintf(fp, "\t//%-6.1f℃\n", (float)(temp_low+i*5-5)/10.0);
		}
		else
			fprintf(fp, "\tkg/m³\n");
			//fputc('\n\tkg/m³', fp);
			
	}
	fflush(fp);
	fclose(fp);
	
	return 0;
}