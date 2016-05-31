#include <stdio.h>
#include <math.h>
#include "LPG.h"

#define TEMP_LOWER 22715
#define TEMP_UPPER 36615

double table_59e(double rho_x, double TF)
{
	rho_x = (int)(rho_x*10.0+0.5)/10.0;
	if (TF > 0)
		TF = (int)(TF*20.0+0.5)/20.0;
	else 
		TF = (int)(TF*20.0-0.5)/20.0;
//	printf ("T59/1: rho_x = %f, TF = %f.\n", rho_x, TF);
	
	double Tx = TF+273.15;
//	printf ("T59/2: Tx = %f\n", Tx);

	double gamma_x = rho_x/999.016;
//	printf ("T59/3: gamma_x = %f.\n", gamma_x);

	if ((int)(Tx*100.0+0.5) < TEMP_LOWER || (int)(Tx*100.0+0.5) > TEMP_UPPER)
	{
		printf ("T59/4: Tx is out of range.\n");
		return -1;
	}
	if (gamma_x < 0.20995 || gamma_x >= 074005)
	{
		printf ("T59/4: gamma_x is out of range.\n");
		return -1;
	}
//	printf ("T59/4: Tx and rho_x are within range, continue.\n");
	
	
	double T_BK = Tx*1.8-459.67;
	//LPG-23E中的double gamma_x、double TF、gamma_array[0][3]不需要精确
	double gamma_60 = table_23e(gamma_x, T_BK);
//	printf ("T59/5: gamma_60 = %f.\n", gamma_60);
	
	double C_TL = calc_ctl(gamma_60, 293.15*1.8-459.67);
	if (C_TL == -1)
	{
//		printf ("T59/7: C_TL == -1, quit.\n");
		return -1;
	}
//	printf ("T59/7: C_TL = %f.\n", C_TL);

	if (gamma_60 == -1)
	{
//		printf ("T59/7: gamma_60 == -1, quit.\n");
		return -1;
	}
//	printf ("T59/7: gamma_60 = %f.\n", gamma_60);
	double gamma_20 = gamma_60*C_TL;
	
	double pho_20 = gamma_20*999.016;
//	printf ("T59/8: pho_20 = %f.\n", pho_20);
	pho_20 = (int)(pho_20*10.0+0.5)/10.0;
//	printf ("T59/8: pho_20 = %.1f, has been rounded.\n", pho_20);
	
	return pho_20;
}

int main(void)
{
/*	table_59e(210.00, -44.500);
	table_59e(532.57, -44.120);
	table_59e(673.66, -23.330);
	table_59e(245.49, 87.77);
	table_59e(499.55, 87.820);
	table_59e(395.09, 15.43);
	table_59e(449.59, 93.020);
	table_59e(600.74, 80.650);
	table_59e(736.80, -44.230);
	table_59e(224.56, 30.680);
	table_59e(339.63, 18.130);
	table_59e(727.86, -33.070);
	table_59e(209.74, 11.530);
	table_59e(739.32, 11.530);
	table_59e(645.62, -46.030);
	table_59e(645.62, 93.070);
	
*/	
	FILE * fp;
	//int density_low = 355, density_high = 687;
	int density_low = 210, density_high = 735;
	int length1 = (density_high-density_low)/5+2;
	int temp_low = -400, temp_high = 900;
	//int temp_low = 310, temp_high = 700;
	int length2 = (temp_high-temp_low)/5+2;
	fp = fopen("table_59e.txt","w");
	fprintf(fp, "//该系数为 实际值*100000\n//注意：若取值<0，则无效。\n");
	fprintf(fp, "//密度范围：210～735，单位：kg/m³\n//温度范围：-40.0～+90.0，单位：℃\n");
	
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
				fprintf(fp, "%7.0f, ", (table_59e(density_low+j*5-5, (temp_low+i*5-5)/10.0)*10));
				//fprintf(fp, "%10f", table_59e(density_low+j*5-5, (temp_low+i*5-5)/10.0));
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