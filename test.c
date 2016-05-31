#include <stdio.h>
void main(void)
{
	int num=9; 
	float* pFloat=(&num); 
	printf("num的值为：%d\n", num); 
	printf("*pFloat的值为：%f\n", *pFloat);
	*pFloat=9.0; 
	printf("num=%d\n",num); 
	printf("*pFloat=%f\n",*pFloat); 
	
	float t = 1.1;
	for (int i = 0; i<10; i++)
	{
		printf ("%d: %f.\n", i,t+=1.1);
	}
}