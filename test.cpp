//a simple example of the usage of the LBP class

#include <cv.h>
#include <highgui.h>

#include <cmath>
#include <cstdio>

#include "./source/LBP.h"

#define MAX_LENGTH 256

using namespace std;

int main(int argc, char** argv)
{
	int i, j, k;
	int width, height, tlength;
	
	tlength = 29; //how many images
	char** filename = new char*[tlength];
	for(i=0; i<tlength; i++)
	{
		filename[i] = new char[MAX_LENGTH];
		sprintf(filename[i], "%d.png", i);
	}
	
	uchar*** pppArray = NULL;
	IplImage *pImg = NULL;	
	IplImage *pGray = NULL;
	
	pppArray = new uchar**[tlength];
	if(!pppArray)
	{
		printf("Allocating memory failed.\n");
		return -1;
	}
	
	pImg = cvLoadImage(filename[0]);
	if(!pImg)
	{
		printf("Can not load image.\n");
		return -2;
	}	
	width = pImg->width;
	height = pImg->height;
	printf("The size of the input image is: %d x %d\n\n", width, height);
	
	pGray = cvCreateImage(cvGetSize(pImg), 8, 1);
	
	//load image(s) and write the data into a 3d array
	for(i=0; i<tlength; i++)
	{	
		//load image	
		pImg = cvLoadImage(filename[i]);
		if(!pImg)
		{
			printf("Cannot load image.\n");
			return -2;
		}
		
		if(i>=1)
		{
			if((pImg->width!=width) || (pImg->height!=height))
			{
				printf("All the images in the sequence must be of the same size!\n");
				return -3;
			}
		}
		
		//convert color image to gray image
		cvCvtColor(pImg, pGray, CV_RGB2GRAY);
			
		pppArray[i] = new uchar*[height];
		if(!pppArray[i])
		{
			printf("Allocating memory failed.\n");
			return -1;
		}		
			
		for(j=0; j<height; j++)
		{
			pppArray[i][j] = new uchar[width];
			if(!pppArray[i][j])
			{
				printf("Allocating memory failed.\n");
				return -1;
			}
			
			for(k=0; k<width; k++)
			{
				//write data into the 3d array
				pppArray[i][j][k] = (uchar)CV_IMAGE_ELEM(pGray, char, j, k);  
			}
		}		
	}
	
	
	//calculate the histograms
	LBP lbp;	
	//==========================================================================//
	
	lbp.width = width;         lbp.height = height;    	    lbp.tlength = 1;        	
	lbp.R.xR = 1;              lbp.R.yR = 1;	
	lbp.SN.xy = 8;     	
	lbp.uni = 1;               lbp.interp = 1;              lbp.norm = 1;
	lbp.CreateHistogram(pppArray, 0);        //lbp based histogram
	printf("======== LBP ========\n");
	printf("The histogram is:\n");
	for(i=0; i<lbp.uni_bin.xy+1; i++)
	{
		printf("%f ", lbp.uni_hist.pHist_xy[i]);
	}
	printf("\n");
	printf("The bin number of uniform patterns is: %d\n", lbp.uni_bin.xy);	
	printf("\n\n");

	
	//======================================================================//	
	/*	
	lbp.width = width;     lbp.height = height;    	lbp.tlength = tlength;        
	lbp.R.xR = 1;          lbp.R.yR = 1;            lbp.R.tR = 2;		
	lbp.SN.xy = 8;         lbp.SN.xt = 8;           lbp.SN.yt = 8;
	lbp.uni = 1;           lbp.interp = 1;          lbp.norm = 1;			
	lbp.CreateHistogram(pppArray, 1);        //lbp-top based histogram		
	printf("======== LBP-TOP ========\n");
	printf("The histogram of X-Y plane:\n");
	for(i=0; i<lbp.uni_bin.xy+1; i++)
	{
		printf("%f ", lbp.uni_hist.pHist_xy[i]);
	}
	printf("\n");
	printf("The histogram of X-T plane:\n");
	for(i=0; i<lbp.uni_bin.xt+1; i++)
	{
		printf("%f ", lbp.uni_hist.pHist_xt[i]);
	}
	printf("\n");
	printf("The histogram of Y-T plane:\n");
	for(i=0; i<lbp.uni_bin.yt+1; i++)
	{
		printf("%f ", lbp.uni_hist.pHist_yt[i]);
	}
	printf("\n");	
	printf("The bin numbers of uniform patterns is: %d %d %d\n", lbp.uni_bin.xy, lbp.uni_bin.xt, lbp.uni_bin.yt);	
	printf("\n\n");

	*/
	//==================================================================//
	/*
	lbp.width = width;     lbp.height = height;    	lbp.tlength = tlength;        
	lbp.R.xR = 1;          		                    lbp.R.tR = 2;
	lbp.SN.xy = 4;         
	lbp.uni = 0;           lbp.interp = 1;          lbp.norm = 1;		
	lbp.CreateHistogram(pppArray, 2, 1);   //rilbp based histogram		
	printf("======== RIVLBP ========\n");
	int cnt=0;
	for(i=0; i<lbp.ri_bin; i++)
	{	
		if(lbp.pHist_ri[i])
		{
			cnt++;
		}
		
		printf("%f ", lbp.pHist_ri[i]);
	}		
	printf("\n");
	printf("The bin number is: %d \n", lbp.ri_bin);	
	printf("The number of non-zero rivlbp patterns is: %d ", cnt);
	printf("\n\n");	
	//=============================================================//
	lbp.SN.xy = 4;
	lbp.uni = 1;
	lbp.CreateHistogram(pppArray, 2, 1);   //uniform rilbp based histogram		
	printf("======== uniform RIVLBP ========\n");		
	cnt=0;
	for(i=0; i<lbp.ri_bin; i++)
	{	
		if(lbp.pHist_ri[i])
		{
			cnt++;
		}
		
		printf("%f ", lbp.pHist_ri[i]);
	}		
	printf("\n");
	printf("The bin number of uniform rivlbp is: %d \n", lbp.ri_bin);	
	printf("The number of non-zero uniform rivlbp patterns is: %d ", cnt);
	printf("\n\n");
	*/
	//================================================================//
	//release memory
	if(pppArray)
	{
		for(i=0; i<tlength; i++)
		{
			if(pppArray[i])
			{
				for(j=0; j<height; j++)
				{
					if(pppArray[i][j])
					{
						delete []pppArray[i][j];
					}
				}
				
				delete []pppArray[i];
			}
		}
		
		delete []pppArray;					
	}
	
	cvReleaseImage(&pImg);
	cvReleaseImage(&pGray);

	return 0;
}
