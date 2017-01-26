#ifndef GABOR_H
#define GABOR_H
#include <stdio.h>
//#include <direct.h>
#include "cxcore.h"
#include <iostream>

#include "constants.h"
#include "cv.h"
#include "highgui.h"
#include "cvaux.h"



//#define M_PI 3.14159265359 // already defined in math.h
#define IMAGE_SIZE 128
#define DSCALE 8

extern bool is_gabor_lbp;

class GABOR
{
public:

  GABOR(const char* filename);
 

  ~GABOR();
  int gabor_filter(char* base_fld);//declare gabor filter
 
 int gabor_filtered_data(const char ch,IplImage* img, IplImage* avg_img  );
  
  
private:

  char* gabor_f;
private:
  int save_mat_image(CvMat* mat, char* name);
  double MeanVector(double* v, int vSize);
  void ZeroMeanUnitLength( double* v, int vSize);
  CvMat*  GetMat2( FILE *fh, bool isVector );
  CvMat*  GetMat( FILE *fh, bool isVector );

  int WriteMat(CvMat* m, FILE *fh, bool isVector ) ;
  int WriteMat2(CvMat* m, FILE *fh, bool isVector ) ;
  CvMat** LoadGaborFFT(char* fldname);
  void UnloadGaborFFT(CvMat** mGabor);
  int Mulfft3( const CvArr* srcAarr, const CvArr* srcBarr, CvArr* dstarr );

  int writeData(double* v, int length, FILE* fh);
  double* getData(FILE* fh, int &length);

  int gabor_extraction(IplImage* img,double* object,CvMat** mGabor);

  double* extract_features(const IplImage* img,int& nsize);
  
  int gabor_kernel(int scale,int orientation,int mask_size,double kmax,double sigma,char* filename);

  void calc_average_img(IplImage* avg_img);
  
public:
  void printFilename(){
	  std::cout<<"Gabor Filename: " << std::endl;
	  for (size_t i=0; i<6; i++){
		  std::cout<<gabor_f[i]<<std::endl;
	  }
  }
  

};


  #endif
