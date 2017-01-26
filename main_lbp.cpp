#include "cv.h"
#include "highgui.h"
#include<iostream>
#include<cstring>
//#include<cmath>
#include "./source/constants.h"
#include "./source/LBP.h"
#include "./source/face_features.h"
//this code uses only basic local binary pattern- no time information is used here.
//Declaration of functions here
//#include "./source/gabor_extract.h"
bool is_gabor_lbp = true; // global variable 	

int main ()
{

  
  // int  width=0, height=0;
  int  frame_no =0; // NR is number of row=width, NC is number of column=height
  // uchar** ppArray = NULL;
  IplImage* img = NULL;
  //IplImage* gray_img = NULL;
  bool video = true; // true if video is used 
  bool lbp_whole_img = true; // true if whole face image is considered
 		     
  //whole face and block wise cannot be true together
  bool is_block_wise = false;// true if each image is divided into sub
  // for lbp i.e only 59 dimensional feature
  // vector
  // blocks 
  face_features f_features;

  // store all the file path in an array of string
  
  size_t no_file =1; // number of file 
  char** filepath = NULL;

  filepath = new char* [no_file];
  //===============================================================//
  //===============================================================//
  
  for (int t= 0; t< no_file; t++)

    {
      filepath[t] = new char [100] ;
      if(!filepath[t]){
	printf("Filepath memory allocation failed ! \n");
	return -1;
      }

    }
  //=======================================================================//

  //===========MMI database videos
  //========================================//
  
  for( int i=0; i<no_file; i++);
  {
    //============happiness data ========//
    /*
    strcpy (filepath[0], "../../videos/emotion/happiness/S001-008.avi");
    strcpy (filepath[1],"../../videos/emotion/happiness/S001-108.avi");
    strcpy (filepath[2],"../../videos/emotion/happiness/S001-109.avi");
    strcpy (filepath[3], "../../videos/emotion/happiness/S002-107.avi");
 //=========sadness data================//	
    strcpy (filepath[4],"../../videos/emotion/sadness/S001-114.avi");
    strcpy (filepath[5], "../../videos/emotion/sadness/S001-115.avi");
    strcpy (filepath[6],"../../videos/emotion/sadness/S002-113.avi");
    strcpy (filepath[7], "../../videos/emotion/sadness/S002-114.avi");
 //===============Disgust data ==========//
    strcpy (filepath[8], "../../videos/emotion/disgust/S001-104.avi");
    strcpy (filepath[9], "../../videos/emotion/disgust/S001-105.avi");
    strcpy (filepath[10], "../../videos/emotion/disgust/S002-103.avi");    
    strcpy (filepath[11], "../../videos/emotion/disgust/S002-104.avi");
//================Anger data =============//
    strcpy (filepath[12], "../../videos/emotion/anger/S001-097.avi");
    strcpy (filepath[13], "../../videos/emotion/anger/S001-100.avi");
    strcpy (filepath[14], "../../videos/emotion/anger/S001-101.avi");
    strcpy (filepath[15], "../../videos/emotion/anger/S002-099.avi");
    strcpy (filepath[16],	"../../videos/emotion/anger/S002-100.avi");
    //========Surprise data ================//
    strcpy (filepath[17], "../../videos/emotion/surprise/brow_ris.avi");
    strcpy (filepath[18], "../../videos/emotion/surprise/S001-116.avi");
    strcpy (filepath[19],"../../videos/emotion/surprise/S001-117.avi");
    strcpy (filepath[20], "../../videos/emotion/surprise/S002-117.avi");
    strcpy (filepath[21], "../../videos/emotion/surprise/S002-118.avi");
  //==========Fear data ========================//
    strcpy (filepath[22], "../../videos/emotion/fear/S001-106.avi"); 
    strcpy (filepath[23], "../../videos/emotion/fear/S001-107.avi"); 
    */
 
     strcpy (filepath[0], "../../videos/emotion/fear/S002-105.avi");
    /*
   
     
   
   
   

    
    
   
  
  
  
	
   
  
   
   	
   
   
  
   
   
    */
    
   }
   
   
  //=====================================================//


  //===============CK+ database videos ==========================================//
  /*
    for( int i=0; i<no_file; i++);
    {

    strcpy (filepath[0], "../../videos/emotion/ck/happiness.avi");
    strcpy (filepath[1], "../../videos/emotion/ck/sadness.avi");
    strcpy (filepath[2], "../../videos/emotion/ck/disgust.avi");	
	
    strcpy (filepath[3], "../../videos/emotion/ck/anger.avi");



	
    strcpy (filepath[4], "../../videos/emotion/ck/surprise.avi");
    strcpy (filepath[5], "../../videos/emotion/ck/fear.avi");
	
    //	strcpy (filepath[1], "../../videos/emotion/surprise/S001-117.avi");
    //	strcpy (filepath[2], "../../videos/emotion/surprise/S002-117.avi");

	
    // char* filepath[2] = {
    //   "../../videos/emotion/surprise/S001-116.avi",
    //     "../../videos/emotion/surprise/S002-118.avi" };


	
    }
  */
  //=================================================================//

      

 
  //=================================================================//
  // to write output data into a file 
  std::ofstream fyd;
  fyd.open("yd_data.txt", std::ios_base::out | std::ios::app);
  if(fyd== NULL)
    {
      std::cout<<"err in opening file to write output data "<<"yd_data.txt" <<std::endl;
      return -2;
    }
  //=================================================================//
   
  
  //============================================//
  // CvCapture *input_video = cvCaptureFromFile("home/work/code/main/facial_emotion_detection/videos/emotion/fear/S001-106.avi");//S034-004.avi");
  if(video)
    {

      for(int i=0; i<no_file; i++)
	{
	  //CvCapture *input_video = cvCaptureFromFile("../../videos/emotion/surprise/S002-118.avi");//S034-004.avi");

	  CvCapture *input_video = cvCaptureFromFile( filepath[i] );//S034-004.avi");
	  if(!input_video)
	    {
	      /* Either the video didn't exist OR it uses a codec OpenCV
	       * doesn't support.
	       */
	      std::cout<<"video file path "<<filepath[i]<<std::endl;
	      fprintf(stderr, "Error: Can't open video.\n");
	      exit(EXIT_FAILURE);
	    }
    
	  //======================================//
	  while(cvGrabFrame(input_video)) 
	    // while(true)
	    {
	      // Load image
	      // img = cvLoadImage("./images/an.jpg", 1 );

	      img =  cvRetrieveFrame(input_video);
   
	      if(!img)
		{
		  printf("Can not load  file. \n");
		  return -2;
		}
	      //call face detection code here ====assign face width and height as NR and NC
	      if(!lbp_whole_img)
		f_features.features_lbp(img, is_block_wise);
	      else
		{
		  f_features.face_lbp(img);
	    
		}
	      fyd<<i+1<<std::endl; // write the output file 
	  
	      frame_no++;
	      std::cout<<"frame no = "<<frame_no<<std::endl;

	    }// end of while loop

	}// end of for loop
    } // end of if video
  //==========================================//
  if(!video)
    {

      img = cvLoadImage("./images/3.jpeg", 1 );

      if(!img)
	{
	  std::cout<<"image file not found ! " <<std::endl;

	  return -2;
	}
 
      /*
      //=======call the Gabor function here ==========//
      IplImage* gabor_img= cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, 3 );
      cvCopy(img, gabor_img);

 
      GABOR gabor;
      char ch= 'e';
      
      gabor.gabor_filtered_data(ch,  gabor_img );


      cvReleaseImage(&gabor_img);
      //=============//

      */
      
     
      if(!lbp_whole_img)
	f_features.features_lbp(img, is_block_wise);
      else {

	// detect face here === then pass the face region for LBP code
	// =====// 
	// IplImage* gray_imag1 = NULL;
	// IplImage* face_img = NULL;
	// gray_imag1 = cvCreateImage (cvGetSize(img), 8, 1);
 
	// //convert color image to gray image
	// cvCvtColor(img, gray_imag1, CV_RGB2GRAY);
	// CvRect* wroi= new CvRect ();
	// f_features.face_detect(gray_imag1, wroi);

	// face_img = cvCreateImage (cvSize(wroi->width, wroi->height), 8, 1);
	// f_features.get_ROI_img(gray_imag1, wroi, face_img); 
	// // wroi->x = 0;
	// // wroi->y = 0;
	// // wroi->width = img->width;
	// // wroi->height = img->height;


	
	// char* filename_a = new char [50];

	// strcpy(filename_a, "whole_img_hist.txt");
       
	// f_features.histogram_img_roi(face_img, wroi, filename_a );

	// delete [] filename_a;
	// delete wroi;

	// cvReleaseImage (&gray_imag1);
	// cvReleaseImage (&face_img);

	f_features.face_lbp(img);
      }

     
    }

   
  //============release all  memories==============//
  fyd.close();
  // delete [] filepath;

      
  //del (filepath, no_file);

  for (int t= 0; t < no_file; t++)

    {
      delete [] filepath[t];

    }

  delete [] filepath;
      
  cvReleaseImage (&img);
 
  // cvReleaseCapture(&input_video); //release the video capture memory
  // giving segmentation fault
  cvDestroyAllWindows();
  return 0;
 
}// End of main 
//========================================================================
