#include <iostream>
#include <cstring>
#include <cmath>
#include <string.h>
#include "face_features.h"
#include "LBP.h"
#include "GABOR.h"


//bool is_gabor_lbp = true;

face_features::face_features()
{
	
  eye1r0_prev = new CvRect ();
  eye2r1_prev = new CvRect ();
  nose_prev = new CvRect ();

  lip_rect_prev= new CvRect ();
  face_rect_prev= new CvRect();
  lip_position_prev= new CvPoint2D32f [4];
 
}


face_features::~face_features()
{
  delete lip_position_prev;
  delete lip_rect_prev;
  delete face_rect_prev;
  delete eye1r0_prev;
  delete eye2r1_prev;
  delete nose_prev;

}


void face_features::face_detect(IplImage *img, CvRect *const face_rect)
{
 
  CvHaarClassifierCascade *cascade_f; //= new CvHaarClassifierCascade ();
  CvMemStorage		*storage; //= new CvMemStorage ();
  char file1[] = "./xml_files/haarcascade_frontalface_alt.xml";
  /* load the face classifier */
  cascade_f = (CvHaarClassifierCascade*)cvLoad(file1, 0, 0, 0);
  /* setup memory storage, needed by the object detector */
  storage = cvCreateMemStorage(0);
  CvSeq *faces = cvHaarDetectObjects(
 				     img, cascade_f, storage,
 				     1.1, 3, CV_HAAR_DO_CANNY_PRUNING, cvSize( 40, 40 ) );
  if (faces->total == 0)
    {
      std::cout<<"no face is detected in this frame!" <<std::endl; 
      return ;
    }
  /* draw a rectangle */
  CvRect *r= (CvRect*)cvGetSeqElem(faces, 0); // check how to release memory
  //==================//
  face_rect->x = r->x;
  face_rect->y = r->y;
  face_rect->height = r->height;
  face_rect->width = r->width;
  //======================//

  cvClearMemStorage(storage); // clear the memory of face data 
  //======releasing memories==========//
  cvReleaseHaarClassifierCascade( &cascade_f );
 
  cvDestroyAllWindows();
  // delete r; 
  return;

} // only face detection is done here in the above function 

//===========================================//



void face_features::face_sub_regions_detect(IplImage *img, CvRect *const face_rect, CvRect *const eye1_r, CvRect *const eye2_r, CvRect *const nose, CvPoint2D32f *const lip_point,  struct CvRect   * const r2)
{
 
  CvHaarClassifierCascade *cascade_f; //= new CvHaarClassifierCascade ();
  CvHaarClassifierCascade *cascade_e; //= new CvHaarClassifierCascade ();
  CvMemStorage		*storage; //= new CvMemStorage ();

  //====================================//
  int extra_eye_width = 10;
  int extra_eye_height = 0;
  int lip_width_extra = 25;
  int chk_sign;
  CvPoint cnt;
  CvPoint2D32f eye_cnt[2];
  // CvRect* r_ik = new CvRect ();
  // CvRect *r =  new CvRect ();
  //CvRect* r_it = new CvRect ();
  //======================================//
  char file1[] = "./xml_files/haarcascade_frontalface_alt.xml";
  char file2[] = "./xml_files/haarcascade_eye.xml";

  /* load the face classifier */
  cascade_f = (CvHaarClassifierCascade*)cvLoad(file1, 0, 0, 0);

  /* load the eye classifier */
  cascade_e = (CvHaarClassifierCascade*)cvLoad(file2, 0, 0, 0);

  /* setup memory storage, needed by the object detector */
  storage = cvCreateMemStorage(0);


  CvSeq *faces = cvHaarDetectObjects(
 				     img, cascade_f, storage,
 				     1.1, 3, CV_HAAR_DO_CANNY_PRUNING, cvSize( 40, 40 ) );

 
  
  if (faces->total == 0)
    {
      std::cout<<"no face is detected in this frame!" <<std::endl; 
      return ;
    }
  /* draw a rectangle */
 

  CvRect *r= (CvRect*)cvGetSeqElem(faces, 0); // check how to release memory
  //==================//
  face_rect->x = r->x;
  face_rect->y = r->y;
  face_rect->height = r->height;
  face_rect->width = r->width;
  //======================//
  //=================//
  r2->x = r->x;
  r2->y = r->y;
  r2->height = r->height;
  r2->width = r->width;
  //====================//

  // to be used for mouth detection
  // std::cout<<"face rect  data are =  "<<face_rect->x <<"\t" <<face_rect->y <<"\t" <<face_rect->height <<"\t"<<face_rect->width<<std::endl;
 

  cvClearMemStorage(storage); // clear the memory of face data 

  //==========================================//
 
  // Set the Region of Interest: estimate the eyes' position //
 
  cvSetImageROI(img, cvRect(r->x, r->y + (r->height/5.5), r->width, r->height/3.0));
 
  // cvSaveImage("a.jpeg", img);
  
  //=========face detection over =========================//
  float scale_factor=0.0;
  if (face_rect->width<250)
    scale_factor=1.10;
  else
    scale_factor=1.15;


  // detect eyes //
  CvSeq* eyes = cvHaarDetectObjects( 
				    img, cascade_e, storage,
				    scale_factor, 3, 0, cvSize(25, 15));
  // cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;

  // if(!eyes)
  //   {

  //     std::cout<< "no eyes are detected in this frame " <<std::endl;
  //     return ;
  //   }
  
  if(eyes->total>=2)
    {

      //copy the 2nd eye seq roi into r_i1
      //copy the 1st eye seq roi into r_i0
  
      CvRect *r_ik = (CvRect*)cvGetSeqElem(eyes,1);
    
      CvRect *r_it = (CvRect *)cvGetSeqElem(eyes, 0 );
  
   
      if((r_ik->x-r_it->x)<0)
	{
	  //	right_eye_first=true;
	  chk_sign=-1;//i.e right eye is detected first
	  eye1_r->x = r_ik->x;
	  eye1_r->y = r_ik->y;
	  eye1_r->height = r_ik->height;
	  eye1_r->width = r_ik->width;//left eye_roi



	  eye2_r->x = r_it->x;//right eye roi
	  eye2_r->y = r_it->y;//right eye roi
	  eye2_r->height = r_it->height;//right eye roi
	  eye2_r->width = r_it->width;//right eye roi
    
	} 
      else
	{
	  chk_sign=1;//i.e left eye is detected first
	  eye1_r->x = r_it->x;//left eye roi
	  eye1_r->y = r_it->y;
	  eye1_r->height = r_it->height;
	  eye1_r->width = r_it->width;


    
	  eye2_r->x = r_ik->x;// right eye roi
	  eye2_r->y = r_ik->y;
	  eye2_r->height = r_ik->height;
	  eye2_r->width = r_ik->width;
    
	}

   
      eye1_r->x = eye1_r->x - extra_eye_width;
      eye2_r->x = eye2_r->x - extra_eye_width;

      eye1_r->width =  eye1_r->width + 2*extra_eye_width;
      eye2_r->width = eye2_r->width + 2*extra_eye_width;

      eye1_r->y =  eye1_r->y + extra_eye_height;
      eye2_r->y = eye2_r->y + extra_eye_height;


      eye1_r->height = eye1_r->height + extra_eye_height;
      eye2_r->height = eye2_r->height + extra_eye_height;
      //======================================================//
    
 
      for( int i = 0; i < (eyes ? eyes->total : 0); i++ ) {
	//  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
	//  cout<<i+1<<" total number of eye_detected="<<"\t"<<eyes->total<<endl;
   
	if (i==0)
	  r=r_it;
	else
	  r=r_ik;
   
	cnt = cvPoint (r->x + (r->width)/2,r->y + (r->height)/2);//get_centr(r);//get the eye center
    
	eye_cnt[i].x=cnt.x + face_rect->x; //center of each eye with respect to
	//the whole face image 
	eye_cnt[i].y=cnt.y + face_rect->y;

	//	std::cout<<"eye centers are  :"<< eye_cnt[i].x <<"\t"<<eye_cnt[i].y<<std::endl; 
      }
  
      //============================================================//
  
      //estimate the lip position having known the eye centers
      float lip_scale=1.0;
  
      if (r2->width<250)
	lip_scale=1.0;
      else
	lip_scale=1.4;
      //  int sign_chk=eye_cnt[1].x-eye_cnt[0].x;
 
      if ( chk_sign<0)
	{
	  //  cout<< "info:  Right eye is detected first" <<endl;
     
	  lip_point[0].x=eye_cnt[1].x - lip_width_extra;//-10 is extra 
	  lip_point[0].y=eye_cnt[1].y  + (-1)* 1.4* lip_scale *(eye_cnt[1].x-eye_cnt[0].x);
	  lip_point[1].x=eye_cnt[0].x + lip_width_extra;// +10 is extra
	  lip_point[1].y=eye_cnt[0].y + (-1)*1.4*lip_scale *(eye_cnt[1].x-eye_cnt[0].x);
	}

      else

	{
	  //  cout<< "info:  Left eye is detected first" <<endl;
	  lip_point[0].x=eye_cnt[0].x-lip_width_extra;// -10 is extra
	  lip_point[0].y=eye_cnt[0].y + 1.4*lip_scale *(eye_cnt[1].x-eye_cnt[0].x);
	  lip_point[1].x=eye_cnt[1].x+lip_width_extra;// +10 is extra
	  lip_point[1].y=eye_cnt[1].y + 1.4*lip_scale *(eye_cnt[1].x-eye_cnt[0].x);
	 
	}
  
      //lip position estimation done.

      //clear the storage for next object detection
      cvClearMemStorage(storage);
 

      cvResetImageROI(img);
  
      //detect mouth using haar cascade :used in the function detect Mouth
      
      r2->x = lip_point[0].x;
      r2->y = lip_point[0].y;
      r2->width = (lip_point[1].x- lip_point[0].x);
      // anima == changed on 19th sept================//
      r2->height = face_rect->height/3.0;

      // std::cout<<" img height === " <<img->height <<"\t"<<"lip=="<<r2->y + r2->height<<std::endl;
      if((r2->y + r2->height) > (img->height))
	r2->height = r2->height - ((r2->y + r2->height)- img->height) ;
      //=========changes done =============// 

	 
  
      // //nose rectangle w.r.t to whole image is // as the eye detection sequence is not know we need to check which one is left/right eye. 
 
  
      nose->x=eye1_r->x + eye1_r->width/2 + face_rect->x;
      nose->y=eye1_r->y + eye1_r->height/2 + face_rect->y + face_rect->height/5.5;
      nose->width=eye2_r->x - eye1_r->x;
      nose->height=face_rect->height/3;
 
  
      
      eye1r0_prev->x = eye1_r->x;
      eye1r0_prev->y = eye1_r->y;
      eye1r0_prev->height = eye1_r->height;
      eye1r0_prev->width = eye1_r->width;

  
      eye2r1_prev->x = eye2_r->x;
      eye2r1_prev->y = eye2_r->y;
      eye2r1_prev->height = eye2_r->height;
      eye2r1_prev->width = eye2_r->width;
  
      nose_prev->x = nose->x;
      nose_prev->y = nose->y;
      nose_prev->height = nose->height;
      nose_prev->width = nose->width;
  
      
      lip_position_prev->x = lip_point->x;
      lip_position_prev->y = lip_point->y;
     
      lip_rect_prev->x = r2->x;
      lip_rect_prev->y = r2->y;
      lip_rect_prev->height = r2->height;
      lip_rect_prev->width = r2->width;
           
    } // end of if condition 

  
  else

    {cout<<"No eyes are detected in the frame!! previous eye's location is set here  "<<endl;

      eye1_r->x = eye1r0_prev->x;
      eye1_r->y = eye1r0_prev->y;
      eye1_r->height = eye1r0_prev->height;
      eye1_r->width = eye1r0_prev->width;
      
      eye2_r->x = eye2r1_prev->x;
      eye2_r->y = eye2r1_prev->y;
      eye2_r->height = eye2r1_prev->height;
      eye2_r->width = eye2r1_prev->width;
      
      nose->x = nose_prev->x;
      nose->y = nose_prev->y;
      nose->height = nose_prev->height;
      nose->width = nose_prev->width;
      
      
      lip_point->x = lip_position_prev->x;
      lip_point->y = lip_position_prev->y; 

      r2->x = lip_rect_prev->x;
      r2->y = lip_rect_prev->y;
      r2->height = lip_rect_prev->height;
      r2->width = lip_rect_prev->width;
  
    }



  // std::cout<<"mem r2="<<r2<<std::endl;

  // std::cout<<"r2 a= "<<r2->x  <<"\t" <<r2->y <<"\t" <<r2->height <<"\t"<<r2->width<<std::endl;

  
  //   //========clear memories=================//

 
  cvResetImageROI(img);

 
  // std::cout<<"eye rect "<<eye1_r->x<<"\t" <<eye1_r->y<<"\t"<<eye1_r->height<<"\t"<<eye1_r->width <<std::endl;
  //====== get the ROI with respect to face ROI===========//

  //rectangle enclosing the 1st eye. The rectangle is w.r.t the whole image. 
  eye1_r->x = face_rect->x + eye1_r->x;
  eye1_r->y = face_rect->y + face_rect->height/5.5+eye1_r->y;
  eye1_r->width = eye1_r->width;//width of the eye should be same.
  eye1_r->height = eye1_r->height;//height of the eye should remain same.
  //=====================================//

  //rectangle enclosing the 2nd eye w.r.t the original image
  eye2_r->x = face_rect->x + eye2_r->x;
  eye2_r->y = face_rect->y + face_rect->height/5.5 + eye2_r->y;
  eye2_r->width = eye2_r->width;//width of the eye should be same.
  eye2_r->height = eye2_r->height;//height of the eye should remain same.

  //std::cout<<"eye rect after"<<eye1_r->x<<"\t" <<eye1_r->y<<"\t"<<eye1_r->height<<"\t"<<eye1_r->width <<std::endl;
  //========================================================//  
  
  //======releasing memories==========//
  cvReleaseHaarClassifierCascade( &cascade_f );
  cvReleaseHaarClassifierCascade( &cascade_e );
  cvClearMemStorage(storage);
   
  //delete r;
  //  delete r_ik;
  //delete r_it;
   
  //std::cout<<"ok"<<std::endl;

  cvDestroyAllWindows();
  return;
}
//====================face detection is done =================//

void face_features::get_ROI_img(IplImage* img, CvRect* r, IplImage* ROI_img )
{
  /* Set the Region of Interest: estimate the eyes' position */
  cvSetImageROI(img, cvRect(r->x, r->y , r->width, r->height));

  //  std::cout<<"I am in get roi img function "<<std::endl;
  // cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;

  //std::cout<<"r= "<<r->height <<"\t"<<r->width<<std::endl;
  // std::cout<<"roi img height and width are :"<<ROI_img->height <<"\t"<<ROI_img->width<<std::endl;
  cvCopy(img, ROI_img);
  //std::cout<<"saving a1 image"<<std::endl;
  cvSaveImage ("a1.jpeg", ROI_img);
  cvResetImageROI(img);
  
  return;
}

void face_features::histogram_img_roi(IplImage* img, CvRect* const ROI, char* const filename )
{

  CvPoint* img_size= new CvPoint ();
  img_size->x = ROI->width;
  img_size->y = ROI->height;
  uchar** ppArray = NULL;
  int height = 0, width = 0;
  IplImage* imag1 = NULL;
 
  width = img_size->x;
  height = img_size->y;
  imag1 = cvCreateImage (cvGetSize(img), 8, 1);
  //  std::cout<<"w = " <<width<<"\t h= "<<height<<std::endl;
  //=====================================//

  LBP lbp;

  lbp.width = width;         lbp.height = height;    	    lbp.tlength = 1;        	
  lbp.R.xR = 1;              lbp.R.yR = 1;	
  lbp.SN.xy = 8;     	
  lbp.uni = 1;               lbp.interp = 1;              lbp.norm = 1;


  //=====================================//
  /*
    cvNamedWindow("test", CV_WINDOW_AUTOSIZE);
    cvShowImage("test", img);
    cvWaitKey(0);
  */
 
  ppArray= new uchar* [height];
  for (int i=0; i<height; i++){
    
    ppArray[i] = new uchar[width]; // each row is having NC column
    if(!ppArray[i]){
      printf("Memory allocation failed ! \n");
      return ;
    }
  }

  for (size_t m=0; m<height; m++){

    for (size_t n=0; n<width; n++){
      //write data into 2D array
      ppArray[m][n] = (uchar)CV_IMAGE_ELEM(img, char, m, n);


      cvSetReal2D(imag1, m, n, ppArray[m][n] );//threshold image
    }
  
  }

  cvSaveImage("a11.jpg", imag1);
  
  // char* filename =  new char[50]; // file name should not be more than
  // 50 character 
  // strcpy_s(filename, 20, "basic_hist.txt");
  // strcpy(filename, "basic_hist.txt");
  
  lbp.CreateHistogram(ppArray, 0, 0,  filename);        //lbp based histogram

  /*
    printf("======== LBP ========\n");
    printf("The histogram is:\n");
    for(int i=0; i<lbp.uni_bin.xy+1; i++){
    printf("%f ", lbp.uni_hist.pHist_xy[i]);
    }
    printf("\n");
    printf("The bin number of uniform patterns is: %d\n", lbp.uni_bin.xy);	
    printf("\n\n");
  */

  //======release memories =======//
  // delete [] filename;
  del (ppArray, height);
  delete img_size;
  cvReleaseImage (&imag1);
  cvDestroyAllWindows();
  return;
}


void face_features::face_lbp(IplImage* img)
{
  //===this part of the code gets the lbp code for face region of an
  //image 
  IplImage* gabor_avgImg = NULL;
  IplImage* gray_imag1 = NULL;
  IplImage* face_img = NULL;
  gray_imag1 = cvCreateImage (cvGetSize(img), 8, 1);
 
  //convert color image to gray image
  cvCvtColor(img, gray_imag1, CV_RGB2GRAY);
  CvRect* wroi= new CvRect ();
  face_detect(gray_imag1, wroi);

  face_img = cvCreateImage (cvSize(wroi->width, wroi->height), 8, 1);
  get_ROI_img(gray_imag1, wroi, face_img); 
  // wroi->x = 0;
  // wroi->y = 0;
  // wroi->width = img->width;
  // wroi->height = img->height;


  if(is_gabor_lbp == true) // is_gabor_lbp is a global variable
    {
  //=======call the Gabor function here ==========//
  IplImage* gabor_img= cvCreateImage(cvSize(face_img->width, face_img->height), IPL_DEPTH_8U, 1 );
  cvCopy(face_img, gabor_img);

 
   GABOR gabor("gabor");
  // GABOR gabor;
  char ch= 'e';

  gabor_avgImg = cvCreateImage(cvSize(IMAGE_SIZE,IMAGE_SIZE), IPL_DEPTH_8U, 1 );

  // gabor.printFilename();
  gabor.gabor_filtered_data(ch,  gabor_img, gabor_avgImg );

//  gabor.printFilename();
  cvSaveImage("avg.jpg", gabor_avgImg);
  

  cvReleaseImage(&gabor_img);
  std::cout<<"check GABOR: "<<std::endl;
  //=============//
    }
  
  CvRect* gabor_roi = new CvRect ();	
  char* filename_a = new char [50];
  char* filename_b = new char [50];
   
  float fdata =0.0;
  strcpy(filename_a, "whole_img_hist.txt");
  strcpy(filename_b, "face_hist_appended.txt");
 
       
  histogram_img_roi(face_img, wroi, filename_a );

  
  if(is_gabor_lbp == true) // is_gabor_lbp is a global variable 
    {
      char* filename_gabor = new char [50];
      char* filename_gabor_app = new char [50];
      strcpy(filename_gabor, "gabor_avg_hist.txt");
      strcpy(filename_gabor_app, "gabor_avg_hist_append.txt");


      //=================================================================//
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
      //===== generate histogram for each of the gabor filtered image
      //================changed on 23rd december ==================//
      
      IplImage* source_img = NULL;
      CvRect* source_roi = new CvRect();
      char* src_filename = new char [50];
      
      char* filename_gsource = new char [50];
      
      strcpy(filename_gsource, "gabor_source_hist.txt");
       char* filename_gsource_app = new char [50];
       strcpy(filename_gsource_app, "gabor_source_hist_app_whole_face.txt");
       
       std::ofstream f_source_w;
       std::ifstream f_source_r;
       int s_data_count =0; 
       float s_data =0.0;
       //====open the file in append mode ===========//
        f_source_w.open(filename_gsource_app, std::ios_base::out | std::ios::app);
      if(f_source_w== NULL)
	{
	  std::cout<<"err in opening file to write "<<filename_gsource <<std::endl;
	  return;
	}
      //====================================//
      
      for(int i =0; i<5; i++) // outer for loop 1
	{

	  for(int j =0; j<8; j++) //inner for loop 1
	    {
	   sprintf(src_filename,"./data/gwt%d_%d.bmp",i,j);
	  // sprintf(filename,"./data/gwt1_6.bmp",i,j);
	  // read each file // get data // keep accumulating data and find
	  // average image 

	  
	   source_img = cvLoadImage(src_filename, 1 );

	   
	   source_roi->x =  0;
	   source_roi->y =  0;
	   source_roi->width =  source_img->width;
	   source_roi->height =  source_img->height;
	   
	  

	   // extract each roi for eyebow1, eyebrow2, lip and nose... //
	   // concatenate the histograms for each of the 4 roi regions --// 
	   
	   histogram_img_roi(source_img, source_roi, filename_gsource);


	   open2read (f_source_r, filename_gsource);

	   while (s_data_count != 59)
	     {

	       f_source_r >> s_data;
	       f_source_w <<s_data <<"\t";
	       s_data_count ++;
      
	     }

	   s_data_count =0;
	   //std::cout<<"gabor_file count = " <<i *j <<std::endl;    

	   f_source_r.close();
	   
	   
	}

    }
    
 f_source_w<<std::endl;
      //==================remove the memory==============//
 
      f_source_w.close();
      delete source_roi;
      delete []filename_gsource;
      delete []filename_gsource_app;
      delete []src_filename;
      cvReleaseImage (&source_img);

      // the above part calculates the lbp features for each of the 40
      // gabor filters and stories the data in contatenation 
//================================================================//
      
       gabor_roi->x = 0;
       gabor_roi->y = 0;
       gabor_roi->height = gabor_avgImg->height;
       gabor_roi->width = gabor_avgImg->width;

       //===read an image of specific scale and orientation ===//
       //==give it as an input to the lbp histogram ======//

       
       histogram_img_roi(gabor_avgImg, gabor_roi, filename_gabor );

     //  std::cout<<"i am here "<<std::endl;

      //=========read file containing histogram of gabor avg image==
      //write the data in append mode in another file =============//

      //==========================================================//
      float fdata_gabor=0.0;
      std::ofstream fw_gabor;
      std::ifstream fr_gabor;
      fw_gabor.open(filename_gabor_app, std::ios_base::out | std::ios::app);
      if(fw_gabor == NULL)
	{
	  std::cout<<"err in opening file to write "<<filename_gabor_app <<std::endl;


	 
	  return;


	}


      open2read (fr_gabor, filename_gabor);

      int gabor_hist_data_count = 0;

      while (gabor_hist_data_count != 59)
	{

	  fr_gabor >> fdata_gabor;
	  fw_gabor <<fdata_gabor <<"\t";
	  gabor_hist_data_count ++;
      
	}

      fw_gabor<<std::endl;




      //==========================================================//
      //releasing memories =======//
       delete gabor_roi;
      fw_gabor.close();
      fr_gabor.close();
      delete [] filename_gabor;
      delete [] filename_gabor_app;
      
    } // end of if condition for gabor based lbp

 
  
  //============this part of the reads data from whole_img_hist.txt ===//
  // ===and write back to another file with each frame data in one
  // row. write in append mode so that all frames data gets collected.
  std::ofstream fw_face;
  std::ifstream fr_face;
  fw_face.open(filename_b, std::ios_base::out | std::ios::app);
  if(fw_face == NULL)
    {
      std::cout<<"err in opening file to write "<<filename_b <<std::endl;

      return;

    }


  open2read (fr_face, filename_a);

  int data_count = 0;

  while (data_count != 59)
    {

      fr_face >> fdata;
      fw_face <<fdata <<"\t";
      data_count ++;
      
    }

  fw_face<<std::endl;
  
  //========remove memory=========//
 
  fw_face.close();
  fr_face.close();
  delete [] filename_a;
  delete [] filename_b;
 
  delete wroi;
  cvReleaseImage (&gray_imag1);
  cvReleaseImage (&face_img);
  cvReleaseImage (&gabor_avgImg);
  return;

}




void face_features::features_lbp(IplImage* img, const bool is_block_wise )
{
 
  IplImage* gray_img = NULL;
  int width =0, height =0;
  // uchar** ppArray = NULL;
  
  struct CvRect *face_rect;
  CvRect *eye1_r;
  CvRect  *eye2_r;
  CvRect *nose;
  face_rect = new CvRect ();
  eye1_r= new CvRect ();
  eye2_r = new CvRect ();
  nose = new CvRect ();
  bool video = false;
  
  // CvPoint * img_size = new CvPoint ();
  CvPoint2D32f *lip_point = new CvPoint2D32f [4]; struct CvRect *lip_rect = new struct CvRect ;	

  IplImage* eye1_img =NULL;
  IplImage* eye2_img =NULL;
  IplImage* lip_img  =NULL;
  IplImage* nose_img =NULL;
 
  char * filename_block = new char [50];
  char * filename_overlap_block = new char [50];
  char * filename1 = new char [50];
  char * filename2 = new char [50];
  char * filename3 = new char [50];
  char * filename4 = new char [50];

  char * filename = new char [50];
 

  width = img->width;//column in an image

  height = img->height;//row in an image

  
  gray_img = cvCreateImage (cvGetSize(img), 8, 1);
 
  //convert color image to gray image
  cvCvtColor(img, gray_img, CV_RGB2GRAY);

  // cvCopyImage(img, gray_img);// crushes here 
  // cvNamedWindow("test", 0); 
  // cvShowImage("test", gray_img) ;

  // cvWaitKey(0);
 
  
  //========detect the face here======================//
  //detect eyes and locate the nose, lips and eyebrow regions==//
  //Apply LBP over each of the regions separately===//
 
  face_sub_regions_detect(gray_img,  face_rect,  eye1_r,  eye2_r,  nose, lip_point, lip_rect);

  // std::cout<<"face rect " <<face_rect->height <<"\t"<<face_rect->width<<std::endl;
  //std::cout<<"i am here "<<std::endl;
  //==apply lbp histogram method over each of the image region separately and write the  
 
  strcpy(filename_block, "block_wise_data_bs_60.txt");
  strcpy(filename_overlap_block, "overlap_block_wise_data_bs_24.txt");
  if(is_block_wise == true)
    {
      // pass only the face image instead of whole image 


      IplImage* face_img = NULL;
     
      //========this is just to try==============//
      //remove some portion of face which are unnecessary noise//
      //std::cout<< "old face rect = "<<face_rect->x<<"\t"<<face_rect->y <<"\t" <<face_rect->height <<"\t" <<face_rect->width <<std::endl;
      face_rect->x = face_rect->x + face_rect->width/10;
      face_rect->y = face_rect->y + face_rect->height/20;
      face_rect->width = face_rect->width - face_rect->width/5;
      face_rect->height = face_rect->height - face_rect->height/10;
      
      //  std::cout<< " face rect = "<<face_rect->x<<"\t"<<face_rect->y <<"\t" <<face_rect->height <<"\t" <<face_rect->width <<std::endl;

      //=========================================//
      face_img = cvCreateImage (cvSize(face_rect->width, face_rect->height), 8, 1);

      
      get_ROI_img(gray_img, face_rect, face_img );

      cvSaveImage("face_img.jpg", face_img);
      IplImage* resiz_img = NULL;
      int new_width = 240;
      int new_height = 240;

      CvRect* ROI_resize = NULL;
      ROI_resize = new CvRect ();
      ROI_resize->x =0;
      ROI_resize->y =0;
      ROI_resize->height = new_height;
      ROI_resize->width = new_width;

      resiz_img = cvCreateImage (cvSize(new_width, new_height), 8, 1);
      
      // resize the face image to make it a const size for all frames so
      // that the data dimension remains the same. 
      resize_img(face_img,  new_width,  new_height, resiz_img );

      //pass the fixed sized face image --i.e resized image into the
      //block wise lbp function and overlapping block wise lbp
      //function. pass the ROI of the resized image into those
      //functions. 
      
      // block_wise_lbp( face_img, face_rect, filename_block );

      // overlapping_block_wise_lbp( face_img, face_rect, filename_overlap_block );

      block_wise_lbp( resiz_img, ROI_resize, filename_block );

      overlapping_block_wise_lbp( resiz_img, ROI_resize, filename_overlap_block );

      cvReleaseImage( &face_img);
      cvReleaseImage (&resiz_img);
      delete ROI_resize;
      
    }
 

  // This commented part is necessary for face sub-part lbp code
  // extraction // 

  strcpy(filename1, "eye1_hist.txt");
  strcpy(filename2, "eye2_hist.txt");
  strcpy(filename3, "nose_hist.txt");
  strcpy(filename4, "lip_hist.txt");
  strcpy(filename, "histogram_data.txt");
  
  if(!is_block_wise==true)
    {
      //=============eye1 feature ==========================//
      eye1_img = cvCreateImage (cvSize(eye1_r->width, eye1_r->height), 8, 1);
      get_ROI_img(gray_img, eye1_r, eye1_img );
      histogram_img_roi(eye1_img, eye1_r, filename1);

      
      //  cvSaveImage("eye1.jpg", eye1_img);
      //=================================================//
      //  std::cout<<" eye1 roi done "<<std::endl;
      //=============eye2 feature ==========================//
      eye2_img = cvCreateImage (cvSize(eye2_r->width, eye2_r->height), 8, 1);
      get_ROI_img(gray_img, eye2_r, eye2_img );
      histogram_img_roi(eye2_img, eye2_r, filename2);
      // cvSaveImage("eye2.jpg", eye2_img);
      //=================================================//
      // std::cout<<" eye2 roi done "<<std::endl;
      //=============nose feature ==========================//
      nose_img = cvCreateImage (cvSize(nose->width, nose->height), 8, 1);
      get_ROI_img(gray_img, nose, nose_img );
      histogram_img_roi(nose_img, nose, filename3);
      // cvSaveImage("nose.jpg", nose_img);
      //=================================================//
      //  std::cout<<" nose roi done "<<std::endl;
      //=============lip feature ==========================//
      //  std::cout<<"lip rect "<<lip_rect->x<<"\t"<<lip_rect->y<<"\t"<<lip_rect->width <<"\t" <<lip_rect->height <<std::endl;
      // std::cout<<"gray img "<<gray_img->width <<"\t" <<gray_img->height <<std::endl;
      lip_img = cvCreateImage (cvSize(lip_rect->width, lip_rect->height), 8, 1);
      get_ROI_img(gray_img, lip_rect, lip_img );
      histogram_img_roi(lip_img, lip_rect, filename4);
      // cvSaveImage("lip.jpg", lip_img);

      //======anima on 22nd october 2013 =======//
      //============find gabor lbp for the subregion
      //======================//
      if( is_gabor_lbp == true)
	{
	  IplImage* gabor_eye1avgImg = NULL;
	  IplImage* gabor_eye2avgImg = NULL;
	  IplImage* gabor_noseavgImg = NULL;
	  IplImage* gabor_lipavgImg = NULL;
	  
	  gabor_eye1avgImg = cvCreateImage(cvSize(IMAGE_SIZE,IMAGE_SIZE), IPL_DEPTH_8U, 1 );
	  gabor_eye2avgImg = cvCreateImage(cvSize(IMAGE_SIZE,IMAGE_SIZE), IPL_DEPTH_8U, 1 );
	  gabor_noseavgImg = cvCreateImage(cvSize(IMAGE_SIZE,IMAGE_SIZE), IPL_DEPTH_8U, 1 );
	  gabor_lipavgImg = cvCreateImage(cvSize(IMAGE_SIZE,IMAGE_SIZE), IPL_DEPTH_8U, 1 );
	  


	  char*  filename_gabor_eye1 = new char [50];
	  char*  filename_gabor_eye2 = new char [50];
	  char*  filename_gabor_nose = new char [50];
	  char*  filename_gabor_lip = new char [50];
	  
	  
	  char*  filename_gabor_eye1_app = new char [50];
	  char*  filename_gabor_eye2_app = new char [50];
	  char*  filename_gabor_nose_app = new char [50];
	  char*  filename_gabor_lip_app = new char [50];
	  char* filename_gabor_all_app = new char [50];
	  
	  strcpy(filename_gabor_eye1, "eye1_gabor_hist.txt");
	  strcpy(filename_gabor_eye2, "eye2_gabor_hist.txt");
	  strcpy(filename_gabor_nose, "nose_gabor_hist.txt");
	  strcpy(filename_gabor_lip, "lip_gabor_hist.txt");
	  
	  strcpy(filename_gabor_eye1_app, "eye1_gabor_hist_app.txt");
	  strcpy(filename_gabor_eye2_app, "eye2_gabor_hist_app.txt");
	  strcpy(filename_gabor_nose_app, "nose_gabor_hist_app.txt");
	  strcpy(filename_gabor_lip_app, "lip_gabor_hist_app.txt");
	  strcpy(filename_gabor_all_app, "gabor_hist_app_4_regions_236_dimn.txt");




	  gabor_based_lbp(eye1_img,  gabor_eye1avgImg,
			  filename_gabor_eye1, filename_gabor_eye1_app );
	  
	  gabor_based_lbp(eye2_img,  gabor_eye2avgImg,
			  filename_gabor_eye2, filename_gabor_eye2_app );
	  
	  gabor_based_lbp(nose_img,  gabor_noseavgImg,
			  filename_gabor_nose, filename_gabor_nose_app );
	  
	  gabor_based_lbp(lip_img,  gabor_lipavgImg,
			  filename_gabor_lip, filename_gabor_lip_app );
	  
	  //this function reads data from the files and writes the data
	  //into a single file ... concatenates the data ====//
	  lbp_file_concatenate(filename_gabor_all_app,
			       filename_gabor_eye1,
			       filename_gabor_eye2,
			       filename_gabor_nose,
			       filename_gabor_lip );


      //=======releasing memories =================//
      cvReleaseImage(&gabor_eye1avgImg);
      cvReleaseImage(&gabor_eye2avgImg );
      cvReleaseImage(&gabor_noseavgImg );
      cvReleaseImage(&gabor_lipavgImg );
      
      delete filename_gabor_eye1;
      delete filename_gabor_eye2;
      delete filename_gabor_nose;
      delete filename_gabor_lip;
      
      delete filename_gabor_eye1_app;
      delete filename_gabor_eye2_app;
      delete filename_gabor_nose_app;
      delete filename_gabor_lip_app;
      delete filename_gabor_all_app;
	}
      //======================================//

  lbp_file_concatenate(filename, filename1, filename2, filename3, filename4 );
      
  /*
  //==========This part of the code is to concatenate the files =========//

  std::ofstream fhist_all;
  std::ifstream fhist_eye1;
  std::ifstream fhist_eye2;
  std::ifstream fhist_nose;
  std::ifstream fhist_lip;

  float x = 0.0;
      
      //  open2write(fhist_all, filename);
      //open to write the file in append mode 
      fhist_all.open(filename, std::ios_base::out | std::ios::app);
      if(fhist_all== NULL)
	{
	  std::cout<<"err in opening file to write "<<filename <<std::endl;
	  return;
	}
 
      open2read (fhist_eye1, filename1);
      open2read (fhist_eye2, filename2);
      open2read (fhist_nose, filename3);
      open2read (fhist_lip, filename4);
   
      //========read all the files and write all data into a file =========//
      //===59 is the size of uniform bin histogram data===// It is done to
      //remove junk data coming in-between -- because of blank space after
      //each file ending. 
      int count =0; 
      while (count != 59)//!fhist_eye1.eof())
	{

	  fhist_eye1>>x;// read the data
	  fhist_all<<x<<"\t"; // write to file 
	  count++;  
	}
  
      //  while (!fhist_eye2.eof())
      count=0;
      while (count != 59)
	{

	  fhist_eye2>>x;
       
	  fhist_all<<x<<"\t";
	  count++;
	}

      count =0;
      //  while (!fhist_nose.eof())
      while (count != 59)
	{

	  fhist_nose>>x;
       
	  fhist_all<<x<<"\t";

	  count++;
	}

      // while (!fhist_lip.eof())
      count =0; 
      while (count != 59)
	{

	  fhist_lip>>x;
       
	  fhist_all<<x<<"\t";

	  count++;
	}
   
    
      fhist_all<<std::endl;


 // ======== close all the files =========//
  fhist_all.close();
  fhist_eye1.close();
  fhist_eye2.close();
  fhist_nose.close();
  fhist_lip.close();
  */
  //---------------------------------------------//

    }// end of if condition // not block wise
  //===comment upto this for extraction only blockwise and overlapping
  //block wise lbp //
  
  //=====release memory====//

  delete []lip_point;
  delete lip_rect;
  delete face_rect;
  delete eye1_r;
  delete eye2_r;
  delete nose;

  delete [] filename1;
  delete [] filename2;
  delete [] filename3;
  delete [] filename4;
  delete [] filename;
  delete [] filename_block;
  delete [] filename_overlap_block;
  //cvReleaseImage (&img);
  cvReleaseImage (&gray_img);
  
  cvReleaseImage (&eye1_img);
  cvReleaseImage (&eye2_img);
  cvReleaseImage (&lip_img);
  cvReleaseImage (&nose_img);

  return;
}

//========================================//


void face_features::block_wise_lbp(IplImage* img, CvRect* const ROI, char* const filename )
{
 
  // say we are going to take block size 8x8
  //check if the ROI is multiple of 8x 8

  uchar data =0;
  int block_size = 60;
  int w_a=0;//for width
  int h_a=0; // for height
  size_t new_width=0, new_height =0;
  IplImage * new_img= NULL;
  IplImage* block_img =NULL;
  CvRect *block_rect = NULL;
  block_rect = new CvRect ();
  uchar** ppArray = NULL;
  //==========memory allocation for the two d array ==============//
  ppArray= new uchar* [block_size];
  for (int i=0; i<block_size; i++){
    
    ppArray[i] = new uchar[block_size]; // each row is having NC column
    if(!ppArray[i]){
      printf("Memory allocation failed ! \n");
      return ;
    }
  }
  //========memory allocation done ==================//

  std::ofstream fhist_block_wise;
  fhist_block_wise.open(filename, std::ios_base::out | std::ios::app);
  if(fhist_block_wise== NULL)
    {
      std::cout<<"err in opening file to write "<<filename <<std::endl;
      return;
    }
 
  //===== This part of the code is used for appending the border pixels
  //to make the image of size that is multiple of block_size
  w_a = ROI->width % block_size; // modulus operator

  h_a = ROI->height % block_size;
  // std::cout<<" ROI = "<<ROI->width <<"\t"<<ROI->height<<std::endl;
  if( w_a !=0)
   
  new_width = ROI->width + (block_size - w_a);

  else
    new_width = ROI->width;
  
  //==================================//
  if( h_a !=0)
 
    new_height = ROI->height + (block_size - h_a);
  else
    new_height = ROI->height;

  new_img =  cvCreateImage (cvSize(new_width, new_height), 8, 1);

  block_img = cvCreateImage (cvSize (block_size, block_size), 8, 1);
  
  for (size_t m=0; m<ROI->height; m++){

    for (size_t n=0; n<ROI->width; n++){
      //read from img and write data into 2D array of new_img
      data = (uchar)CV_IMAGE_ELEM(img, char, m, n);


      cvSetReal2D(new_img, m, n, data );
    }
  
  }
  cvSaveImage("new_img1.jpg", new_img);
  //==================================================//
  for (size_t m= ROI->height; m<new_height; m++){

    for (size_t n= ROI->width; n<new_width; n++){
      //write data into 2D array
      data = (uchar)CV_IMAGE_ELEM(img, char, ROI->height, ROI->width); 

      // data =255;
      cvSetReal2D(new_img, m, n, data );//threshold image
    }
  
  }
  
  //  std::cout<<"new img size = " <<new_img->width <<"\t" <<new_img->height<<std::endl;
  

  //within in two for loops == call the lbp function to read
  //histogram== //

  //=====================================//

  LBP lbp;

  lbp.width = block_size;  lbp.height = block_size;    lbp.tlength = 1;        	
  lbp.R.xR = 1;              lbp.R.yR = 1;	
  lbp.SN.xy = 8;     	
  lbp.uni = 1;               lbp.interp = 1;              lbp.norm = 1;

  double * hist_data = NULL;
  size_t uni_bin_size = 59;
  hist_data = new double[uni_bin_size]; // 58 uni bin and 1 extra
  // bin for non-uniform pattern

  // std::cout<<" u bin = "<< lbp.uni_bin.xy +1 <<std::endl;
    
  //=====================================//
  int cnt =0;
  for( size_t i =0; i< new_img->height;)
    //=================== L1============================    
    {
	
      for( size_t j =0; j<new_img->width;)
	// ===================L2==============================
	{

	  block_rect->x = j;
	  block_rect->y = i;
	  block_rect->width = block_size;
	  block_rect->height = block_size;
	    
	  get_ROI_img(new_img, block_rect, block_img );
	  //========write the data into 2d array== for lbp===//


	  //======================= inner Loop 1 =========================
	  for (size_t m=0; m<block_size; m++){

	    for (size_t n=0; n<block_size; n++){
	      //write data into 2D array
	      ppArray[m][n] = (uchar)CV_IMAGE_ELEM(img, char, m, n);
	    }
	  }
	  //====================== inner loop1 ends======================
	  //====apply lbp for each block =====//

	  //lbp.CreateHistogram(ppArray, 0, 0,  filename);  
	  lbp.Create_block_wise_Histogram( ppArray, 0, 0,  hist_data);
	  
	  // cvNamedWindow("test", CV_WINDOW_AUTOSIZE);
	  // cvShowImage("test", block_img);
	  // cvWaitKey(0);

	  // ================== inner loop2 =========================

	  for( size_t u =0; u<uni_bin_size; u++)
	    {
	      //  std::cout<<"data = "<<hist_data[u] <<std::endl;
	      fhist_block_wise<<hist_data[u]<<"\t"; 
	      
	    }
	  //========================== inner loop 2 ends=========================

	  // j will jump to the next point where there is no overlapping 
	  
	  j= j+block_size;
	  cnt++;
	    
	}
      //=========================== L2 ENDS =======================================
      //i will jump to the next point where there is no overlapping 
      i = i+block_size;
    }
  // =========================== L1 ENDs=========================================


  fhist_block_wise<<std::endl;
  // std::cout<<"cnt = "<<cnt<<std::endl;
  //=========check if the reminder is greater than half of the block
  //size, if so then increase the image size by (block_size - reminder)
  // else decrease by (block_size - reminder)
  
  /*
    
    
    CvPoint* img_size= new CvPoint ();
    img_size->x = ROI->width;
    img_size->y = ROI->height;
    uchar** ppArray = NULL;
    int height = 0, width = 0;
    IplImage* imag1 = NULL;
    cvSaveImage("a.jpg", img);
    width = img_size->x;
    height = img_size->y;
    imag1 = cvCreateImage (cvGetSize(img), 8, 1);
  
    //=====================================//

    LBP lbp;

    lbp.width = width;         lbp.height = height;    	    lbp.tlength = 1;        	
    lbp.R.xR = 1;              lbp.R.yR = 1;	
    lbp.SN.xy = 8;     	
    lbp.uni = 1;               lbp.interp = 1;              lbp.norm = 1;


    //=====================================//
 
 
    ppArray= new uchar* [height];
    for (int i=0; i<height; i++){
    
    ppArray[i] = new uchar[width]; // each row is having NC column
    if(!ppArray[i]){
    printf("Memory allocation failed ! \n");
    return ;
    }
    }

    for (size_t m=0; m<height; m++){

    for (size_t n=0; n<width; n++){
    //write data into 2D array
    ppArray[m][n] = (uchar)CV_IMAGE_ELEM(img, char, m, n);


    cvSetReal2D(imag1, m, n, ppArray[m][n] );//threshold image
    }
  
    }


  
    // char* filename =  new char[50]; // file name should not be more than
    // 50 character 
    // strcpy_s(filename, 20, "basic_hist.txt");
    // strcpy(filename, "basic_hist.txt");
  
    lbp.CreateHistogram(ppArray, 0, 0,  filename);        //lbp based histogram

  

    //======release memories =======//
    // delete [] filename;
    del (ppArray, height);
    delete img_size;
    cvReleaseImage (&imag1);

  */

  fhist_block_wise.close();
    
  delete hist_data;
  delete block_rect;
  del (ppArray, block_size);  
  cvReleaseImage (&new_img);
  cvReleaseImage (&block_img);
  cvDestroyAllWindows();
  return;
}

//==========================================//

void face_features::overlapping_block_wise_lbp(IplImage* img, CvRect* const ROI, char* const filename )
{

  // To make it no overlapping just use the jump = block size
  uchar data =0;
  int block_size = 24;
  int block_jump = 24; // block jump should be less than block size

  int loop_c=0; // Number of times the loop will move 
  int N = ROI->width; // width of the image
  int AP=0; // last allowed point  AP = ROI->width - (block_size-1)
  
  std::cout<<" I am in overlapping block " <<std::endl; 
  int w_a=0;
  int h_a=0;
  size_t new_width=0, new_height =0;
  IplImage * new_img= NULL;
  IplImage* block_img =NULL;
  CvRect *block_rect = NULL;
  block_rect = new CvRect ();
  uchar** ppArray = NULL;
  //==========memory allocation for the two d array ==============//
  ppArray= new uchar* [block_size];
  for (int i=0; i<block_size; i++){
    
    ppArray[i] = new uchar[block_size]; // each row is having NC column
    if(!ppArray[i]){
      printf("Memory allocation failed ! \n");
      return ;
    }
  }
  //========memory allocation done ==================//

  std::ofstream fhist_block_wise;
  fhist_block_wise.open(filename, std::ios_base::out | std::ios::app);
  if(fhist_block_wise== NULL)
    {
      std::cout<<"err in opening file to write "<<filename <<std::endl;
      return;
    }
 
  //========this part calculates the amount of block to be added to make
  //the overlapping window fully slides the image.

  AP = N - (block_size -1);
  loop_c = (AP-1)/block_jump +1 ;

  // check if [(N-1)-{(loop_c -1)x block_jump + (block_size -1)}] == 0 

  if( ( (N-1)- ( (loop_c -1)* block_jump + (block_size -1) ) ) == 0 )

    {
      // no need to add any additional block. i.e w_a = 0

      w_a = 0;

    }

  else

    {

      w_a =  (N-1)- ( (loop_c -1)* block_jump + (block_size -1) );
      
    }
  //=====================================================//
  //====This part makes image a multiple of block size =====//
  // w_a = ROI->width % block_size; // modulus operator

  h_a = ROI->height % block_size;
 
  if (w_a !=0)
    new_width = ROI->width + (block_size - w_a);
  else
    new_width = ROI->width;

  if(h_a !=0)
    new_height = ROI->height + (block_size - h_a);
  else
    new_height = ROI->height;

  //=========above part makes the image multiple of block size=======//

  new_width = new_width + (block_jump -1);
 
  new_img =  cvCreateImage (cvSize(new_width, new_height), 8, 1);

  block_img = cvCreateImage (cvSize (block_size, block_size), 8, 1);
  
  for (size_t m=0; m<ROI->height; m++){

    for (size_t n=0; n<ROI->width; n++){
      //write data into 2D array
      data = (uchar)CV_IMAGE_ELEM(img, char, m, n);


      cvSetReal2D(new_img, m, n, data );//threshold image
    }
  
  }
  cvSaveImage("new_img.jpg", new_img);
  
  //==========================================================//

  LBP lbp;

  lbp.width = block_size;  lbp.height = block_size;    lbp.tlength = 1;        	
  lbp.R.xR = 1;              lbp.R.yR = 1;	
  lbp.SN.xy = 8;     	
  lbp.uni = 1;               lbp.interp = 1;              lbp.norm = 1;

  double * hist_data = NULL;
  size_t uni_bin_size = 59;
  hist_data = new double[uni_bin_size]; // 58 uni bin and 1 extra
  // bin for non-uniform pattern

  // std::cout<<" u bin = "<< lbp.uni_bin.xy +1 <<std::endl;
  int ncount =0;
  // int d_length = new_img->width - (block_size-1);
  // std::cout<<"new img width = "<<new_img->width <<std::endl;
  //=====================================//
   
  for( size_t i =0; i< new_img->height-(block_size-1);){
    for( size_t j =0; j<new_img->width-(block_size-1);){
  
      block_rect->x = j;
      block_rect->y = i;
      block_rect->width = block_size;
      block_rect->height = block_size;
	    
      get_ROI_img(new_img, block_rect, block_img );
      //========write the data into 2d array== for lbp===//

      for (size_t m=0; m<block_size; m++){

	for (size_t n=0; n<block_size; n++){
	  //write data into 2D array
	  ppArray[m][n] = (uchar)CV_IMAGE_ELEM(img, char, m, n);
	}
      }
      //====apply lbp for each block =====//

      //lbp.CreateHistogram(ppArray, 0, 0,  filename);  
      lbp.Create_block_wise_Histogram( ppArray, 0, 0,  hist_data);
	  
      // cvNamedWindow("test", CV_WINDOW_AUTOSIZE);
      // cvShowImage("test", block_img);
      // cvWaitKey(0);

      for( size_t u =0; u<uni_bin_size; u++){

	
	//	  std::cout<<"data = "<<hist_data[u] <<std::endl;
	fhist_block_wise<<hist_data[u]<<"\t";
	ncount++;   
      }

      j= j+block_jump;
	 
    }

    i = i+block_size;
  }
    
  fhist_block_wise<<std::endl;
  // std::cout<<"count = "<<ncount<<std::endl;
    
  //=========check if the reminder is greater than half of the block
  //size, if so then increase the image size by (block_size - reminder)
  // else decrease by (block_size - reminder)
  
  /*
    
    
    CvPoint* img_size= new CvPoint ();
    img_size->x = ROI->width;
    img_size->y = ROI->height;
    uchar** ppArray = NULL;
    int height = 0, width = 0;
    IplImage* imag1 = NULL;
    cvSaveImage("a.jpg", img);
    width = img_size->x;
    height = img_size->y;
    imag1 = cvCreateImage (cvGetSize(img), 8, 1);
  
    //=====================================//

    LBP lbp;

    lbp.width = width;         lbp.height = height;    	    lbp.tlength = 1;        	
    lbp.R.xR = 1;              lbp.R.yR = 1;	
    lbp.SN.xy = 8;     	
    lbp.uni = 1;               lbp.interp = 1;              lbp.norm = 1;


    //=====================================//
 
 
    ppArray= new uchar* [height];
    for (int i=0; i<height; i++){
    
    ppArray[i] = new uchar[width]; // each row is having NC column
    if(!ppArray[i]){
    printf("Memory allocation failed ! \n");
    return ;
    }
    }

    for (size_t m=0; m<height; m++){

    for (size_t n=0; n<width; n++){
    //write data into 2D array
    ppArray[m][n] = (uchar)CV_IMAGE_ELEM(img, char, m, n);


    cvSetReal2D(imag1, m, n, ppArray[m][n] );//threshold image
    }
  
    }


  
    // char* filename =  new char[50]; // file name should not be more than
    // 50 character 
    // strcpy_s(filename, 20, "basic_hist.txt");
    // strcpy(filename, "basic_hist.txt");
  
    lbp.CreateHistogram(ppArray, 0, 0,  filename);        //lbp based histogram

  

    //======release memories =======//
    // delete [] filename;
    del (ppArray, height);
    delete img_size;
    cvReleaseImage (&imag1);

  */

  fhist_block_wise.close();

  delete block_rect;
  delete hist_data;
  del (ppArray, block_size);  
  cvReleaseImage (&new_img);
  cvReleaseImage (&block_img);
  cvDestroyAllWindows();
  return;
}

//====end of overlapping block wise LBP=====//

void face_features::resize_img(IplImage* src_img, const int new_width, const int new_height,IplImage* const& des_img )
{

  cvResize(src_img,des_img,CV_INTER_LINEAR);
  //std::cout<<"des img =" << des_img <<std::endl;
  //  cvSaveImage("img.jpg", des_img);
  return;
}



/*
  void face_features::write_output_data(const int i )
  {


  std::ofstream fyd;
  fyd.open("yd_data.txt", std::ios_base::out | std::ios::app);
  if(fyd== NULL)
  {
  std::cout<<"err in opening file to write output data "<<"yd_data.txt" <<std::endl;
  return;
  }
 

  fyd<<i+1<<std::endl;
  
      
  fyd.close();
  
  return;
  }
*/
//===========================================================//
//===========================================================//

void face_features::gabor_based_lbp(IplImage* src_img, IplImage* gabor_avgImg, char* const filename_gabor, char* const filename_gabor_app )
{
  IplImage* gabor_img= cvCreateImage(cvSize(src_img->width, src_img->height), IPL_DEPTH_8U, 1 );
  cvCopy(src_img, gabor_img);
  GABOR gabor("gabor");
  char ch= 'e';

  gabor.gabor_filtered_data(ch,  gabor_img, gabor_avgImg );
  
  CvRect* gabor_roi = new CvRect ();	
 
  // char* filename_gabor = new char [50];
  // char* filename_gabor_app = new char [50];
  // strcpy(filename_gabor, "gabor_avg_hist.txt");
  // strcpy(filename_gabor_app, "gabor_avg_hist_append.txt");
    
  gabor_roi->x = 0;
  gabor_roi->y = 0;
  gabor_roi->height = gabor_avgImg->height;
  gabor_roi->width = gabor_avgImg->width;
       
  histogram_img_roi(gabor_avgImg, gabor_roi, filename_gabor );

  //=========read file containing histogram of gabor avg image==
  //write the data in append mode in another file =============//
  float fdata_gabor=0.0;
  std::ofstream fw_gabor;
  std::ifstream fr_gabor;
  fw_gabor.open(filename_gabor_app, std::ios_base::out | std::ios::app);
  if(fw_gabor == NULL)
    {
      std::cout<<"err in opening file to write "<<filename_gabor_app <<std::endl;
      return;}

  open2read (fr_gabor, filename_gabor);

  int gabor_hist_data_count = 0;

  while (gabor_hist_data_count != 59)
    {

      fr_gabor >> fdata_gabor;
      fw_gabor <<fdata_gabor <<"\t";
      gabor_hist_data_count ++;
      
    }

  fw_gabor<<std::endl;

  //==========================================================//
  //releasing memories =======//
  delete gabor_roi;
  fw_gabor.close();
  fr_gabor.close();
  // delete [] filename_gabor;
  // delete [] filename_gabor_app;
  
  cvReleaseImage(&gabor_img);
  return;
  
}

//======================================================//

void face_features::lbp_file_concatenate(char* const filename, char* const filename1, char* const filename2, char* const filename3, char* const filename4 )
{

  std::ofstream fhist_all;
  std::ifstream fhist_eye1;
  std::ifstream fhist_eye2;
  std::ifstream fhist_nose;
  std::ifstream fhist_lip;

  float x = 0.0;
      
      //  open2write(fhist_all, filename);
      //open to write the file in append mode 
      fhist_all.open(filename, std::ios_base::out | std::ios::app);
      if(fhist_all== NULL)
	{
	  std::cout<<"err in opening file to write "<<filename <<std::endl;
	  return;
	}
 
      open2read (fhist_eye1, filename1);
      open2read (fhist_eye2, filename2);
      open2read (fhist_nose, filename3);
      open2read (fhist_lip, filename4);
   
      //========read all the files and write all data into a file =========//
      //===59 is the size of uniform bin histogram data===// It is done to
      //remove junk data coming in-between -- because of blank space after
      //each file ending. 
      int count =0; 
      while (count != 59)//!fhist_eye1.eof())
	{

	  fhist_eye1>>x;// read the data
	  fhist_all<<x<<"\t"; // write to file 
	  count++;  
	}
  
      //  while (!fhist_eye2.eof())
      count=0;
      while (count != 59)
	{

	  fhist_eye2>>x;
       
	  fhist_all<<x<<"\t";
	  count++;
	}

      count =0;
      //  while (!fhist_nose.eof())
      while (count != 59)
	{

	  fhist_nose>>x;
       
	  fhist_all<<x<<"\t";

	  count++;
	}

      // while (!fhist_lip.eof())
      count =0; 
      while (count != 59)
	{

	  fhist_lip>>x;
       
	  fhist_all<<x<<"\t";

	  count++;
	}
   
    
      fhist_all<<std::endl;


 // ======== close all the files =========//
  fhist_all.close();
  fhist_eye1.close();
  fhist_eye2.close();
  fhist_nose.close();
  fhist_lip.close();
 
  //---------------------------------------------//



  return;
}


//=============EOF====================//


