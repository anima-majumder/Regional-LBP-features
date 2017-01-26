#ifndef FACE_FEATURES
#define FACE_FEATURES

#include "highgui.h"
#include "cv.h"
#include "constants.h"
#include "cvaux.h"

using namespace std;
class face_features

{
public:
  face_features();
  ~face_features();
  void face_sub_regions_detect(IplImage *img, CvRect *const face_rect, CvRect *const eye1_r, CvRect *const eye2_r, CvRect *const nose, CvPoint2D32f *const lip_point,  struct  CvRect  * const  r2);

  void face_detect(IplImage *img, CvRect *const face_rect);
  void resize_img(IplImage* src_img, const int new_width, const int new_height, IplImage* const& des_img );
  
  void get_ROI_img(IplImage* img, CvRect* r, IplImage* ROI_img );
  
  void features_lbp(IplImage* img, const bool is_block_wise);
  void face_lbp(IplImage* img);
  void block_wise_lbp(IplImage* img, CvRect* const ROI, char* const filename );
  void overlapping_block_wise_lbp(IplImage* img, CvRect* const ROI, char* const filename );
  void histogram_img_roi(IplImage* img, CvRect* const ROI, char* const filename );

  void gabor_based_lbp(IplImage* src_img, IplImage* gabor_avgImg, char* const filename_gabor, char* const filename_gabor_app );
  void lbp_file_concatenate(char* const filename, char* const filename1, char* const filename2, char* const filename3, char* const filename4 );
  
private:

  CvRect* lip_rect_prev;
  CvPoint2D32f* lip_position_prev;
  CvRect* face_rect_prev;

  CvRect* eye1r0_prev;
  CvRect* eye2r1_prev;
  CvRect* nose_prev;

};
#endif
