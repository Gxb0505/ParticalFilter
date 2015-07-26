/*
 用于观测方程的函数定义

  @author Rob Hess
  @version 1.0.0-20060306
*/

#include "defs.h"
#include "utils.h"
#include "observation.h"


/*
    BGR图像转到HSV图像

   @param img 待转图像

   @return 返回 3通道32位的HSV图像 色调H范围为  [0,360]
     饱和度S 和亮度V 取值范围为 [0,1]  
*/
IplImage* bgr2hsv( IplImage* bgr )
{
  IplImage* bgr32f, * hsv;

  bgr32f = cvCreateImage( cvGetSize(bgr), IPL_DEPTH_32F, 3 );
  hsv = cvCreateImage( cvGetSize(bgr), IPL_DEPTH_32F, 3 );
  cvConvertScale( bgr, bgr32f, 1.0 / 255.0, 0 );
  cvCvtColor( bgr32f, hsv, CV_BGR2HSV );
  cvReleaseImage( &bgr32f );
  return hsv;
}


/*
  Calculates the histogram bin into which an HSV entry falls
  
   计算HSV直方图
   
   @param h 色调Hue
   @param s 饱和度Saturation
   @param v 亮度Value
   
   @return 返回 HSV相应的bin index
*/
int histo_bin( float h, float s, float v )
{
  int hd, sd, vd;

  /* 如果 S or V 低于 阈值, 返回 "colorless" bin */
  vd = MIN( (int)(v * NV / V_MAX), NV-1 );
  if( s < S_THRESH  ||  v < V_THRESH )
    return NH * NS + vd;
  
  /* 反之 确定 "colorful" bin */
  hd = MIN( (int)(h * NH / H_MAX), NH-1 );
  sd = MIN( (int)(s * NS / S_MAX), NS-1 );
  return sd * NH + hd;
}



/*
    计算多个图像的累计直方图
   
   @param imgs 待处理图像，图像都已经转成HSV
   @param n   图像数量   
   @return 返回未归一化的直方图
*/
histogram* calc_histogram( IplImage** imgs, int n )
{
  IplImage* img;
  histogram* histo;
  IplImage* h, * s, * v;
  float* hist;
  int i, r, c, bin;

  histo = malloc( sizeof(histogram) );
  histo->n = NH*NS + NV;
  hist = histo->histo;
  memset( hist, 0, histo->n * sizeof(float) );

  for( i = 0; i < n; i++ )
    {
      /* 从图像中提取 HSV  */
      img = imgs[i];
      h = cvCreateImage( cvGetSize(img), IPL_DEPTH_32F, 1 );
      s = cvCreateImage( cvGetSize(img), IPL_DEPTH_32F, 1 );
      v = cvCreateImage( cvGetSize(img), IPL_DEPTH_32F, 1 );
      cvCvtPixToPlane( img, h, s, v, NULL );
      
      /* increment appropriate histogram bin for each pixel */
      for( r = 0; r < img->height; r++ )
	for( c = 0; c < img->width; c++ )
	  {
	    bin = histo_bin( /*pixval32f( h, r, c )*/((float*)(h->imageData + h->widthStep*r) )[c],
			     ((float*)(s->imageData + s->widthStep*r) )[c],
			     ((float*)(v->imageData + v->widthStep*r) )[c] );
	    hist[bin] += 1;
	  }
      cvReleaseImage( &h );
      cvReleaseImage( &s );
      cvReleaseImage( &v );
    }
  return histo;
}



/*
   归一化直方图  和为1
   
   @param histo 直方图
*/
void normalize_histogram( histogram* histo )
{
  float* hist;
  float sum = 0, inv_sum;
  int i, n;

  hist = histo->histo;
  n = histo->n;

  /* 计算bins 通过sum的逆进行连乘 */
  for( i = 0; i < n; i++ )
    sum += hist[i];
  inv_sum = 1.0 / sum;
  for( i = 0; i < n; i++ )
    hist[i] *= inv_sum;
}



/*
  计算直方图之间的平方距离和Battacharyya系数 .

   @param h1 未归一化的直方图
   @param h2 未归一化的直方图
   
   @return Rerns 返回直方图之间的平方距离 
*/
float histo_dist_sq( histogram* h1, histogram* h2 )
{
  float* hist1, * hist2;
  float sum = 0;
  int i, n;

  n = h1->n;
  hist1 = h1->histo;
  hist2 = h2->histo;

  /*
    计算Battacharyya系数,
    
    D = \sqrt{ 1 - \sum_1^n{ \sqrt{ h_1(i) * h_2(i) } } }
  */
  for( i = 0; i < n; i++ )
    sum += sqrt( hist1[i]*hist2[i] );
  return 1.0 - sum;
}



/*
计算目标的概率密度分布

@param img HSV图像
@param r 窗口中心点的行位置
@param c 窗口中心点的列位置
@param w 目标区域宽
@param h 目标区域高
@param ref_histo 归一化参考直方图

@return Returns 返回中心点未为（r,c）的区域的概率密度分布
*/
float likelihood( IplImage* img, int r, int c,
		  int w, int h, histogram* ref_histo )
{
  IplImage* tmp;
  histogram* histo;
  float d_sq;

  /* 提取区域 (r,c) 计算和归一化直方图 */
  cvSetImageROI( img, cvRect( c - w / 2, r - h / 2, w, h ) );
  tmp = cvCreateImage( cvGetSize(img), IPL_DEPTH_32F, 3 );
  cvCopy( img, tmp, NULL );
  cvResetImageROI( img );
  histo = calc_histogram( &tmp, 1 );
  cvReleaseImage( &tmp );
  normalize_histogram( histo );

  /* 计算概率密度分布 e^{\lambda D^2(h, h^*)} */
  d_sq = histo_dist_sq( histo, ref_histo );    /*返回Battacharyya系数，代表任意粒子的直方图与参考直方图的相似度*/
  free( histo );
  return exp( -LAMBDA * d_sq );
}



/*
   返回图像每个像素的点位置归一化的概率密度分布
   
   @param img HSV图像
   @param w 目标区域宽
   @param h 目标区域高
   @param ref_histo 归一化参考直方图
     
   @return 返回 单通道, 32位 浮点图像每个像素点位置的归一化的概率密度分布
*/
IplImage* likelihood_image( IplImage* img, int w, int h, histogram* ref_histo )
{
  IplImage* l, *l2;
  CvScalar sum;
  int i, j;

  l = cvCreateImage( cvGetSize( img ), IPL_DEPTH_32F, 1 );
  for( i = 0; i < img->height; i++ )
    for( j = 0; j < img->width; j++ )
      //setpix32f( l, i, j, likelihood( img, i, j, w, h, ref_histo ) );
	( (float*)(l->imageData + l->widthStep*i) )[j] = likelihood( img, i, j, w, h, ref_histo ) ;

  sum = cvSum( l );
  cvScale( l, l, 1.0 / sum.val[0], 0 );
  return l;
}



/*
  输出直方图到文件中.  文件格式如下(基于gnuplot的使用说明:

   0 <h_0>
   ...
   <i> <h_i>
   ...
   <n> <h_n>

    n 是直方图bins数目
    h_i, i = 1..n 是浮点bin 值

   @param histo 待保存的直方图
   @param filename 文件名

   @return 返回1 成功 or 0 失败
*/
int export_histogram( histogram* histo, char* filename )
{
  int i, n;
  float* h;
  FILE* file = fopen( filename, "w" );

  if( ! file )
    return 0;
  n = histo->n;
  h = histo->histo;
  for( i = 0; i < n; i++ )
    fprintf( file, "%d %f\n", i, h[i] );
  fclose( file );
  return 1;
}
