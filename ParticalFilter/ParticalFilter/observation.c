/*
 ���ڹ۲ⷽ�̵ĺ�������

  @author Rob Hess
  @version 1.0.0-20060306
*/

#include "defs.h"
#include "utils.h"
#include "observation.h"


/*
    BGRͼ��ת��HSVͼ��

   @param img ��תͼ��

   @return ���� 3ͨ��32λ��HSVͼ�� ɫ��H��ΧΪ  [0,360]
     ���Ͷ�S ������V ȡֵ��ΧΪ [0,1]  
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
  
   ����HSVֱ��ͼ
   
   @param h ɫ��Hue
   @param s ���Ͷ�Saturation
   @param v ����Value
   
   @return ���� HSV��Ӧ��bin index
*/
int histo_bin( float h, float s, float v )
{
  int hd, sd, vd;

  /* ��� S or V ���� ��ֵ, ���� "colorless" bin */
  vd = MIN( (int)(v * NV / V_MAX), NV-1 );
  if( s < S_THRESH  ||  v < V_THRESH )
    return NH * NS + vd;
  
  /* ��֮ ȷ�� "colorful" bin */
  hd = MIN( (int)(h * NH / H_MAX), NH-1 );
  sd = MIN( (int)(s * NS / S_MAX), NS-1 );
  return sd * NH + hd;
}



/*
    ������ͼ����ۼ�ֱ��ͼ
   
   @param imgs ������ͼ��ͼ���Ѿ�ת��HSV
   @param n   ͼ������   
   @return ����δ��һ����ֱ��ͼ
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
      /* ��ͼ������ȡ HSV  */
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
   ��һ��ֱ��ͼ  ��Ϊ1
   
   @param histo ֱ��ͼ
*/
void normalize_histogram( histogram* histo )
{
  float* hist;
  float sum = 0, inv_sum;
  int i, n;

  hist = histo->histo;
  n = histo->n;

  /* ����bins ͨ��sum����������� */
  for( i = 0; i < n; i++ )
    sum += hist[i];
  inv_sum = 1.0 / sum;
  for( i = 0; i < n; i++ )
    hist[i] *= inv_sum;
}



/*
  ����ֱ��ͼ֮���ƽ�������Battacharyyaϵ�� .

   @param h1 δ��һ����ֱ��ͼ
   @param h2 δ��һ����ֱ��ͼ
   
   @return Rerns ����ֱ��ͼ֮���ƽ������ 
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
    ����Battacharyyaϵ��,
    
    D = \sqrt{ 1 - \sum_1^n{ \sqrt{ h_1(i) * h_2(i) } } }
  */
  for( i = 0; i < n; i++ )
    sum += sqrt( hist1[i]*hist2[i] );
  return 1.0 - sum;
}



/*
����Ŀ��ĸ����ܶȷֲ�

@param img HSVͼ��
@param r �������ĵ����λ��
@param c �������ĵ����λ��
@param w Ŀ�������
@param h Ŀ�������
@param ref_histo ��һ���ο�ֱ��ͼ

@return Returns �������ĵ�δΪ��r,c��������ĸ����ܶȷֲ�
*/
float likelihood( IplImage* img, int r, int c,
		  int w, int h, histogram* ref_histo )
{
  IplImage* tmp;
  histogram* histo;
  float d_sq;

  /* ��ȡ���� (r,c) ����͹�һ��ֱ��ͼ */
  cvSetImageROI( img, cvRect( c - w / 2, r - h / 2, w, h ) );
  tmp = cvCreateImage( cvGetSize(img), IPL_DEPTH_32F, 3 );
  cvCopy( img, tmp, NULL );
  cvResetImageROI( img );
  histo = calc_histogram( &tmp, 1 );
  cvReleaseImage( &tmp );
  normalize_histogram( histo );

  /* ��������ܶȷֲ� e^{\lambda D^2(h, h^*)} */
  d_sq = histo_dist_sq( histo, ref_histo );    /*����Battacharyyaϵ���������������ӵ�ֱ��ͼ��ο�ֱ��ͼ�����ƶ�*/
  free( histo );
  return exp( -LAMBDA * d_sq );
}



/*
   ����ͼ��ÿ�����صĵ�λ�ù�һ���ĸ����ܶȷֲ�
   
   @param img HSVͼ��
   @param w Ŀ�������
   @param h Ŀ�������
   @param ref_histo ��һ���ο�ֱ��ͼ
     
   @return ���� ��ͨ��, 32λ ����ͼ��ÿ�����ص�λ�õĹ�һ���ĸ����ܶȷֲ�
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
  ���ֱ��ͼ���ļ���.  �ļ���ʽ����(����gnuplot��ʹ��˵��:

   0 <h_0>
   ...
   <i> <h_i>
   ...
   <n> <h_n>

    n ��ֱ��ͼbins��Ŀ
    h_i, i = 1..n �Ǹ���bin ֵ

   @param histo �������ֱ��ͼ
   @param filename �ļ���

   @return ����1 �ɹ� or 0 ʧ��
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
