/** @file
    �۲ⷽ�̵Ĳ�������
    
    @author Rob Hess
    @version 1.0.0-20060306
*/

#ifndef OBSERVATION_H
#define OBSERVATION_H

/******************************* ���� *********************************/

/* HSV��ɫֱ��ͼbins����Ŀ */
#define NH 10
#define NS 10
#define NV 10

/* HSV ���ֵ */
#define H_MAX 360.0
#define S_MAX 1.0
#define V_MAX 1.0

/* low thresholds on saturation and value for histogramming */
#define S_THRESH 0.1
#define V_THRESH 0.2

/* �ֲ����� */
#define LAMBDA 20

/******************************** Structures *********************************/

/**
   An HSV histogram represented by NH * NS + NV bins.  Pixels with saturation
   and value greater than S_THRESH and V_THRESH fill the first NH * NS bins.
   Other, "colorless" pixels fill the last NV value-only bins.
*/
typedef struct histogram {
  float histo[NH*NS + NV];   /**< ֱ��ͼ���� */
  int n;                     /**< ���鳤�� */
} histogram;


/*************************** �������� ****************************/

/**
   BGRͼ��ת��HSVͼ��

   @param img ��תͼ��

   @return ���� 3ͨ��32λ��HSVͼ�� ɫ��H��ΧΪ  [0,360]
     ���Ͷ�S ������V ȡֵ��ΧΪ [0,1]  
*/
IplImage* bgr2hsv( IplImage* img );


/**
  ����HSVֱ��ͼ
   
   @param h ɫ��Hue
   @param s ���Ͷ�Saturation
   @param v ����Value
   
   @return ���� HSV��Ӧ��bin index
*/
int histo_bin( float h, float s, float v );


/**
   ������ͼ����ۼ�ֱ��ͼ
   
   @param imgs ������ͼ��ͼ���Ѿ�ת��HSV
   @param n   ͼ������   
   @return ����δ��һ����ֱ��ͼ
*/
histogram* calc_histogram( IplImage** imgs, int n );


/**
  ��һ��ֱ��ͼ  ��Ϊ1
   
   @param histo ֱ��ͼ
*/
void normalize_histogram( histogram* histo );


/**
   ����ֱ��ͼ֮���ƽ�������Battacharyyaϵ�� .

   @param h1 δ��һ����ֱ��ͼ
   @param h2 δ��һ����ֱ��ͼ
   
   @return Rerns ����ֱ��ͼ֮���ƽ������ 
*/
float histo_dist_sq( histogram* h1, histogram* h2 );


/**
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
		  int w, int h, histogram* ref_histo );


/**
   ����ͼ��ÿ�����صĵ�λ�ù�һ���ĸ����ܶȷֲ�
   
   @param img HSVͼ��
   @param w Ŀ�������
   @param h Ŀ�������
   @param ref_histo ��һ���ο�ֱ��ͼ
     
   @return ���� ��ͨ��, 32λ ����ͼ��ÿ�����ص�λ�õĹ�һ���ĸ����ܶȷֲ�
*/
IplImage* likelihood_image( IplImage* img, int w, int h, histogram* ref_histo);


/**
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
int export_histogram( histogram* histo, char* filename );



#endif
