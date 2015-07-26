/** @file
    Definitions related to tracking with particle filtering

    @author Rob Hess
    @version 1.0.0-20060307
*/

#ifndef PARTICLES_H
#define PARTICLES_H

#include "observation.h"

/******************************* ���� *********************************/

/* ״̬ת�Ʒ��������ڸ�˹�����ı�׼�� */
#define TRANS_X_STD 1.0
#define TRANS_Y_STD 0.5
#define TRANS_S_STD 0.001

/* ״̬ת�Ʒ����Իع鶯̬���� */
#define A1  2.0
#define A2 -1.0
#define B0  1.0000

/******************************* Structures **********************************/

/**
   �����Ǽ��ϵͳ��ʵʱ״̬����  ���ӵļ���ʵ����ϵͳ������ʵ���ɢ��
*/
typedef struct particle {
  float x;          /**< ��ǰx���� */
  float y;          /**< ��ǰy���� */
  float s;          /**< �߶� */
  float xp;         /**< ǰһ��x���� */
  float yp;         /**< ǰһ��y����*/
  float sp;         /**< ֮ǰ�ĳ߶� */
  float x0;         /**< ��ʼx���� */
  float y0;         /**< ��ʼy���� */
  int width;        /**< ��ʼĿ�������� */
  int height;       /**< ��ʼĿ������߶� */
  histogram* histo; /**< ��������Ŀ������ο�ֱ��ͼ */
  float w;          /**< Ȩ�� */
} particle;


/**************************** ����ԭ�� ****************************/

/**
   Creates an initial distribution of particles by sampling from a Gaussian
   window around each of a set of specified locations
   
   
   @param regions ����Ŀ������
   @param histos Ŀ������ֱ��ͼ
   @param n ������Ŀ
   @param p ������Ŀ
   
   @return ���ز�������
*/
particle* init_distribution(CvRect* regions, histogram** histos, int n, int p);


/**
   ״̬ת�Ʒ���

   @param p ����
   @param w ֡��
   @param h ֡��
   @param rng a random number generator from which to sample

   @return ���� <EM>p</EM>'״̬ת�Ʒ��̷����µ����� 
*/
particle transition( particle p, int w, int h);


/**
   ��һ��Ȩ�� ��Ϊ1

   @param particles��һ��Ȩ�ص�����
   @param n ������Ŀ
*/
void normalize_weights( particle* particles, int n );


/**
   �ز��������µ����� 

   @param particles �Ѿ�Ȩ�ع�һ��������
   @param n ������Ŀ
  
   @return ����δȨ�ع�һ��������
*/
particle* resample( particle* particles, int n );


/**
     ����Ȩ�رȽ�����.  

  @param p1 ����ָ��
  @param p2 ����ָ��

  @return ���� -1: p1 Ȩ��С�� p2; 1 :p1 Ȩ�ش��� p2; 0: Ȩ����ͬ.
*/
int particle_cmp(const  void* p1,const  void* p2 );


/**
  ��ͼ������ʾ����

   @param img ͼ��
   @param p ����
   @param color ��ɫ
*/
void display_particle( IplImage* img, particle p, CvScalar color);

double GaussRand(void);

#endif
