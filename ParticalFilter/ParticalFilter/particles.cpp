/*
Functions for object tracking with a particle filter

@author Rob Hess
@version 1.0.0-20060310
*/

#include "stdafx.h"
#include "defs.h"
#include "utils.h"
#include "particles.h"
#include<opencv/cv.h>

/*************************** �������� ****************************/

/*

����ָ��λ�õ����ӳ�ʼ�ֲ�
@param regions Ŀ������
@param histos ����Ŀ������ֱ��ͼ
@param n Ŀ��������Ŀ
@param p ��������

@return ����������Χ�ķֲ�����
*/
particle* init_distribution(CvRect* regions, histogram** histos, int n, int p)
{
  particle* particles;
  int np;
  float x, y;
  int i, j, width, height, k = 0;

  particles = (particle*)malloc(p * sizeof(particle));
  np = p / n;

  /* ����n��ģ������ĵ����� */
  for (i = 0; i < n; i++)
  {
    width = regions[i].width;
    height = regions[i].height;
    x = regions[i].x + width / 2;                          /*Ŀ�����ĵ������*/
    y = regions[i].y + height / 2;                         /*Ŀ�����ĵ�������*/
    for (j = 0; j < np; j++)
    {
      particles[k].x0 = particles[k].xp = particles[k].x = x;     //���ﶼ�ǳ�ʼ����һ����ĸ�������
      particles[k].y0 = particles[k].yp = particles[k].y = y;
      particles[k].sp = particles[k].s = 1.0;
      particles[k].width = width;
      particles[k].height = height;
      particles[k].histo = histos[i];
      particles[k++].w = 0;
    }
  }

  /* ���� p �����ӣ�������p���� */
  i = 0;
  while (k < p)
  {
    width = regions[i].width;
    height = regions[i].height;
    x = regions[i].x + width / 2;
    y = regions[i].y + height / 2;
    particles[k].x0 = particles[k].xp = particles[k].x = x;
    particles[k].y0 = particles[k].yp = particles[k].y = y;
    particles[k].sp = particles[k].s = 1.0;
    particles[k].width = width;
    particles[k].height = height;
    particles[k].histo = histos[i];
    particles[k++].w = 0;
    i = (i + 1) % n;
  }

  return particles;
}



/*
ת�Ʒ���

@param p ����
@param w ֡��
@param h ֡��
@param rng a random number generator from which to sample

@return ���ػ��� <EM>p</EM>'s ת�Ʒ��̲����õ���������
*/
particle transition(particle p, int w, int h)
{
  float x, y, s;
  particle pn;

  /* ʹ�ö��׶�̬�Իع�õ��µĲ���״̬ */
  x = A1 * (p.x - p.x0) + A2 * (p.xp - p.x0) +
    B0 * GaussRand()* TRANS_X_STD + p.x0;
  pn.x = MAX(0.0, MIN((float)w - 1.0, x));
  y = A1 * (p.y - p.y0) + A2 * (p.yp - p.y0) +
    B0 * GaussRand()* TRANS_Y_STD + p.y0;
  y = A1 * (p.y - p.y0) + A2 * (p.yp - p.y0) +
    B0 *  +p.y0;
  pn.y = MAX(0.0, MIN((float)h - 1.0, y));
  s = A1 * (p.s - 1.0) + A2 * (p.sp - 1.0) +
    B0 * GaussRand()* TRANS_S_STD + 1.0;
  pn.s = MAX(0.1, s);
  pn.xp = p.x;
  pn.yp = p.y;
  pn.sp = p.s;
  pn.x0 = p.x0;
  pn.y0 = p.y0;
  pn.width = p.width;
  pn.height = p.height;
  pn.histo = p.histo;
  pn.w = 0;

  return pn;
}



/*
��һ������Ȩ�� ���Ϊ1

@param particles ����
@param n ������Ŀ
*/
void normalize_weights(particle* particles, int n)
{
  float sum = 0;
  int i;

  for (i = 0; i < n; i++)
    sum += particles[i].w;
  for (i = 0; i < n; i++)
    particles[i].w /= sum;
}



/*
����Ȩ�طֲ������ز��� ����������

@param particles ������
@param n ������Ŀ
@return ����������
*/
particle* resample(particle* particles, int n)
{
  particle* new_particles;
  int i, j, np, k = 0;

  qsort(particles, n, sizeof(particle), &particle_cmp);     /*����*/
  new_particles = (particle*)malloc(n * sizeof(particle));
  for (i = 0; i < n; i++)
  {
    np = cvRound(particles[i].w * n);     //np�ǵ�ǰ����Ȩ������������Ŀ�ĳ˻�
    for (j = 0; j < np; j++)
    {
      new_particles[k++] = particles[i];  //ǰnp�����Ӷ��ǵ�ǰi�����ӣ����ֵ�˼����������Ȩ�ظߵ������������Ŀ����Ϊ����
      if (k == n)
        goto exit;
    }
  }
  while (k < n)
    new_particles[k++] = particles[0];   //���k��ȻС��n����ʣ������Ӷ����1������������ͬ 
                                         //�����n��������Ŀ����ǰ��ĸ���Ŀ����Ŀ�ǲ�ͬ�ĸ���

exit:
  return new_particles;
}



/*
����Ȩ�رȽ�����.

@param p1 ����ָ��
@param p2 ����ָ��

@return ���� -1: p1 Ȩ��С�� p2; 1 :p1 Ȩ�ش��� p2; 0: Ȩ����ͬ.
*/
int particle_cmp(const void* p1, const void* p2)
{
  particle* _p1 = (particle*)p1;
  particle* _p2 = (particle*)p2;

  if (_p1->w > _p2->w)
    return -1;
  if (_p1->w < _p2->w)
    return 1;
  return 0;
}



/*
��ͼ������ʾ����

@param img ͼ��
@param p ����
@param color ��ɫ
*/
void display_particle(IplImage* img, particle p, CvScalar color)
{
  int x0, y0, x1, y1;

  x0 = cvRound(p.x - 0.5 * p.s * p.width);
  y0 = cvRound(p.y - 0.5 * p.s * p.height);
  x1 = x0 + cvRound(p.s * p.width);
  y1 = y0 + cvRound(p.s * p.height);

  cvRectangle(img, cvPoint(x0, y0), cvPoint(x1, y1), color, 1, 8, 0);
}
/*������˹�����*/
double GaussRand(void)
{
  static double V2, fac;
  static int phase = 0;
  double S, Z, U1, U2, V1;
  if (phase){
    Z = V2 * fac;
  }
  else {
    do {
      U1 = (double)rand() / RAND_MAX;
      U2 = (double)rand() / RAND_MAX;
      V1 = 2 * U1 - 1;
      V2 = 2 * U2 - 1;
      S = V1 * V1 + V2 * V2;
    } while (S >= 1);
    fac = sqrt(-2 * log(S) / S);
    Z = V1 * fac;
  }
  phase = 1 - phase;
  return Z;
}
