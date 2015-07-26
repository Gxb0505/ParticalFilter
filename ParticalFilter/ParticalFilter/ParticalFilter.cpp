// ParticalFilter.cpp : �������̨Ӧ�ó������ڵ㡣
/*
@author Rob Hess
@version 1.0.0-20060306
*/
/************************************************************************/
/* ����Դ�룬��ӭ����http://www.cnblogs.com/yangyangcv                                                                     */
/************************************************************************/ 
#include "stdafx.h"
#include "defs.h"
#include "utils.h"
#include "particles.h"
#include "observation.h"
#include  "opencv/highgui.h"
#include  "opencv/cv.h"
#include  <opencv2/opencv.hpp>
/******************************** ���� ********************************/


/*������Ŀ*/
#define PARTICLES 100     

/* Ĭ��·�������ͼ����չ�� */
#define EXPORT_BASE "./frames/frame_"
#define EXPORT_EXTN ".png"

/*���֡��*/
#define MAX_FRAMES 2048  

/********************************* �ṹ ********************************/


#define MAX_OBJECTS 1       /*������Ŀ����Ŀ*/

typedef struct params {
  CvPoint loc1[MAX_OBJECTS];
  CvPoint loc2[MAX_OBJECTS];
  IplImage* objects[MAX_OBJECTS];
  char* win_name;
  IplImage* orig_img;
  IplImage* cur_img;
  int n;
} params;


/***************************** ����ԭ�� ***************************/

void usage(char*);                                                        /*��ʾ������ʹ�÷�ʽ��ע�ͣ�*/
//void arg_parse( int, char** );
int get_regions(IplImage*, CvRect**);                                     /*ѡ��Ŀ�꺯��*/
void mouse(int, int, int, int, void*);                                      /*����������*/
void reg(char*, CvRect, void*);                                             /*��txt�ļ���Ϣ��ȡregions*/
histogram** compute_ref_histos(IplImage*, CvRect*, int);                  /*����Ŀ�������ֱ��ͼ*/
int export_ref_histos(histogram**, int);                                  /*���ֱ��ͼ���ļ���*/
//int export_frame( IplImage*, int );


/********************************** Globals **********************************/

char* pname;                                      /*��Ŀ����*/
char* vid_file = "10mf.avi";                    /*��Ƶ����*/
int num_particles = PARTICLES;                    /*������Ŀ*/
int show_all = 0;                                 /*��ʾ���ӿ���flag*/
int export = FALSE;                               /*�����������flag*/


/*********************************** Main ************************************/
int main(int argc, char** argv)
{
  IplImage* frame, *hsv_frame, *frames[MAX_FRAMES];
  IplImage** hsv_ref_imgs;
  histogram** ref_histos;
  CvCapture* video;
  particle* particles, *new_particles;
  CvScalar color;
  CvRect* regions;
  int num_objects = 0;
  float s;
  int i, j, k, w, h, x, y, nn;



  video = cvCaptureFromFile(vid_file);                              /*��ȡ��Ƶ*/ \
    CvVideoWriter *pWriter = cvCreateVideoWriter("results.avi", CV_FOURCC('M', 'J', 'P', 'G'), 12, cvSize(1280, 1024), 1);  /*�����Ƶ*/
  //video = cvCaptureFromCAM(0);
  if (!video)
    fatal_error("couldn't open video file %s", vid_file);

  i = 0;
  while (frame = cvQueryFrame(video))
  {
    hsv_frame = bgr2hsv(frame);                                   /*ת��RGB��HSV��ɫ�ռ�*/
    frames[i] = cvCloneImage(frame);   
    /*�ڵ�һ֡ѡȡ����Ŀ��*/
    if (i == 0)
    {
      w = frame->width;
      h = frame->height;
      fprintf(stderr, "Select object region to track\n");
      while (num_objects == 0)
      {
        num_objects = get_regions(frame, &regions);              /*ѡȡĿ��*/
        if (num_objects == 0)
          fprintf(stderr, "Please select a object\n");
      }

      /* ����ο�ֱ��ͼ �� ��ʼ�����ӷֲ� */
      ref_histos = compute_ref_histos(hsv_frame, regions, num_objects); /*����ο�ֱ��ͼ*/
      if (export)
        export_ref_histos(ref_histos, num_objects);                     /*���ֱ��ͼ���ļ���*/
      particles = init_distribution(regions, ref_histos,                 /*��ʼ�����ӷֲ�*/
        num_objects, num_particles);
    }
    else
    {
      /* ��ÿ������״̬����Ԥ��͹۲� */
      for (j = 0; j < num_particles; j++)
      {
        particles[j] = transition(particles[j], w, h);                /*״̬ת��*/
        s = particles[j].s;
        particles[j].w = likelihood(hsv_frame, cvRound(particles[j].y),   /*��ø����ܶȷֲ�,�õ�Ȩ��*/
          cvRound(particles[j].x),
          cvRound(particles[j].width * s),
          cvRound(particles[j].height * s),
          particles[j].histo);
      }

      /* ��һ��Ȩ�� �� �ز��� */
      normalize_weights(particles, num_particles);                        /*��һ��Ȩ��*/
      new_particles = resample(particles, num_particles);                 /*�ز���*/
      free(particles);
      particles = new_particles;
    }

    /* ��ʾ�������� */
    qsort(particles, num_particles, sizeof(particle), &particle_cmp);    /*��Ȩ�شӸߵ�������*/

    if (show_all)
      for (j = num_particles - 1; j > 0; j--)
      {
        color = CV_RGB(255, 255, 0);
        display_particle(frames[i], particles[j], color);
      }

    /* ��ʾȨ���������� */
    color = CV_RGB(255, 0, 0);
    display_particle(frames[i], particles[0], color);

    nn = cvWriteFrame(pWriter, frames[i]);
    // fprintf("%d/n",nn);
    cvNamedWindow("Video", 1);
    cvShowImage("Video", frames[i]);
    if (cvWaitKey(5) == 27)
      break;
    cvReleaseImage(&hsv_frame);

    i++;
  }
  cvReleaseCapture(&video);
  cvReleaseVideoWriter(&pWriter);
}


/************************** �������� *****************************/

/* ʹ��˵�� */
void usage(char* name)
{
  fprintf(stderr, "%s: track a single object using particle filtering\n\n",
    name);
  fprintf(stderr, "Usage: %s [options] <vid_file>\n\n", name);
  fprintf(stderr, "Arguments:\n");
  fprintf(stderr, "  <vid_file>          A clip of video in which " \
    "to track an object\n");
  fprintf(stderr, "\nOptions:\n");
  fprintf(stderr, "  -h                  Display this message and exit\n");
  fprintf(stderr, "  -a                  Display all particles, not just " \
    "the most likely\n");
  fprintf(stderr, "  -o                  Output tracking sequence frames as " \
    "%s*%s\n", EXPORT_BASE, EXPORT_EXTN);
  fprintf(stderr, "  -p <particles>      Number of particles (default %d)\n",
    PARTICLES);
}



/*
�û����Խ���ѡȡĿ��

@param frame ��Ƶ֡
@param regions ָ��Ŀ����ο��ָ��
@return ����Ŀ����Ŀ
*/
int get_regions(IplImage* frame, CvRect** regions)
{
  char* win_name = "First frame";
  params p;
  CvRect* r;
  int i, x1, y1, x2, y2, w, h,s[4];
  CvRect b;
  bool gotbb_mouse = false;
  bool gotbb_file = false;

  int data, num = 0;
  FILE *fp = fopen("dataInfo4.txt", "r");       /*��txt�ļ��ж�ȡĿ��λ����Ϣ*/
  if (!fp)
  {
    printf("can't open file\n");
    return -1;
  }
 else  gotbb_file = true;
  while (!feof(fp))
  {
    fscanf(fp, "%d\n", &data);
    printf("%d\n", data);           //�ֱ��ȡx,y,w,h�ĸ���Ϣ
    s[num] = data;
    num++;
  }
  printf("\n");
  fclose(fp);

  p.n = 0;
  CvPoint* loc;
  int n = p.n;
  p.loc1[n].x = s[0];
  p.loc1[n].y = s[1];
  w = s[2];
  h = s[3];

  b.x = s[0];
  b.y = s[1];
  b.width = s[2];
  b.height = s[3];
 
  p.win_name = win_name;
  p.orig_img = cvCloneImage(frame);
  p.cur_img = NULL;
  p.n = 0;
  cvNamedWindow(win_name, 1);
  cvShowImage(win_name, frame);
  if (gotbb_file)
  reg(win_name, b, &p);      /*���ļ��л�ȡĿ��regions����*/
  else
  {
    cvSetMouseCallback(win_name, &mouse, &p);   /* ʹ�����ص���������Ŀ������ѡȡ */
    gotbb_mouse = true;
  }
  cvWaitKey(0);
  cvDestroyWindow(win_name);
  cvReleaseImage(&(p.orig_img));
  if (p.cur_img)
    cvReleaseImage(&(p.cur_img));

  /* ��ȡ����; ����Ϊ���ο� */
  if (p.n == 0)
  {
    *regions = NULL;
    return 0;
  }
  r = (CvRect*)malloc(n*sizeof(CvRect));     //����CvRect*���͵�ָ��
  for (i = 0; i < p.n; i++)
  {
    if (gotbb_mouse)
    {
      x1=MIN(p.loc1[i].x, p.loc2[i].x);
      x2=MAX(p.loc1[i].x, p.loc2[i].x);
      y1=MIN(p.loc1[i].y, p.loc2[i].y);
      y2=MAX(p.loc1[i].y, p.loc2[i].y);
       w = x2 - x1;
       h = y2 - y1;
    }
    else if (gotbb_file)
    {
      x1 = p.loc1[i].x - w / 2;
      x2 = p.loc1[i].x + w / 2;
      y1 = p.loc1[i].y - h / 2;
      y2 = p.loc1[i].y + h / 2;
    }    
    /*ȷ����͸�Ϊ���� */
    w = (w % 2) ? w : w + 1;  //���ﲻ��Ϊ��һ��Ҫ����
    h = (h % 2) ? h : h + 1;
    r[i] = cvRect(x1, y1, w, h);
    cvRectangle(p.orig_img, cvPoint(x1, y1), cvPoint(x2, y2), CV_RGB(255, 0, 0), 1, 8, 0); //���ƴ���
  }
  *regions = r;        //��ʼĿ��
  return p.n;
}



/*
���ļ��л�ȡregions�ĺ���
*/
void reg(char* win_name, CvRect s, void* param)
{
  params* p = (params*)param;
  CvPoint* loc;
  int n;
  IplImage* tmp;
  static int pressed = FALSE;
  int w, h;
  /* Ŀ��������Ͻ� */
  w = s.width;
  h = s.height;
  n = p->n;

  loc = p->loc1;
  loc[n].x = s.x - w / 2;
  loc[n].y = s.y- h / 2;


  /* ѡȡ���½ǲ����Ŀ�� */

  loc = p->loc2;
  loc[n].x = s.x + w / 2;
  loc[n].y = s.y + h / 2;
  cvReleaseImage(&(p->cur_img));
  p->cur_img = NULL;
  cvRectangle(p->orig_img, p->loc1[n], p->loc2[n], CV_RGB(255, 0, 0), 1, 8, 0);
  cvShowImage(p->win_name, p->orig_img);
  p->n++;
}

/*
���ص�����
*/
void mouse(int event, int x, int y, int flags, void* param)
{
  params* p = (params*)param;
  CvPoint* loc;
  int n;
  IplImage* tmp;
  static int pressed = FALSE;

  /* on left button press, remember first corner of rectangle around object */
  if (event == CV_EVENT_LBUTTONDOWN)
  {
    n = p->n;
    if (n == MAX_OBJECTS)
      return;
    loc = p->loc1;
    loc[n].x = x;
    loc[n].y = y;
    pressed = TRUE;
  }

  /* on left button up, finalize the rectangle and draw it in black */
  else if (event == CV_EVENT_LBUTTONUP)
  {
    n = p->n;
    if (n == MAX_OBJECTS)
      return;
    loc = p->loc2;
    loc[n].x = x;
    loc[n].y = y;
    cvReleaseImage(&(p->cur_img));
    p->cur_img = NULL;
    cvRectangle(p->orig_img, p->loc1[n], loc[n], CV_RGB(0, 0, 0), 1, 8, 0);
    cvShowImage(p->win_name, p->orig_img);
    pressed = FALSE;
    p->n++;
  }

  /* on mouse move with left button down, draw rectangle as defined in white */
  else if (event == CV_EVENT_MOUSEMOVE  &&  flags & CV_EVENT_FLAG_LBUTTON)
  {
    n = p->n;
    if (n == MAX_OBJECTS)
      return;
    tmp = cvCloneImage(p->orig_img);
    loc = p->loc1;
    cvRectangle(tmp, loc[n], cvPoint(x, y), CV_RGB(255, 255, 255), 1, 8, 0);
    cvShowImage(p->win_name, tmp);
    if (p->cur_img)
      cvReleaseImage(&(p->cur_img));
    p->cur_img = tmp;
  }
}


/*
����ѡȡ����Ŀ��ο�ֱ��ͼ

@param frame   ͼ��
@param regions Ŀ������
@param n       ��������
@param export if TRUE, ͼ�������־

@return        ֱ��ͼ
*/
histogram** compute_ref_histos(IplImage* frame, CvRect* regions, int n)              /*����Ŀ��ο�ֱ��ͼ*/
{
  histogram** histos = (histogram**)malloc(n * sizeof(histogram*));                            /*�����ڴ�*/
  IplImage* tmp;
  int i;

  /* extract each region from frame and compute its histogram */
  for (i = 0; i < n; i++)
  {
    cvSetImageROI(frame, regions[i]);                                        /*����ROI����*/
    tmp = cvCreateImage(cvGetSize(frame), IPL_DEPTH_32F, 3);
    cvCopy(frame, tmp, NULL);
    cvResetImageROI(frame);
    histos[i] = calc_histogram(&tmp, 1);                                    /*����ֱ��ͼ*/
    normalize_histogram(histos[i]);                                         /*ֱ��ͼ��һ��*/
    cvReleaseImage(&tmp);
  }

  return histos;
}



/*
����ο�ֱ��ͼ���ļ�

@param ref_histos �ο�ֱ��ͼ
@param n          ֱ��ͼ��Ŀ
@return           ���� 1  �ɹ� or 0  ʧ��
*/
int export_ref_histos(histogram** ref_histos, int n)
{
  char name[32];
  char num[3];
  FILE* file;
  int i;

  for (i = 0; i < n; i++)
  {
    sprintf(num, "%02d", i);
    strcpy(name, "hist_");
    strcat(name, num);
    strcat(name, ".dat");
    if (!export_histogram(ref_histos[i], name))
      return 0;
  }

  return 1;
}



