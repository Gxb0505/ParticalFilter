// ParticalFilter.cpp : 定义控制台应用程序的入口点。
/*
@author Rob Hess
@version 1.0.0-20060306
*/
/************************************************************************/
/* 更多源码，欢迎访问http://www.cnblogs.com/yangyangcv                                                                     */
/************************************************************************/ 
#include "stdafx.h"
#include "defs.h"
#include "utils.h"
#include "particles.h"
#include "observation.h"
#include  "opencv/highgui.h"
#include  "opencv/cv.h"
#include  <opencv2/opencv.hpp>
/******************************** 定义 ********************************/


/*粒子数目*/
#define PARTICLES 100     

/* 默认路径和输出图像扩展名 */
#define EXPORT_BASE "./frames/frame_"
#define EXPORT_EXTN ".png"

/*最大帧数*/
#define MAX_FRAMES 2048  

/********************************* 结构 ********************************/


#define MAX_OBJECTS 1       /*最大跟踪目标数目*/

typedef struct params {
  CvPoint loc1[MAX_OBJECTS];
  CvPoint loc2[MAX_OBJECTS];
  IplImage* objects[MAX_OBJECTS];
  char* win_name;
  IplImage* orig_img;
  IplImage* cur_img;
  int n;
} params;


/***************************** 函数原型 ***************************/

void usage(char*);                                                        /*显示本工程使用方式（注释）*/
//void arg_parse( int, char** );
int get_regions(IplImage*, CvRect**);                                     /*选择目标函数*/
void mouse(int, int, int, int, void*);                                      /*鼠标操作函数*/
void reg(char*, CvRect, void*);                                             /*从txt文件信息获取regions*/
histogram** compute_ref_histos(IplImage*, CvRect*, int);                  /*计算目标区域的直方图*/
int export_ref_histos(histogram**, int);                                  /*输出直方图到文件中*/
//int export_frame( IplImage*, int );


/********************************** Globals **********************************/

char* pname;                                      /*项目名称*/
char* vid_file = "10mf.avi";                    /*视频名称*/
int num_particles = PARTICLES;                    /*粒子数目*/
int show_all = 0;                                 /*显示粒子控制flag*/
int export = FALSE;                               /*输出跟踪序列flag*/


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



  video = cvCaptureFromFile(vid_file);                              /*读取视频*/ \
    CvVideoWriter *pWriter = cvCreateVideoWriter("results.avi", CV_FOURCC('M', 'J', 'P', 'G'), 12, cvSize(1280, 1024), 1);  /*输出视频*/
  //video = cvCaptureFromCAM(0);
  if (!video)
    fatal_error("couldn't open video file %s", vid_file);

  i = 0;
  while (frame = cvQueryFrame(video))
  {
    hsv_frame = bgr2hsv(frame);                                   /*转换RGB到HSV颜色空间*/
    frames[i] = cvCloneImage(frame);   
    /*在第一帧选取跟踪目标*/
    if (i == 0)
    {
      w = frame->width;
      h = frame->height;
      fprintf(stderr, "Select object region to track\n");
      while (num_objects == 0)
      {
        num_objects = get_regions(frame, &regions);              /*选取目标*/
        if (num_objects == 0)
          fprintf(stderr, "Please select a object\n");
      }

      /* 计算参考直方图 和 初始化粒子分布 */
      ref_histos = compute_ref_histos(hsv_frame, regions, num_objects); /*计算参考直方图*/
      if (export)
        export_ref_histos(ref_histos, num_objects);                     /*输出直方图到文件中*/
      particles = init_distribution(regions, ref_histos,                 /*初始化粒子分布*/
        num_objects, num_particles);
    }
    else
    {
      /* 对每个粒子状态进行预测和观测 */
      for (j = 0; j < num_particles; j++)
      {
        particles[j] = transition(particles[j], w, h);                /*状态转移*/
        s = particles[j].s;
        particles[j].w = likelihood(hsv_frame, cvRound(particles[j].y),   /*获得概率密度分布,得到权重*/
          cvRound(particles[j].x),
          cvRound(particles[j].width * s),
          cvRound(particles[j].height * s),
          particles[j].histo);
      }

      /* 归一化权重 和 重采样 */
      normalize_weights(particles, num_particles);                        /*归一化权重*/
      new_particles = resample(particles, num_particles);                 /*重采样*/
      free(particles);
      particles = new_particles;
    }

    /* 显示所有粒子 */
    qsort(particles, num_particles, sizeof(particle), &particle_cmp);    /*按权重从高到底排列*/

    if (show_all)
      for (j = num_particles - 1; j > 0; j--)
      {
        color = CV_RGB(255, 255, 0);
        display_particle(frames[i], particles[j], color);
      }

    /* 显示权重最大的粒子 */
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


/************************** 函数定义 *****************************/

/* 使用说明 */
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
用户可以交互选取目标

@param frame 视频帧
@param regions 指向目标矩形框的指针
@return 返回目标数目
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
  FILE *fp = fopen("dataInfo4.txt", "r");       /*从txt文件中读取目标位置信息*/
  if (!fp)
  {
    printf("can't open file\n");
    return -1;
  }
 else  gotbb_file = true;
  while (!feof(fp))
  {
    fscanf(fp, "%d\n", &data);
    printf("%d\n", data);           //分别读取x,y,w,h四个信息
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
  reg(win_name, b, &p);      /*从文件中获取目标regions矩形*/
  else
  {
    cvSetMouseCallback(win_name, &mouse, &p);   /* 使用鼠标回调函数进行目标区域选取 */
    gotbb_mouse = true;
  }
  cvWaitKey(0);
  cvDestroyWindow(win_name);
  cvReleaseImage(&(p.orig_img));
  if (p.cur_img)
    cvReleaseImage(&(p.cur_img));

  /* 提取区域; 保存为矩形框 */
  if (p.n == 0)
  {
    *regions = NULL;
    return 0;
  }
  r = (CvRect*)malloc(n*sizeof(CvRect));     //返回CvRect*类型的指针
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
    /*确保宽和高为奇数 */
    w = (w % 2) ? w : w + 1;  //这里不懂为何一定要奇数
    h = (h % 2) ? h : h + 1;
    r[i] = cvRect(x1, y1, w, h);
    cvRectangle(p.orig_img, cvPoint(x1, y1), cvPoint(x2, y2), CV_RGB(255, 0, 0), 1, 8, 0); //绘制窗口
  }
  *regions = r;        //初始目标
  return p.n;
}



/*
从文件中获取regions的函数
*/
void reg(char* win_name, CvRect s, void* param)
{
  params* p = (params*)param;
  CvPoint* loc;
  int n;
  IplImage* tmp;
  static int pressed = FALSE;
  int w, h;
  /* 目标矩形左上角 */
  w = s.width;
  h = s.height;
  n = p->n;

  loc = p->loc1;
  loc[n].x = s.x - w / 2;
  loc[n].y = s.y- h / 2;


  /* 选取右下角并获得目标 */

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
鼠标回调函数
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
计算选取跟踪目标参考直方图

@param frame   图像
@param regions 目标区域
@param n       区域数量
@param export if TRUE, 图像输出标志

@return        直方图
*/
histogram** compute_ref_histos(IplImage* frame, CvRect* regions, int n)              /*计算目标参考直方图*/
{
  histogram** histos = (histogram**)malloc(n * sizeof(histogram*));                            /*分配内存*/
  IplImage* tmp;
  int i;

  /* extract each region from frame and compute its histogram */
  for (i = 0; i < n; i++)
  {
    cvSetImageROI(frame, regions[i]);                                        /*设置ROI区域*/
    tmp = cvCreateImage(cvGetSize(frame), IPL_DEPTH_32F, 3);
    cvCopy(frame, tmp, NULL);
    cvResetImageROI(frame);
    histos[i] = calc_histogram(&tmp, 1);                                    /*计算直方图*/
    normalize_histogram(histos[i]);                                         /*直方图归一化*/
    cvReleaseImage(&tmp);
  }

  return histos;
}



/*
输出参考直方图到文件

@param ref_histos 参考直方图
@param n          直方图数目
@return           返回 1  成功 or 0  失败
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



