/** @file
    Definitions related to tracking with particle filtering

    @author Rob Hess
    @version 1.0.0-20060307
*/

#ifndef PARTICLES_H
#define PARTICLES_H

#include "observation.h"

/******************************* 定义 *********************************/

/* 状态转移方程中用于高斯采样的标准差 */
#define TRANS_X_STD 1.0
#define TRANS_Y_STD 0.5
#define TRANS_S_STD 0.001

/* 状态转移方程自回归动态参数 */
#define A1  2.0
#define A2 -1.0
#define B0  1.0000

/******************************* Structures **********************************/

/**
   粒子是监测系统的实时状态变量  粒子的集合实质是系统后验概率的离散化
*/
typedef struct particle {
  float x;          /**< 当前x坐标 */
  float y;          /**< 当前y坐标 */
  float s;          /**< 尺度 */
  float xp;         /**< 前一刻x坐标 */
  float yp;         /**< 前一刻y坐标*/
  float sp;         /**< 之前的尺度 */
  float x0;         /**< 初始x坐标 */
  float y0;         /**< 初始y坐标 */
  int width;        /**< 初始目标区域宽度 */
  int height;       /**< 初始目标区域高度 */
  histogram* histo; /**< 描述跟踪目标区域参考直方图 */
  float w;          /**< 权重 */
} particle;


/**************************** 函数原型 ****************************/

/**
   Creates an initial distribution of particles by sampling from a Gaussian
   window around each of a set of specified locations
   
   
   @param regions 跟踪目标区域
   @param histos 目标区域直方图
   @param n 区域数目
   @param p 粒子数目
   
   @return 返回采样粒子
*/
particle* init_distribution(CvRect* regions, histogram** histos, int n, int p);


/**
   状态转移方程

   @param p 粒子
   @param w 帧宽
   @param h 帧高
   @param rng a random number generator from which to sample

   @return 基于 <EM>p</EM>'状态转移方程返回新的粒子 
*/
particle transition( particle p, int w, int h);


/**
   归一化权重 和为1

   @param particles归一化权重的粒子
   @param n 粒子数目
*/
void normalize_weights( particle* particles, int n );


/**
   重采样产生新的粒子 

   @param particles 已经权重归一化的粒子
   @param n 粒子数目
  
   @return 返回未权重归一化的粒子
*/
particle* resample( particle* particles, int n );


/**
     基于权重比较粒子.  

  @param p1 粒子指针
  @param p2 粒子指针

  @return 返回 -1: p1 权重小于 p2; 1 :p1 权重大于 p2; 0: 权重相同.
*/
int particle_cmp(const  void* p1,const  void* p2 );


/**
  在图像中显示粒子

   @param img 图像
   @param p 粒子
   @param color 颜色
*/
void display_particle( IplImage* img, particle p, CvScalar color);

double GaussRand(void);

#endif
