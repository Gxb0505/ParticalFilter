/** @file
    观测方程的参数定义
    
    @author Rob Hess
    @version 1.0.0-20060306
*/

#ifndef OBSERVATION_H
#define OBSERVATION_H

/******************************* 定义 *********************************/

/* HSV颜色直方图bins的数目 */
#define NH 10
#define NS 10
#define NV 10

/* HSV 最大值 */
#define H_MAX 360.0
#define S_MAX 1.0
#define V_MAX 1.0

/* low thresholds on saturation and value for histogramming */
#define S_THRESH 0.1
#define V_THRESH 0.2

/* 分布参数 */
#define LAMBDA 20

/******************************** Structures *********************************/

/**
   An HSV histogram represented by NH * NS + NV bins.  Pixels with saturation
   and value greater than S_THRESH and V_THRESH fill the first NH * NS bins.
   Other, "colorless" pixels fill the last NV value-only bins.
*/
typedef struct histogram {
  float histo[NH*NS + NV];   /**< 直方图数组 */
  int n;                     /**< 数组长度 */
} histogram;


/*************************** 函数定义 ****************************/

/**
   BGR图像转到HSV图像

   @param img 待转图像

   @return 返回 3通道32位的HSV图像 色调H范围为  [0,360]
     饱和度S 和亮度V 取值范围为 [0,1]  
*/
IplImage* bgr2hsv( IplImage* img );


/**
  计算HSV直方图
   
   @param h 色调Hue
   @param s 饱和度Saturation
   @param v 亮度Value
   
   @return 返回 HSV相应的bin index
*/
int histo_bin( float h, float s, float v );


/**
   计算多个图像的累计直方图
   
   @param imgs 待处理图像，图像都已经转成HSV
   @param n   图像数量   
   @return 返回未归一化的直方图
*/
histogram* calc_histogram( IplImage** imgs, int n );


/**
  归一化直方图  和为1
   
   @param histo 直方图
*/
void normalize_histogram( histogram* histo );


/**
   计算直方图之间的平方距离和Battacharyya系数 .

   @param h1 未归一化的直方图
   @param h2 未归一化的直方图
   
   @return Rerns 返回直方图之间的平方距离 
*/
float histo_dist_sq( histogram* h1, histogram* h2 );


/**
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
		  int w, int h, histogram* ref_histo );


/**
   返回图像每个像素的点位置归一化的概率密度分布
   
   @param img HSV图像
   @param w 目标区域宽
   @param h 目标区域高
   @param ref_histo 归一化参考直方图
     
   @return 返回 单通道, 32位 浮点图像每个像素点位置的归一化的概率密度分布
*/
IplImage* likelihood_image( IplImage* img, int w, int h, histogram* ref_histo);


/**
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
int export_histogram( histogram* histo, char* filename );



#endif
