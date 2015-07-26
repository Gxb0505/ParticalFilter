/*
  This file contains general program definitions.
  
  @author Rob Hess
  @version 1.0.0-20060306
*/

#ifndef DEFS_H
#define DEFS_H

/********************************* Includes **********************************/

/*  C��׼��ͷ�ļ� */
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <stdarg.h>
#include <stdlib.h>
//#include "unistd.h"

/* OpenCV��ͷ�ļ� */
#include "opencv/cv.h"
#include "opencv/cxcore.h"
#include "opencv/highgui.h"

/******************************* Defs and macros *****************************/

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef MIN                           
#define MIN(x,y) ( ( x < y )? x : y )
#endif
#ifndef MAX
#define MAX(x,y) ( ( x > y )? x : y )
#endif
#ifndef ABS
#define ABS(x) ( ( x < 0 )? -x : x )
#endif

/********************************** Structures *******************************/

#endif
