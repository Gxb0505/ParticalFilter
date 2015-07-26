/*
  utils.c

  该c文件定义了各式各样实用的函数.
*/

#include "utils.h"
#include "defs.h"


/***************************** 内联函数 ******************************/

/* 返回32位浮点图像的一个像素值 */
// float pixval32f(IplImage* img, int r, int c)
// {
//   return ( (float*)(img->imageData + img->widthStep*r) )[c];
// }


/* 设置64位浮点图像像素值  */
// void setpix32f(IplImage* img, int r, int c, float val)
// {
//   ( (float*)(img->imageData + img->widthStep*r) )[c] = val;
// }

/*************************** 函数定义 ****************************/

/*
  fatal_error() 打印格式错误信息 并退出
*/
void fatal_error(char* format, ...)
{
  va_list ap;  // 格式参数列表
  
  fprintf( stderr, "Error: ");

  // assemble argument list and print format using it
  va_start( ap, format );
  vfprintf( stderr, format, ap );
  va_end(ap);
  printf("\n");
  exit(1);
}



/*
  replace_extension() 用新的扩展名替换文件扩展名 并返回结果.  
*/
char* replace_extension(const char* file, const char* new_extn)
{
  char* new_file = malloc( (strlen(file)+strlen(new_extn)+2) * sizeof(char) );
  char* lastdot;   // 扩展位置

  // 复制文件名，找到扩展位置
  strcpy( new_file, file );
  lastdot = strrchr( new_file, '.' );

  // 如果文件包含扩展名, 去除掉,附上新的扩展
  if( lastdot != NULL )
    *(lastdot + 1) = '\0';
  else
    strcat( new_file, "." );
  strcat( new_file, new_extn );

  return new_file;
}



/*
  prepend_path() 创造文件全路径. 返回路径名
*/
char* prepend_path(const char* path, const char* file)
{
  int name_length = strlen(path) + strlen(file) + 2;
  char* pathname = (char*)malloc( name_length * sizeof(char) );

  sprintf( pathname, name_length, "%s/%s", path, file );

  return pathname;
}



/*
  remove_path() 去除路径名的路径并返回去除路径的文件名 
*/
char* remove_path(const char* pathname)
{
  char* last_slash = strrchr( pathname, '/' );  // last / in pathname
  char* filename;                               // 文件名

  // 如果路径里没有 / , 复制进去
  if( last_slash == NULL )
    {
      filename = (char*)malloc( ( strlen(pathname) + 1 ) * sizeof(char) );
      strcpy( filename, pathname );
    }

  // 反之, 直接在最后加上 /
  else
    {
      filename = (char*)malloc( strlen(last_slash++) * sizeof(char) );
      strcpy( filename, last_slash );
    }

  return filename;
}



/*
  is_image_file() 返回 TRUE 如果文件是图像文件 反之 FALSE
   基于文件扩展名而得
*/
int is_image_file(char* file)
{
  // 找到文件扩展名
  char* lastdot = strrchr(file, '.');

  if( ! lastdot )
    return FALSE;
  
  // 如果含有图像扩展名, 返回 TRUE
  if( ( strcmp(lastdot, ".png") == 0 ) ||
      ( strcmp(lastdot, ".jpg") == 0 ) ||
      ( strcmp(lastdot, ".jpeg") == 0 ) ||
      ( strcmp(lastdot, ".pbm") == 0 ) ||
      ( strcmp(lastdot, ".pgm") == 0 ) ||
      ( strcmp(lastdot, ".ppm") == 0 ) ||
      ( strcmp(lastdot, ".bmp") == 0 ) ||
      ( strcmp(lastdot, ".tif") == 0 ) ||
      ( strcmp(lastdot, ".tiff") == 0 ) )
    {
      return TRUE;
    }
  
  // 反之 FASLE
  return FALSE;
}



/*
  draw_x() 绘制 X 在图像点pt中  The X has 半径 r, 权重 w,
  and 颜色 c.
*/
void draw_x(IplImage* img, CvPoint pt, int r, int w, CvScalar color)
{
  cvLine( img, pt, cvPoint(pt.x+r, pt.y+r), color, w, 8, 0 );
  cvLine( img, pt, cvPoint(pt.x-r, pt.y+r), color, w, 8, 0 );
  cvLine( img, pt, cvPoint(pt.x+r, pt.y-r), color, w, 8, 0 );
  cvLine( img, pt, cvPoint(pt.x-r, pt.y-r), color, w, 8, 0 );
}



/*
  progress() 绘制旋转 (|, /, -, \, |) .
  If done is TRUE, prints "done", otherwise, increments the progress
  thing.
*/
void progress(int done)
{
  static int state = -1;
  
  if( state == -1 )
    fprintf( stdout, "  " );

  if( done )
    {
      fprintf( stdout, "\b\bdone\n");
      state = -1;
    }

  else
    {
      switch( state )
	{
	case 0:
	  fprintf( stdout, "\b\b| ");
	  break;
	  
	case 1:
	  fprintf( stdout, "\b\b/ ");
	  break;
	  
	case 2:
	  fprintf( stdout, "\b\b- ");
	  break;
	  
	default:
	  fprintf( stdout, "\b\b\\ ");
	  break;
	}
      
      fflush(stdout);
      state = (state + 1) % 4;
    }
}



/*
  去除stream一定数目的字符.

  @param stream the stream from which to erase characters
  @param n the number of characters to erase
*/
void erase_from_stream( FILE* stream, int n )
{
  int j;
  for( j = 0; j < n; j++ )
    fprintf( stream, "\b" );
  for( j = 0; j < n; j++ )
    fprintf( stream, " " );
  for( j = 0; j < n; j++ )
    fprintf( stream, "\b" );
}



/*
  数组x2
  
  @param array 数组指针
  @param n 数组元素数量
  @param size 数组大小
  
  @return Returns 返回新的数组元素数量  If no
    memory is available, sets errno to ENOMEM and returns 0.
*/
int array_double( void** array, int n, int size )
{
  void* tmp;

  tmp = realloc( *array, 2 * n * size );
  if( ! tmp )
    {
      if( *array )
	{
	  free( *array );
	  errno = ENOMEM;
	  return 0;
	}
    }
  *array = tmp;
  return n*2;
}
