/*
  utils.c

  ��c�ļ������˸�ʽ����ʵ�õĺ���.
*/

#include "utils.h"
#include "defs.h"


/***************************** �������� ******************************/

/* ����32λ����ͼ���һ������ֵ */
// float pixval32f(IplImage* img, int r, int c)
// {
//   return ( (float*)(img->imageData + img->widthStep*r) )[c];
// }


/* ����64λ����ͼ������ֵ  */
// void setpix32f(IplImage* img, int r, int c, float val)
// {
//   ( (float*)(img->imageData + img->widthStep*r) )[c] = val;
// }

/*************************** �������� ****************************/

/*
  fatal_error() ��ӡ��ʽ������Ϣ ���˳�
*/
void fatal_error(char* format, ...)
{
  va_list ap;  // ��ʽ�����б�
  
  fprintf( stderr, "Error: ");

  // assemble argument list and print format using it
  va_start( ap, format );
  vfprintf( stderr, format, ap );
  va_end(ap);
  printf("\n");
  exit(1);
}



/*
  replace_extension() ���µ���չ���滻�ļ���չ�� �����ؽ��.  
*/
char* replace_extension(const char* file, const char* new_extn)
{
  char* new_file = malloc( (strlen(file)+strlen(new_extn)+2) * sizeof(char) );
  char* lastdot;   // ��չλ��

  // �����ļ������ҵ���չλ��
  strcpy( new_file, file );
  lastdot = strrchr( new_file, '.' );

  // ����ļ�������չ��, ȥ����,�����µ���չ
  if( lastdot != NULL )
    *(lastdot + 1) = '\0';
  else
    strcat( new_file, "." );
  strcat( new_file, new_extn );

  return new_file;
}



/*
  prepend_path() �����ļ�ȫ·��. ����·����
*/
char* prepend_path(const char* path, const char* file)
{
  int name_length = strlen(path) + strlen(file) + 2;
  char* pathname = (char*)malloc( name_length * sizeof(char) );

  sprintf( pathname, name_length, "%s/%s", path, file );

  return pathname;
}



/*
  remove_path() ȥ��·������·��������ȥ��·�����ļ��� 
*/
char* remove_path(const char* pathname)
{
  char* last_slash = strrchr( pathname, '/' );  // last / in pathname
  char* filename;                               // �ļ���

  // ���·����û�� / , ���ƽ�ȥ
  if( last_slash == NULL )
    {
      filename = (char*)malloc( ( strlen(pathname) + 1 ) * sizeof(char) );
      strcpy( filename, pathname );
    }

  // ��֮, ֱ���������� /
  else
    {
      filename = (char*)malloc( strlen(last_slash++) * sizeof(char) );
      strcpy( filename, last_slash );
    }

  return filename;
}



/*
  is_image_file() ���� TRUE ����ļ���ͼ���ļ� ��֮ FALSE
   �����ļ���չ������
*/
int is_image_file(char* file)
{
  // �ҵ��ļ���չ��
  char* lastdot = strrchr(file, '.');

  if( ! lastdot )
    return FALSE;
  
  // �������ͼ����չ��, ���� TRUE
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
  
  // ��֮ FASLE
  return FALSE;
}



/*
  draw_x() ���� X ��ͼ���pt��  The X has �뾶 r, Ȩ�� w,
  and ��ɫ c.
*/
void draw_x(IplImage* img, CvPoint pt, int r, int w, CvScalar color)
{
  cvLine( img, pt, cvPoint(pt.x+r, pt.y+r), color, w, 8, 0 );
  cvLine( img, pt, cvPoint(pt.x-r, pt.y+r), color, w, 8, 0 );
  cvLine( img, pt, cvPoint(pt.x+r, pt.y-r), color, w, 8, 0 );
  cvLine( img, pt, cvPoint(pt.x-r, pt.y-r), color, w, 8, 0 );
}



/*
  progress() ������ת (|, /, -, \, |) .
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
  ȥ��streamһ����Ŀ���ַ�.

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
  ����x2
  
  @param array ����ָ��
  @param n ����Ԫ������
  @param size �����С
  
  @return Returns �����µ�����Ԫ������  If no
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
