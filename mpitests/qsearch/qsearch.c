#include <stdio.h>

int search_position_logn(int *array, int size, int val, int *pos);

int main(int argc, char **argv)
{
  int array[] = { 1, 3, 4, 6, 8, 9, 10, 11, 12, 20 };
  int n=10, p, val = atoi(argv[1]);

  search_position_logn(array, n, val, &p);

  if(p>=0){
    printf("the value %d was found on position %d\n",val,p); 
  }
  else{
    printf("the value %d was not found\n",val); 
  }

  return 1;
}

int search_position_logn(int *array, int size, int val, int *pos)
{
  /* Returns: 
   * a) the position <pos> of <val> inside <array> (size <size>)
   * b) <pos> = -1 if <val> does not exist 
   *
   * Note: the array should be sorted
   *
   */

  int  left = 0, right = size-1, middle;

  while(left <= right){

    middle = (right + left)/2; 
    if(array[middle] == val){
      *pos = middle;
      return 0;
    }
    if(array[middle] < val){
      left = middle + 1;
    }
    else{
      right = middle - 1;
    }
  }
  *pos = -1;
  return 0;
}
