
/*
   Giant util functions
 */

int cmp_int(void *a, void *b)
{
    if ( *(int*)a > *(int*)b){
	return 1;
    }else if(*(int*)a == *(int*)b){
	return 0;
    }else{
	return -1;
    }
}
