#include "myio.h"

#define NRM  "\x1B[0m"
#define RED  "\x1B[31m"
#define GRN  "\x1B[32m"

#define TEST(word, n_test, offset_cmp) { \
  strcpy(str_search, word); \
  myio_file_get_offset_line_start_word(file_name, str_search, &offset); \
  fseek(fm, offset, SEEK_SET); \
  fgets(buf, 256, fm); \
  str_token = strtok(buf, " \n"); \
  printf(NRM "test %d %s" NRM "\n", n_test, (strcmp(str_search, str_token) == 0 && offset == offset_cmp) ? GRN"OK" : RED"FAIL"); \
}

int main(void){

  int offset;
  const char file_name[] = "data/cube_2d.msh";
  char buf[256];
  char *str_token;
  char str_search[128];

  FILE *fm = fopen(file_name, "r"); if(fm == NULL) return 1;

  printf("unity_1\n");

  TEST("$Elements", 1, 4202)
  TEST("$Nodes", 2, 152)
  TEST("$PhysicalNames", 3, 35)
  //  printf("%s %d",str_token, offset);

  fclose(fm);

  return 0;
}
