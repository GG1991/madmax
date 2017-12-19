#ifndef LIST_H
#define LIST_H

typedef int  (*fcmp) (void *, void *);

typedef struct _node_list_t{

    void   *  data;
    struct _node_list_t *  next;

}node_list_t;


typedef struct{

    node_list_t *  head;
    node_list_t *  tail;
    int            sizedata;
    int            sizelist;
    fcmp           cmp;

}list_t;

int list_init(list_t * list, int sizedata, fcmp cmp);
int list_insert_se(list_t * list, void *data);
int list_insertlast(list_t * list, void *data);
int list_delfirst(list_t * list);
int list_del(list_t *list, node_list_t* pNod);
int list_clear(list_t *list);


#endif
