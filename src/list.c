/* list_t operations
 * 
 * 
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include "malloc.h"
#include "string.h"
#include "list.h"


int list_init(list_t * list, int sizedata, fcmp cmp){

    list->head     = NULL;
    list->tail     = NULL;
    list->sizedata = sizedata;
    list->sizelist = 0;
    list->cmp      = cmp;

    return 0;

}

int list_insert_se(list_t * list, void *data){

    /* Insert sort exclusive an element
     * 
     * 
     * return 0  if success
     *        1  if repeated
     *       -1  error
     * 
     */

    node_list_t * node, *onode;
    void * aux;

    if(!list)
        return 1;
    node = (node_list_t *)malloc(sizeof(node_list_t));
    if(!node)
        return -1;
    if (list->sizedata){
        node->data = (void*)malloc(list->sizedata);
        if(!node->data)
            return 1;
        memcpy(node->data,data,list->sizedata);
    }else{
        node->data = data;
    }
    node->next=NULL;
    if(!list->cmp)
        return -1;
    if(list->sizelist==0){
        list->sizelist ++;
        list->head = list->tail = node;
        return 0;
    }
    onode = list->head;
    while((*list->cmp)(onode->data,data)<0 && onode->next!=NULL){
        onode = onode->next;
    }
    if((*list->cmp)(onode->data,data)==0){
        free(node);
        return 1;
    }
    list->sizelist ++;
    node->next = onode->next;
    onode->next = node;
    if((*list->cmp)(onode->data,data)>0){
        aux = onode->data;
        onode->data = node->data;
        node->data = aux;
    }    
    if(onode==list->tail)
        list->tail=node;
    return 0;

}

int list_insertlast(list_t * list, void *data){

    /* Insert sort exclusive an element
     * 
     * return 0  if success
     *        1  error 
     */

    node_list_t * node;

    if(!list)
        return 1;
    node = (node_list_t *)malloc(sizeof(node_list_t));
    if(!node)
        return 1;
    if(list->sizedata){
        node->data = (void*)malloc(list->sizedata);
        if(!node->data)
            return 1;
        memcpy(node->data,data,list->sizedata);
    }else{
        node->data = data;
    }
    node->next=NULL;
    if(list->sizelist==0){
        list->sizelist ++;
        list->head = list->tail = node;
        return 0;
    }
    list->sizelist ++;
    list->tail->next = node;
    list->tail = node;    
    return 0;

}

int list_delfirst(list_t * list){

    /* Deletes the first element in the list
     * 
     * return 0 sucess.
     *        1 if list=NULL or if there is no element to delete.
     */

    node_list_t *aux;
    if(!list)
        return 1;
    if(list->sizelist==0) 
        return 1;
    if(list->sizelist==1){
        free(list->head);
        list->head=NULL;
        list->tail=NULL;
        list->sizelist--;
        return 0;
    }
    aux = list->head->next;
    free(list->head);
    list->head=aux;    
    list->sizelist--;
    return 0;
}

int list_del(list_t *list, node_list_t* pNod){

    /* Removes the element if 
     * exists 
     */

    node_list_t *pNodA;
    if(!list)
        return 1;
    if(!pNod)
        return 1;
    if(!list->sizelist)
        return 0;
    pNodA=list->head;
    if(pNodA==pNod){
        list->sizelist--;
        list->head=pNod->next;
        if(list->sizelist==0){
            list->tail=NULL;
        }
        free(pNod);
        return 0;
    }
    while(pNodA->next){
        if(pNodA->next==pNod)
            break;
        pNodA=pNodA->next;
    }
    if(!pNodA->next)
        return 0;
    list->sizelist--;
    pNodA->next=pNod->next;
    if(!pNod->next)
        list->tail=pNodA;
    free(pNod);
    return 0;
}

int list_free(list_t *list){

    /* Frees the memory allocated by the list 
    */

    if(!list)
        return 0;
    while(list->sizelist){
        list_delfirst(list);
    }
    free(list);
    return 0;
}


