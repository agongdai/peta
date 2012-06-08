/* mystack.h -- Stack declaration and function prototypes:  */
#ifndef STACK_H_
#define STACK_H_

#define STACK_SIZE 1024
#define NOT_FOUND = -999;

typedef struct
{
    int v[STACK_SIZE];
    int top;
    int cursor;
} list;

void push(list *s, int val);
int pop(list *s);
int pop2(list *s);
void init(list *s);
int full(list *s);
int empty (list *s);
void p_stack(list *s);
void clear(list *s);
int find(list *s, int value);
int get_size(list *s);

#endif
