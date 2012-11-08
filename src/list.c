#include<stdio.h>
#include "list.h"

void push(list *s, int val) {
	if (full(s)) {
		fprintf(stderr, "Stack overflows!");
		exit(1);
	}
	s->v[++s->top] = val;
}

int pop(list *s) {
	if (empty(s)) {
		fprintf(stderr, "Stack empty!");
		return -1;
	}
	return (s->v[s->top--]);
}

int pop2(list *s) {
	if (empty(s)) {
		fprintf(stderr, "Stack empty!");
		return -1;
	}
	return (s->v[s->cursor--]);
}

void init(list *s) {
	s->top = 0;
	s->cursor = 0;
}

int full(list *s) {
	return (s->top >= STACK_SIZE - 1);
}

int empty(list *s) {
	return (s->top == 0);
}

int get_size(list *s) {
	return s->top;
}

int find(list *s, int value) {
	int i = s->top;
	while (i) {
		if (value == s->v[i])
			return i;
		i--;
	}
	return 0;
}

void clear(list *s) {
	while (s->top--)
		s->v[s->top + 1] = 0;
	s->cursor = 0;
}

void p_stack(list *s) {
	int i;
	s->cursor = s->top;
	if (s->top == 0)
		printf("Stack is empty.\n");
	else {
		printf("Stack contents: ");
		while (!empty(s)) {
			printf("%d ", pop2(s));
		}
		printf("\n");
	}
}
