#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdint.h>
#include<stdbool.h>
#include<string.h>
#include<math.h>

//typedef the item ptr type in order to make the stack quicker to reuse in different programs
typedef uint8_t *item_ptr;

//stack item type that points to the data and that points to the previous item on the stack
typedef struct stack_item{
    item_ptr data;
    struct stack_item *lower_item;

}stack_item; 

//creates a new item on the stack
stack_item *new_item(item_ptr data, stack_item *lower_item ){

    stack_item *new = (stack_item *)malloc(sizeof(stack_item));
    if(new ==NULL){
        exit(1);
    }
    new->data=data;
    new->lower_item=lower_item;

    return new;
}

//push item to the stack
stack_item *push(item_ptr data,stack_item *stack){

    return new_item(data,stack);
}

//pop item from the stack
item_ptr pop(stack_item **stack){

    item_ptr output = NULL;
    stack_item *lower_item = NULL;

    //check for an empty stack
    if((*stack)!=NULL){
        //extract data for return
        output = (*stack)->data;
        //save address of the lower item (even if it's NULL)
        lower_item = (*stack)->lower_item;
        //delete current item
        free((*stack));
        //save lower item address into the stack pointer (even if it's NULL)
        (*stack) = lower_item;
    }
    return output;
}


int main (){

    stack_item *stack = NULL;
    uint8_t a = 10;
    stack = push(&a,stack);
    uint8_t b = 37;
    stack = push(&b,stack);
    uint8_t e = 26;
    stack = push(&e,stack);
    uint8_t d = 72;
    stack = push(&d,stack);

    printf("Top: %d\n",*(stack->data));

    uint8_t *c = pop(&stack);
    if(c==NULL){
        printf("Nothing to pop!\n");
    }else{
        printf("popped :%d\n",*c);
    }

    c = pop(&stack);
    if(c==NULL){
        printf("Nothing to pop!\n");
    }else{
        printf("popped :%d\n",*c);
    }

    c = pop(&stack);
    if(c==NULL){
        printf("Nothing to pop!\n");
    }else{
        printf("popped :%d\n",*c);
    }

    c = pop(&stack);
    if(c==NULL){
        printf("Nothing to pop!\n");
    }else{
        printf("popped :%d\n",*c);
    }


    uint8_t f = 11;
    stack = push(&f,stack);

    c = pop(&stack);
    if(c==NULL){
        printf("Nothing to pop!\n");
    }else{
        printf("popped :%d\n",*c);
    }

     c = pop(&stack);
    if(c==NULL){
        printf("Nothing to pop!\n");
    }else{
        printf("popped :%d\n",*c);
    }



    return 0;
}

