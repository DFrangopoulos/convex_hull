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

//creates a new stack item
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
void push(item_ptr data,stack_item **stack){

    *stack = new_item(data,*stack);
    return;
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

uint8_t stack_transfer(stack_item **src_stack,stack_item **dst_stack){
    item_ptr temp = NULL;
    temp = pop(src_stack);
    if(temp==NULL){
        return 1;
    }
    push(temp,dst_stack);
    return 0;
}


int main (){

    stack_item *stack1 = NULL;
    stack_item *stack2 = NULL;

    uint8_t a = 10;
    uint8_t b = 37;
    uint8_t e = 26;
    uint8_t d = 72;

    push(&a,&stack1);
    push(&b,&stack1);
    push(&d,&stack2);

    stack_transfer(&stack1,&stack2);

    printf("Top S1: %d\n",*(stack1->data));
    printf("Top S2: %d\n",*(stack2->data));

  



    return 0;
}

