#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdint.h>
#include<stdbool.h>
#include<string.h>
#include<math.h>


#define MAX_LINE_BUFFER 150
#define _USE_MATH_DEFINES

//Point "Object"
typedef struct point{
    double x;
    double y;
    double angle_to_p0;
}point_s;

//Gnuplot Window Ranges
typedef struct plot_ranges{
    double x_min;
    double x_max;
    double y_min;
    double y_max;
}plot_range_s;

//typedef the item ptr type in order to make the stack quicker to reuse in different programs
typedef point_s *item_ptr;

//stack item type that points to the data and that points to the previous item on the stack
typedef struct stack_item{
    item_ptr data;
    struct stack_item *lower_item;
}stack_item; 

//STACK Prototypes
stack_item *new_item(item_ptr data, stack_item *lower_item);
void push(item_ptr data,stack_item **stack);
item_ptr pop(stack_item **stack);
uint8_t stack_transfer(stack_item **src_stack,stack_item **dst_stack);

//MATH Prototypes
stack_item *polar_order(point_s *all_points,uint32_t point_num,uint32_t low_leftmost_idx);
double get_magnitude(double x, double y);
double get_angle(double x, double y);



void genpoints_to_file(char *fileout, uint32_t point_num);
uint8_t read_point_s_file(char *filename,uint32_t point_num, point_s *all_points);
uint32_t lowest_point_idx(point_s *all_points,uint32_t point_num);
plot_range_s get_range(point_s *all_points,uint32_t point_num);

FILE *open_randsrc();
double random_double(FILE *src);
uint32_t random_int(FILE *src);


//arg1 -> filename
//arg2 -> point_num
//arg3 -> 0/1 regenerate pointfile
int main(int argc , char **argv){

    //Check arg count
    if(argc!=4){
        return 1;
    }

    int point_num =0;
    point_num = atoi(argv[2]);

    if(point_num<1){
        exit(1);
    }

    int regen =0;
    regen = atoi(argv[3]);

    //if regen isn't 0 create new file and use it
    if(regen > 0){  
        genpoints_to_file(argv[1],point_num);
    }

    point_s *all_points = NULL;
    all_points = (point_s *)malloc(sizeof(point_s)*point_num);

    if(all_points==NULL){
        exit(1);
    }

    if (read_point_s_file(argv[1], point_num, all_points)>0){
        exit(1);
    }

    uint32_t low_leftmost_idx = lowest_point_idx(all_points, point_num);
    plot_range_s range = get_range(all_points,point_num);

    stack_item *sorted = polar_order(all_points,point_num,low_leftmost_idx);

    //TEST Pop all//

    item_ptr output = NULL;
    while(true){

        output = pop(&sorted);
        if(output==NULL){
            break;
        }
        printf("[%g,%g] theta=%g\n",output->x,output->y,output->angle_to_p0);

    }


    //Open a pipe to gnuplot
    FILE *gnuplot = popen("gnuplot --persist", "w");
    if(gnuplot==NULL){
        exit(1);
    }


    fprintf(gnuplot,"set key off\n");
    fprintf(gnuplot,"set xrange [%g:%g]\n",range.x_min-1,range.x_max+1);
    fprintf(gnuplot,"set yrange [%g:%g]\n",range.y_min-1,range.y_max+1);
    fprintf(gnuplot, "set style line 1 linecolor rgb 'green' linetype 1 linewidth 0 pointtype 1 pointsize 1.5\n");
    fprintf(gnuplot, "set style line 2 linecolor rgb 'red' linetype 1 linewidth 0 pointtype 1 pointsize 1.5\n");
    fprintf(gnuplot, "plot '-' with points pointtype 1 pointsize 1.5 , '-' with linespoints linestyle 2, '-' with linespoints linestyle 1\n");

    for(uint32_t i=0;i<low_leftmost_idx;i++){
        fprintf(gnuplot,"%g %g\n",all_points[i].x,all_points[i].y);
    }
    fprintf(gnuplot,"e\n");

    fprintf(gnuplot,"%g %g\n",all_points[low_leftmost_idx].x,all_points[low_leftmost_idx].y);
    fprintf(gnuplot,"e\n");

    for(uint32_t i=low_leftmost_idx+1;i<point_num;i++){
        fprintf(gnuplot,"%g %g\n",all_points[i].x,all_points[i].y);
    }
    fprintf(gnuplot,"e\n");

    fflush(gnuplot);

    //fprintf(gnuplot, "set xrange [-10:10]\n");
    //fprintf(gnuplot, "set yrange [-10:10]\n");
    //fprintf(gnuplot, "set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 2 pointtype 7 pointsize 1.5\n");
    //fprintf(gnuplot, "plot '-' with linespoints linestyle 1, '-'\n");

    //for(uint32_t i=0;i<5;i++){
    //    fprintf(gnuplot,"%g %g\n",all_points[i].x,all_points[i].y);
    //}
    //fprintf(gnuplot,"e\n");

    //for(uint32_t i=5;i<15;i++){
    //    fprintf(gnuplot,"%g %g\n",all_points[i].x,all_points[i].y);
    //}
    //fprintf(gnuplot,"e\n");

    //fflush(gnuplot);

    free(all_points);

    return 0;
}

uint8_t read_point_s_file(char *filename,uint32_t point_num, point_s *all_points){

    FILE *input;
    input = fopen(filename,"r");
    if(input==NULL){
        return 1;
    }

    char x_buffer[MAX_LINE_BUFFER];
    char y_buffer[MAX_LINE_BUFFER];

    for(uint32_t i=0; i<point_num;i++){
        for(uint8_t j=0;j<MAX_LINE_BUFFER;j++){
            x_buffer[j]=getc(input);
            if(x_buffer[j]==EOF){exit(1);}
            if(x_buffer[j]==0x20){
                x_buffer[j]=0x00;
                break;
            }
        }

        for(uint8_t k=0;k<MAX_LINE_BUFFER;k++){
            y_buffer[k]=getc(input);
            if(y_buffer[k]==EOF){exit(1);}
            if(y_buffer[k]==0x0a){
                y_buffer[k]=0x00;
                break;
            }
        }

        all_points[i].x = atof(x_buffer);
        all_points[i].y = atof(y_buffer);
    }

    fclose(input);
    return 0;
}

void genpoints_to_file(char *fileout, uint32_t point_num){

    FILE *src, *output;

    src = open_randsrc();
    output = fopen(fileout,"w");

    if(src==NULL || output==NULL){
        exit(1);
    }

    for(uint32_t i=0; i<point_num;i++){
        fprintf(output,"%g %g\n",random_double(src),random_double(src));
    }

    fclose(src);
    fclose(output);

    return;
}

FILE *open_randsrc(){
    return fopen("/dev/urandom","r");
}

double random_double(FILE *src){

    return (((double)random_int(src)/(double)0xFFFFFFFF) * (double)5.0) - (double)2.5;

}

uint32_t random_int(FILE *src){

    if(src==NULL){
        return 0;
    }

    uint32_t rdn;
    fread(&rdn,sizeof(rdn),1,src);
    return rdn;
}

uint32_t lowest_point_idx(point_s *all_points,uint32_t point_num){

    point_s lowest = {all_points[0].x,all_points[0].y};
    uint32_t index = 0;

    for(uint32_t i=1;i<point_num;i++){
        //record index of lowest point
        if(lowest.y>all_points[i].y){
            lowest.y = all_points[i].y;
            lowest.x = all_points[i].x;
            index = i;
        //record index of lowest and leftmost point
        }else if((lowest.y==all_points[i].y)&&(lowest.x>all_points[i].x)){
            lowest.x = all_points[i].x;
            index = i; 
        }
    }
    return index;
}

plot_range_s get_range(point_s *all_points,uint32_t point_num){

    plot_range_s range = {all_points[0].x,all_points[0].x,all_points[0].y,all_points[0].y};

    for(uint32_t i=1;i<point_num;i++){

        if(range.x_min > all_points[i].x){
            range.x_min = all_points[i].x;
        }
        if(range.x_max < all_points[i].x){
            range.x_max = all_points[i].x;
        }
        if(range.y_min > all_points[i].y){
            range.y_min = all_points[i].y;
        }
        if(range.y_max < all_points[i].y){
            range.y_max = all_points[i].y;
        }
    }

    return range;
}

//----START STACK FUNCS
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

//transfer item from one stack to another
uint8_t stack_transfer(stack_item **src_stack,stack_item **dst_stack){
    item_ptr temp = NULL;
    temp = pop(src_stack);
    if(temp==NULL){
        return 1;
    }
    push(temp,dst_stack);
    return 0;
}

//----END STACK FUNCS

//----START MATH FUNCS
double get_angle(double x, double y){
    return atan2(y,x);
}

double get_magnitude(double x, double y){
    return sqrt(pow(x,2)+pow(y,2));
}

//Create two stacks for sorting
//Move items from one stack to the  
stack_item *polar_order(point_s *all_points,uint32_t point_num,uint32_t low_leftmost_idx){

    //Create an initial stack with all points but p0

    stack_item *initial_stack = NULL;
    //Push the rest
    for(uint32_t i=0;i<point_num;i++){

        if(i!=low_leftmost_idx){
            //printf("pushed:[%g,%g]\n",all_points[i].x,all_points[i].y);
            push(&all_points[i],&initial_stack);
        }
    }

    //Create two empty stacks for sorting
    stack_item *stack1 = NULL;
    stack_item *stack2 = NULL;

    //point being evaluated
    item_ptr eval = NULL;
    while(true){
        //Pop_item from the initial_stack
        eval = pop(&initial_stack);
        double eval_mag,top_stack_mag;
        if(eval!=NULL){
            //Get the angle wtr to p0
            eval->angle_to_p0 = get_angle(eval->x-all_points[low_leftmost_idx].x,eval->y-all_points[low_leftmost_idx].y);
            //Compare to angle at the top of stack1
            while(true){
                //if the stack is not empty
                if(stack1!=NULL){
                    //if eval's angle is larger than the top of stack1 
                    //pop the top of stack1 and push it onto stack2
                    if(eval->angle_to_p0 > stack1->data->angle_to_p0){
                        if(stack_transfer(&stack1,&stack2)){
                            exit(1);
                        }
                    //if it's equal compare magnitudes (collinear vectors)
                    }else if(eval->angle_to_p0 == stack1->data->angle_to_p0){

                        eval_mag = get_magnitude(eval->x-all_points[low_leftmost_idx].x,eval->y-all_points[low_leftmost_idx].y);
                        top_stack_mag = get_magnitude(stack1->data->x-all_points[low_leftmost_idx].x,stack1->data->y-all_points[low_leftmost_idx].y);

                        //printf("[%g,%g] mag:%g angle:%g\n",eval->x,eval->y,eval_mag,eval->angle_to_p0);
                        //printf("[%g,%g] mag:%g angle:%g\n",stack1->data->x,stack1->data->y,top_stack_mag,stack1->data->angle_to_p0);
                        //less than 90 degress (Q1) largest magnitude goes onto the stack first
                        if((*eval).angle_to_p0 < M_PI_2){
                            //if eval's angle is smaller 
                            if(eval_mag < top_stack_mag){
                                //push eval onto stack1
                                push(eval,&stack1);
                            }
                            else{
                                //push top of stack1 onto stack2
                                if(stack_transfer(&stack1,&stack2)){
                                    exit(1);
                                }
                                //push eval onto stack 1 and add everything from stack 2 on top of it
                                push(eval,&stack1);
                                //push all items from stack2 onto stack1
                                while(true){
                                    if(stack_transfer(&stack2,&stack1)){
                                        //stack2 is empty
                                        break;
                                    }
                                }
                                //move to next item on initial_stack
                                break;
                            }
                        }else{
                            //more than 90 degress (Q2) smallest magnitude goes onto the stack first
                            if(eval_mag < top_stack_mag){
                                //push top of stack1 onto stack2
                                if(stack_transfer(&stack1,&stack2)){
                                    exit(1);
                                }
                                //push eval onto stack 1 and add everything from stack 2 on top of it
                                push(eval,&stack1);
                                //push all items from stack2 onto stack1
                                while(true){
                                    if(stack_transfer(&stack2,&stack1)){
                                        //stack2 is empty
                                        break;
                                    }
                                }
                                //move to next item on initial_stack
                                break;
                            }else{
                                //push eval onto stack1
                                push(eval,&stack1);
                            }
                        }


                    //push eval onto stack 1 and add everything from stack 2 on top of it
                    }else{
                        //if eval's angle is smaller push it onto stack1
                        push(eval,&stack1);
                        //push all items from stack2 onto stack1
                        while(true){
                            if(stack_transfer(&stack2,&stack1)){
                                //stack2 is empty
                                break;
                            }
                        }
                        //move to next item on initial_stack
                        break;
                    }
                }else{
                    //just push it onto stack 1 it's the first element
                    push(eval,&stack1);
                    //move to next item on initial_stack
                    break;
                }
            }
        }else{
            //all items from initial stack have been processed
            break;
        }
    }

    //stack1 contained all points except p0 sorted by smallest angle relative
    //to p0 with collinear points taken into account

    return stack1;
}


