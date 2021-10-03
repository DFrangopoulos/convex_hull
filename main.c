#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdint.h>
#include<stdbool.h>
#include<string.h>
#include<math.h>

#define MAX_LINE_BUFFER 150

//Point "Object"
typedef struct point{
    double x;
    double y;
    double angle_wrt_p0;
}point_s;

//typedef the item ptr type in order to make the stack quicker to reuse in different programs
typedef point_s *item_ptr;

//Gnuplot Window Ranges
typedef struct plot_ranges{
    double x_min;
    double x_max;
    double y_min;
    double y_max;
}plot_range_s;

//stack item type that points to the data and that points to the previous item on the stack
typedef struct stack_item{
    item_ptr data;
    struct stack_item *lower_item;
}stack_item; 

FILE *open_randsrc();
double random_double(FILE *src);
uint32_t random_int(FILE *src);
void genpoints_to_file(char *fileout, uint32_t point_num);
uint8_t read_point_s_file(char *filename,uint32_t point_num, point_s *all_points);

uint32_t lowest_point_idx(point_s *all_points,uint32_t point_num);
plot_range_s get_range(point_s *all_points,uint32_t point_num);

stack_item *new_item(item_ptr data, stack_item *lower_item);
void push(item_ptr data,stack_item **stack);
item_ptr pop(stack_item **stack);
uint8_t stack_transfer(stack_item **src_stack,stack_item **dst_stack);
void dump_stack(stack_item **stack);

double get_angle(double x1, double y1, double p0x, double p0y);
double get_magnitude(double x1, double y1, double p0x, double p0y);
stack_item *polar_order(point_s *all_points,uint32_t point_num,uint32_t low_leftmost_idx);


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

    dump_stack(&sorted);

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

//open /dev/urandom as a source
FILE *open_randsrc(){
    return fopen("/dev/urandom","r");
}
//compute a random double
double random_double(FILE *src){
    return (((double)random_int(src)/(double)0xFFFFFFFF) * (double)5.0) - (double)2.5;
}
//return a random int using a randomness source
uint32_t random_int(FILE *src){
    if(src==NULL)
        return 0;
    uint32_t rdn;
    fread(&rdn,sizeof(rdn),1,src);
    return rdn;
}
//Generate random points and write them to a file
void genpoints_to_file(char *fileout, uint32_t point_num){
    FILE *src, *output;
    src = open_randsrc();
    output = fopen(fileout,"w");
    if(src==NULL || output==NULL)
        exit(1);

    for(uint32_t i=0; i<point_num;i++){
        fprintf(output,"%g %g\n",random_double(src),random_double(src));
    }
    fclose(src);
    fclose(output);
    return;
}
//Read a file containing points
uint8_t read_point_s_file(char *filename,uint32_t point_num, point_s *all_points){
    FILE *input;
    input = fopen(filename,"r");
    if(input==NULL)
        return 1;
    char x_buffer[MAX_LINE_BUFFER];
    char y_buffer[MAX_LINE_BUFFER];
    for(uint32_t i=0; i<point_num;i++){
        for(uint8_t j=0;j<MAX_LINE_BUFFER;j++){
            x_buffer[j]=getc(input);
            if(x_buffer[j]==EOF)
                exit(1);
            if(x_buffer[j]==0x20){
                x_buffer[j]=0x00;
                break;
            }
        }
        for(uint8_t k=0;k<MAX_LINE_BUFFER;k++){
            y_buffer[k]=getc(input);
            if(y_buffer[k]==EOF)
                exit(1);
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

//find the index of the lowest and leftmost point aka p0
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
        } else if((lowest.y==all_points[i].y)&&(lowest.x>all_points[i].x)){
            lowest.x = all_points[i].x;
            index = i; 
        }
    }
    return index;
}
//find the maximum ranges for the plot
plot_range_s get_range(point_s *all_points,uint32_t point_num){
    plot_range_s range = {all_points[0].x,all_points[0].x,all_points[0].y,all_points[0].y};
    for(uint32_t i=1;i<point_num;i++){
        if(range.x_min > all_points[i].x)
            range.x_min = all_points[i].x;
        if(range.x_max < all_points[i].x)
            range.x_max = all_points[i].x;
        if(range.y_min > all_points[i].y)
            range.y_min = all_points[i].y;
        if(range.y_max < all_points[i].y)
            range.y_max = all_points[i].y;
    }
    return range;
}

//creates a new stack item
stack_item *new_item(item_ptr data, stack_item *lower_item){
    stack_item *new = (stack_item *)malloc(sizeof(stack_item));
    if(new ==NULL)
        exit(1);
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
    if(temp==NULL)
        return 1;
    push(temp,dst_stack);
    return 0;
}
//pop and print all elements fro stack
void dump_stack(stack_item **stack){
    item_ptr output = NULL;
    while(true){
        output = pop(stack);
        if(output==NULL)
            break;
        printf("[%g,%g] angle=%g\n",output->x,output->y,output->angle_wrt_p0);
    }
    return;
}

//compute the angle at which the point is wrt to p0
double get_angle(double x1, double y1, double p0x, double p0y){
    return atan2(y1-p0y,x1-p0x);
}
//compute the magnitude of the p0 - point vector
double get_magnitude(double x1, double y1, double p0x, double p0y){
    return sqrt(pow((x1-p0x),2)+pow((y1-p0y),2));
}
//Order all points wrt p0 by angle and return a stack (top of the stack should be the point with the smallest angle)
stack_item *polar_order(point_s *all_points,uint32_t point_num,uint32_t low_leftmost_idx){

    //Define p0
    point_s p0 = all_points[low_leftmost_idx];

    //Create 3 stacks
    stack_item *initial = NULL, *s1 = NULL, *s2 = NULL;

    //Push all points from all_points (- p0) onto initial;
    for (uint32_t i = 0; i < point_num; i++)
    {
        if(i!=low_leftmost_idx){
            push(&(all_points[i]),&initial);
        }
    }

    item_ptr evald = NULL;
    while(true){

        evald = pop(&initial);
        if(evald==NULL)
            break;
        evald->angle_wrt_p0 = get_angle(evald->x,evald->y,p0.x,p0.y);
        double m_evald, m_top;
        
        while(true){
            if(s1==NULL){
                //Case 1: S1 is empty
                push(evald,&s1);
                //push s2 onto s1
                while(!stack_transfer(&s2,&s1));
                break;
            } else if(evald->angle_wrt_p0 > s1->data->angle_wrt_p0){
                //Case 2: the angle of evald is larger than the angle of the top element of s1
                //move top of s1 to s2
                stack_transfer(&s1,&s2);
            } else if(evald->angle_wrt_p0 < s1->data->angle_wrt_p0){
                //Case 3: the angle of evald is smaller than the angle of the top element of s1
                //push eval onto s1
                push(evald,&s1);
                //push s2 onto s1
                while(!stack_transfer(&s2,&s1));
                break;
            } else {
                //Case 4: the angles of evald and the top element of s1 are equal
                //collinear case
                //Compute magnitudes
                m_evald = get_magnitude(evald->x,evald->y,p0.x,p0.y);
                m_top = get_magnitude(s1->data->x,s1->data->y,p0.x,p0.y);

                if(evald->angle_wrt_p0<M_PI_2){
                    //if the angle is strictly less than 90degrees
                    //the largest magnitude goes onto the stack first
                    if(m_evald > m_top){
                        //move top of s1 to s2
                        stack_transfer(&s1,&s2);
                    } else {
                        //push eval onto s1
                        push(evald,&s1);
                        //push s2 onto s1
                        while(!stack_transfer(&s2,&s1));
                        break;
                    }
                } else {
                    //if the angle is 90 degrees or more
                    //the smallest magnitude goes onto the stack first
                    if(m_evald < m_top){
                        //move top of s1 to s2
                        stack_transfer(&s1,&s2);
                    } else {
                        //push eval onto s1
                        push(evald,&s1);
                        //push s2 onto s1
                        while(!stack_transfer(&s2,&s1));
                        break;
                    }
                }
            }
        }
    }
    return s1;    
}
