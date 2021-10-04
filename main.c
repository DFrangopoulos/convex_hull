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
void duplicate_stack(stack_item **in_stack, stack_item **copy_stack);
void rev_duplicate_stack(stack_item **in_stack, stack_item **rev_copy_stack);

double get_angle(double x1, double y1, double p0x, double p0y);
double get_magnitude(double x1, double y1, double p0x, double p0y);
double cp2d(item_ptr a, item_ptr b, item_ptr c);

stack_item *polar_order(point_s *all_points,uint32_t point_num,uint32_t low_leftmost_idx);
bool cwt(item_ptr p, item_ptr c,item_ptr n);
void graham_scan(stack_item **sorted_stack, stack_item **hull_stack, stack_item **inner_stack);
bool is_in_hull(stack_item** hull,point_s test_point);


FILE *setup_gnuplot(plot_range_s range);
void plot_points(stack_item **point_stack, FILE *gnuplot);
void close_gnuplot(FILE *gnuplot);

//arg1 -> filename
//arg2 -> point_num
//arg3 -> 0/1 regenerate pointfile
int main(int argc , char **argv){

    //Check arg count
    if(argc!=4)
        return 1;

    //Get point count
    int point_num =0;
    point_num = atoi(argv[2]);
    if(point_num<1)
        exit(1);

    //Check if the point file needs to be regenerated
    int regen = 0;
    regen = atoi(argv[3]);
    //if regen isn't 0 create new file and use it
    if(regen > 0)
        genpoints_to_file(argv[1],point_num);

    //Load points
    point_s *all_points = NULL;
    all_points = (point_s *)malloc(sizeof(point_s)*point_num);

    if(all_points==NULL)
        exit(1);

    if(read_point_s_file(argv[1], point_num, all_points)>0)
        exit(1);

    //Find p0 index
    uint32_t low_leftmost_idx = lowest_point_idx(all_points, point_num);
    //set angle to 0
    all_points[low_leftmost_idx].angle_wrt_p0=0;
    //calculate range
    plot_range_s range = get_range(all_points,point_num);

    //Sort points into a stack with smallest angle wrt p0 on top
    stack_item *sorted = polar_order(all_points,point_num,low_leftmost_idx);
    //Put p0 on top
    push(&(all_points[low_leftmost_idx]),&sorted);

    //Graham Scan
    stack_item *hull_stack = NULL, *inner_stack = NULL;
    graham_scan(&sorted, &hull_stack, &inner_stack);

    //Reverse Duplicate hull_stack
    stack_item *rev_hull_stack = NULL;
    rev_duplicate_stack(&hull_stack, &rev_hull_stack);

    point_s test_point = {-2.2 , -1.5 , 0};

    if(is_in_hull(&rev_hull_stack, test_point))
        printf("In the hull\n");
    else
        printf("Not in the hull\n");

    //Find if point is in convex hull


    //printf("RDupHull\n");
    //dump_stack(&rev_hull_stack);

    //printf("Inner\n");
    //dump_stack(&inner_stack);

    //push p0 on onto hull_stack for gnuplot to close the circle
    push(&(all_points[low_leftmost_idx]),&hull_stack);

    FILE *gnuplot = setup_gnuplot(range);
    if(gnuplot==NULL)
        exit(1);

    plot_points(&hull_stack,gnuplot);
    plot_points(&inner_stack,gnuplot);
    close_gnuplot(gnuplot);
    
    
    //fprintf(gnuplot, "set style line 1 linecolor rgb 'green' linetype 1 linewidth 0 pointtype 1 pointsize 1.5\n");
    //fprintf(gnuplot, "set style line 2 linecolor rgb 'red' linetype 1 linewidth 0 pointtype 1 pointsize 1.5\n");
    //fprintf(gnuplot, "plot '-' with points pointtype 1 pointsize 1.5 , '-' with linespoints linestyle 2, '-' with linespoints linestyle 1\n");

    free(all_points);

    return 0;
}

//Open gnuplot pipe and setup
FILE *setup_gnuplot(plot_range_s range){
    //Open a pipe to gnuplot
    FILE *gnuplot = popen("gnuplot --persist", "w");
    if(gnuplot==NULL)
        return NULL;
    //Setup plot with ranges
    fprintf(gnuplot,"set key off\n");
    fprintf(gnuplot,"set xrange [%g:%g]\n",range.x_min-1,range.x_max+1);
    fprintf(gnuplot,"set yrange [%g:%g]\n",range.y_min-1,range.y_max+1);
    //Create hull linestyle
    fprintf(gnuplot, "set style line 1 linecolor rgb 'green' linetype 1 linewidth 0 pointtype 1 pointsize 1.5\n");
    //Prepare plotting first plot hull then the inner points
    fprintf(gnuplot, "plot  '-' with linespoints linestyle 1, '-' with points pointtype 1 pointsize 1.5\n");
    fflush(gnuplot);
    return gnuplot;
}
//Plot points
void plot_points(stack_item **point_stack, FILE *gnuplot){
    item_ptr output = NULL;
    while(true){
        output = pop(point_stack);
        if(output==NULL)
            break;
        fprintf(gnuplot,"%g %g\n",output->x,output->y);
    }
    fprintf(gnuplot,"e\n");
    fflush(gnuplot);
    return;
}
//Close pipe
void close_gnuplot(FILE *gnuplot){
    fflush(gnuplot);
    pclose(gnuplot);
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
//duplicate input stack into copy stack
void duplicate_stack(stack_item **in_stack, stack_item **copy_stack){
    stack_item *tmp = NULL;
    //flip all elements
    while(!stack_transfer(in_stack,&tmp));
    //pushback to in_stack and copy
    item_ptr tmp2 = NULL;
    while(true){
        tmp2 = pop(&tmp);
        if(tmp2==NULL)
            break;
        push(tmp2,in_stack);
        push(tmp2,copy_stack);
    }
    return;
}
//reverse duplicate input stack into rev copy stack
void rev_duplicate_stack(stack_item **in_stack, stack_item **rev_copy_stack){
    stack_item *tmp = NULL;
    item_ptr tmp2 = NULL;
    //flip all elements and copy to two stacks
    while(true){
        tmp2 = pop(in_stack);
        if(tmp2==NULL)
            break;
        push(tmp2,rev_copy_stack);
        push(tmp2,&tmp);
    }
    //transfer tmp back to input
    while(!stack_transfer(&tmp,in_stack));
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


//Check for a clockwise turn using the 2D "cross product"
bool cwt(item_ptr p, item_ptr c,item_ptr n){
    //if the 2D "cross-product" is less than 0 then the turn is clockwise
    //also treat collinear vectors as a clockwise turn
    if( ( (c->x - p->x) * (n->y - c->y) - (c->y - p->y) * (n->x - c->x) ) <=0 )
        return true;
    return false;
}
//Graham Scan
void graham_scan(stack_item **sorted_stack, stack_item **hull_stack, stack_item **inner_stack){

    //transfer p0 to hull_stack
    stack_transfer(sorted_stack,hull_stack);

    //add first point (after p0) to the hull stack to prime the loop
    stack_transfer(sorted_stack,hull_stack);
    item_ptr previous, current, next;

    while(true){
        //get next point
        next = pop(sorted_stack);
        if(next == NULL)
            break;
        
        if(*hull_stack!=NULL){
            if((*hull_stack)->lower_item!=NULL){
                current = (*hull_stack)->data;
                previous = (*hull_stack)->lower_item->data;
            }
        }

        while(true){
            //if a cwt is made discard current point
            if(cwt(previous,current,next)){
                stack_transfer(hull_stack,inner_stack);
                //update points
                if(*hull_stack!=NULL){
                    if((*hull_stack)->lower_item!=NULL){
                        current = (*hull_stack)->data;
                        previous = (*hull_stack)->lower_item->data;
                    } else {
                        break;
                    }
                } else {
                    break;
                }
            } else {
                push(next,hull_stack);
                break;
            }
        }
        
    }

    return;
}

//2D cross-product
double cp2d(item_ptr a, item_ptr b, item_ptr c){
    //||BA x BC|| = (ax-bx)*(cy-by)-(ay-by)*(cx-bx)
    return ((a->x - b->x) * (c->y - b->y) - (a->y - b->y) * (c->x - b->x));
}

//find convex hull sector
bool is_in_hull(stack_item** hull,point_s test_point){
    if(hull==NULL)
        exit(1);
    point_s p0 = *pop(hull);
    //find the sector in which the point lies by
    //doing the 2D "cross product" of p0pi x p0test
    //which needs to be positive and by maximizing i
    item_ptr previous = NULL, current = NULL;
    while(true){
        current = pop(hull);
        if(current==NULL)
            break;
        //do the 2D cross product p0pi x p0test
        if(cp2d(current,&p0,&test_point)< 0){
            break;
        }
        previous = current;
    }
    if(current==NULL)
        return false;

    //Check if the point is within the portion of the sector that is inside the convex hull
    //2D "cross product" of pipi+1 x piptest
    //which needs to be positive
    if(cp2d(current,previous,&test_point)>= 0)
        return true;
    else
       return false;
}