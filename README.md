# Dependencies
 - gnuplot

# Compilation
 - gcc main.c -lm -o main

# Usage
 - ./main generate <output_file> <point_num>
 - ./main convex <point_file>
 - ./main check <point_file> <point_x> <point_y>
 - ./main intersect <point_file_1> <point_file_2>

# Quickstart
 - First run the generate command to obtain a set of points (or to view the file input format)
 - convex : draws the convex hull for those points with gnuplot
 - check : returns T/F if the point is inside or outside the hull
 - intersect : checks if 2 hulls intersect