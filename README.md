# Computational-Geometry
Course work - IIT Kharagpur

chan.cpp - Chan's algorithm for computing convex hull.
Input: A set of n points in the plane.
Output: The ordered list of hull points. Also outputs the intermediate result after Graham's Scan technique.

intersec.cpp - line segment intersection algorithm.
Input: A set of n line segments (taken as pair of points).
Output: The intersection point. Also outputs the convex hull formed out of the intersection points.

maximize.cpp - Randomized Incremental 2D Linear Programming algorithm. Clearly describe the data structure you have used. Upload a  write-up supporting your data structure.
Input: A set of N  linear inequalities on the plane and a non-zero objective function C.
Output: A feasible point P that maximizes the dot product between objective function and the point (C.P). C and P both are vectors.

Description of the Incremental 2D Linear Programming algorithm
--------------------------------------------------------------

Half-planes on line. 
Finally our aim is to find the solution of 2D optimization problem in expected linear time.
We need O(n) time incase the solution of i inequalities is different from the first i-1.
If it is same we just require O(1) time. This leads to O(n) time complexity in expected case using Randomized Permutation of Half Plane.

The main datastructure we are using are Point, Halfplane, Projection,  which is shown below:


---------------------------------------------------------------------------------------------------------

	struct Point
	{
		double x, y;

		Point() {}
		Point(double _x, double _y) : x(_x), y(_y) {}

		friend ostream& operator<<( ostream& output, const Point& p )
		{
			char temp[500];
			sprintf(temp,"(%0.2lf,%0.2lf)",p.x,p.y);
			output << temp;
			return output;
		}   
	};

a point on a plane is represented by its x and y co-ordinates

-----------------------------------------------------------------------------------------------------------


	struct Halfplane
	{
		// Half palne is represented by inequality
		//			ax + by <= c
		double a,b,c;
		double m;	// slope
	
		Halfplane(double _a, double _b, double _c) : a(_a), b(_b), c(_c) {
			if (_b != 0) {
				m = -a/b;
			}
		}
	
		friend std::ostream& operator<<(std::ostream& os, const Halfplane& obj)
		{
	   	char temp[100];
	   	sprintf(temp,"%0.2lfx + %0.2lfy <= %0.2lf",obj.a,obj.b,obj.c);
	   	os << temp;
	   	return os;
		}
	
	};

A Half plane is represented by its inequality 'ax + by <= c' <br>
m gives the slope of line ax + by = c


-----------------------------------------------------------------------------------------------------------


	const int LINE = 1;
	const int RAY = 2;
	const int LINESEGMENT = 3;
	const int EMPTY = 4;
	const int UP = 11;
	const int DOWN = 12;


	struct Projection
	{
		Point p,q;	// q should not be accessed if it is a RAY
					// p  should be below q if it is a LINESEGMENT
		double a,b,c;	// ax + by <= c
		double m;	// slope
		int type;	// type can be LINE, RAY, LINESEGMENT, EMPTY
		int dir;	// direction if the projection is a ray
	
		void setabc(Projection proj) {
			a = proj.a;
			b = proj.b;
			c = proj.c;
		}
	
		friend std::ostream& operator<<(std::ostream& os, const Projection& obj)
		{
	   	char temp[100];
	   	sprintf(temp,"%0.2lfx + %0.2lfy <= %0.2lf",obj.a,obj.b,obj.c);
	   	os << temp;
	   	return os;
		}
	
	};


In the above data structures, a,b,c represent the equation of line in the inequality 'ax + by <= c' <br>
m represents the slope of line


We observe there can be four possibilities in Intersection of Projection of all half planes from 0 to i-1 on the line
formed from ith inequality.

Projection can be an entire LINE <br>
Initially, a Projection object is initiated to a LINE and is then projected onto half planes 0 to i-1

Projection can be an EMPTY (NULLSET) <br>
Incase the intersection is empty if the feasible region itself is null

Projection can be an RAY <br>
To represent a ray we have the starting point as p. 
The direction of ray whether upwards of downwards is represented by dir field. 
proj.dir = UP or DOWN

Projection can be an LINESEGMENT <br>
Line segment is represented by two of its extreme points p and q.
here p should be below q, i.e. y-coordinate of p should be less than y-coordinate of q


While taking projection of half-planes on line different cases arise and have been handled.
