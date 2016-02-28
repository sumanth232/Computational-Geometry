/*
BANDI SUMANTH
11CS30006
Computational Geometry 2015
Assignment 5 - Randomized Incremental 2D Linear Programming algorithm
March 27th, 2014
IIT Kharagpur



Assignment 5 - Randomized Incremental 2D Linear Programming algorithm
---------------------------------------------------------------------
*/

#include <vector>
#include <set>
#include <stack>
#include <algorithm>
#include <functional>
#include <numeric>
#include <utility>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>
#include <cstring>
#include <cassert>
#include <climits>
#include <cfloat>
#include <iomanip>

using namespace std;

#define s(n) scanf("%d",&n)
 
#define ll long long
#define all(c) c.begin(), c.end()
#define pb push_back
#define ff first
#define ss second
 
#define rep(i, n) for(int i = 0; i < (n); ++i)
#define tr(c,it) for(typeof(c.begin()) it = c.begin(); it !=c.end(); it++)
#define cpresent(c,x) ((c).find(x) != (c).end()) 

// typedef vector<int> vi;
// typedef vector< vector<int> > vvi;
typedef pair<int,int> pii;
typedef std::vector<pii> vpii;

#define sz(a) int((a).size()) 
#define szar(ar) int(sizeof(ar)/sizeof(ar[0]))


//-------------------------------------------------------------------------------------


const int LINE = 1;
const int RAY = 2;
const int LINESEGMENT = 3;
const int EMPTY = 4;
const int UP = 11;
const int DOWN = 12;

const double h = 1e-4;

const double MAX = 10000;
const double MIN = 1e-4;



struct Halfplane
{
	// Half palne is represented by inequality
	//			ax + by <= c
	double a,b,c;
	double m;

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

typedef std::vector<Halfplane> vhp;

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

struct Projection
{
	Point p,q;	// q should not be accessed if it is a RAY
				// p  should be below q if it is a LINESEGMENT
	double a,b,c;
	double m;	// slope
	int type;
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



// ------------helper functions--------------------------

// A function to generate a random permutation of arr[]
void randomize ( vhp & arr, int n )
{
    // Use a different seed value so that we don't get same
    // result each time we run this program
    // srand ( time(NULL) );
 
    // Start from the last element and swap one by one. We don't
    // need to run for the first element that's why i > 0
    for (int i = n-1; i > 0; i--)
    {
        // Pick a random index from 0 to i
        int j = rand() % (i+1);
 
        // Swap arr[i] with the element at random index
        swap(arr[i], arr[j]);
    }
}

double det(double a,double b,double c,double d)
{
    return a*d-b*c;
}
Point intersect(double a1,double b1,double c1,double a2,double b2,double c2)
{   
	assert(a1/b1 != a2/b2);
    Point r;
    double d,e,f;
    f=det(a1,b1,a2,b2);
    d=det(c1,b1,c2,b2);
    e=det(a1,c1,a2,c2);
    r.x=d/f;
    r.y=e/f;
    return r;
}

Point intersect(const Halfplane& x, const Halfplane& y) {
	return intersect(x.a,x.b,x.c,y.a,y.b,y.c);
}

Point intersect(const Halfplane& x, const Projection& proj) {
	return intersect(x.a,x.b,x.c,proj.a,proj.b,proj.c);
}

bool isLies(Halfplane hp, Point p) {
	// cout << "diff("  << hp << "," << p << ") = " << hp.a*p.x + hp.b*p.y - hp.c << endl;
	return (hp.a*p.x + hp.b*p.y <= hp.c);
}

bool isLiesOnLineSgement(Projection proj, Point p) {
	if(proj.a*p.x + proj.b*p.y <= proj.c) {
		// it should also lie between the two end points
		double miny = min(proj.p.y, proj.q.y);
		double maxy = max(proj.p.y, proj.q.y);
		return (p.y >= miny && p.y <= maxy);
	}
	else return false;
}

Projection Project_halfplane(const Halfplane& hp, Projection proj){
	Projection new_proj;

	if(proj.type == EMPTY) return proj;

	Point ip = intersect(hp,proj);
	char temp[200];

	cout << "ip (" <<  hp << "," << proj << ")= " << ip << endl;

	// for direction
	double m = proj.m;
	double upx,upy,downx,downy;

	assert(m != 0);
	double x = ip.x, y = ip.y;
	if(m >= 0) {
		upx = x + h; upy = y + m*h;
		downx = x - h; downy = y - m*h;
	}
	else {
		upx = x - h; upy = y - m*h;
		downx = x + h; downy = y + m*h;
	}


	if(proj.type == EMPTY) {
		return proj;
	}
	else if(proj.type == LINE) {
		new_proj.p = ip;
		new_proj.setabc(proj);
		new_proj.m = proj.m;
		new_proj.type = RAY;

		if(isLies(hp,Point(upx,upy))) new_proj.dir = UP;
		else new_proj.dir = DOWN;
		return new_proj;
	}
	else if(proj.type == LINESEGMENT) {
		// if intersection point lies on the linesegment
		// if((proj.a*ip.x + proj.b*ip.y == proj.c) && ()) {
		if(isLiesOnLineSgement(proj,ip)) {
			printf("ip lies on LINESEGMENT\n");
			
			new_proj.setabc(proj);
			new_proj.m = proj.m;
			new_proj.type = LINESEGMENT;

			// fiirst blindly copy p,q of old LINESEGMENT
			new_proj.p = proj.p;
			new_proj.q = proj.q;


			// now correct accordingly
			if(isLies(hp,Point(upx,upy))) {
				// new_proj.dir = UP;
				new_proj.p = ip;
			}
			else {
				// new_proj.dir = DOWN;
				new_proj.q = ip;
			}
			return new_proj;
		}
		else {
			/* 	both points lie on same side of halfplane
			   	so the result is either EMPTY or a LINESEGMENT*/
			/* 	so just check one of the points, and if that lies in the halfplane , the result is 
				is the whole LINESEGMENT */
			if(hp.a*proj.p.x + hp.b*proj.p.y <= hp.c) {
				return proj;
			}
			else {
				new_proj.type =EMPTY;
				return new_proj;
			}

		}
	}
	else if(proj.type == RAY) {
		// IF halfplane and proj are pointing in opposite directions
		cout << "proj dir of RAY = " << proj.dir << endl;
		cout << "isLies(hp,Point(downx,downy)) = " << isLies(hp,Point(downx,downy)) << endl;
		cout << "isLies(hp,Point(upx,upy)) = " << isLies(hp,Point(upx,upy)) << endl;
		cout << "isLies(hp,Point(ip.x,ip.y)) = " << isLies(hp,Point(ip.x,ip.y)) << endl;

		// RAY = UP, HP = DOWN
		if(proj.dir == UP && isLies(hp,Point(downx,downy))) {
			//if the down point lies on the RAY
			if(proj.a*downx + proj.b*downy == proj.c) {
				// the result is a LINESEGMENT
				new_proj.type = LINESEGMENT;

				assert(ip.y != proj.p.y);
				if(ip.y <= proj.p.y) {
					new_proj.p = ip;
					new_proj.q = proj.p;
				}
				else {
					new_proj.p = proj.p;
					new_proj.q = ip;
				}
				
				new_proj.setabc(proj);
				new_proj.m = proj.m;
				return new_proj;
			}
			else {
				new_proj.type = EMPTY;
				return new_proj;
			}
		}
		else if(proj.dir == DOWN && isLies(hp,Point(upx,upy))) {
			// if the up point lies on the ray
			if(proj.a*upx + proj.b*upy == proj.c) {
				// the result is a LINESEGMENT
				new_proj.type = LINESEGMENT;

				assert(ip.y != proj.p.y);
				if(ip.y <= proj.p.y) {
					new_proj.p = ip;
					new_proj.q = proj.p;
				}
				else {
					new_proj.p = proj.p;
					new_proj.q = ip;
				}

				new_proj.setabc(proj);
				new_proj.m = proj.m;
				return new_proj;
			}
			else {
				new_proj.type = EMPTY;
				return new_proj;
			}
		}

		// IF halfplane and proj are pointing in same directions

		// both UP or both DOWN
		else if((proj.dir == UP && isLies(hp,Point(upx,upy))) || (proj.dir == DOWN && isLies(hp,Point(downx,downy)))) {
			// if entire ray lies in the halfplane
			if(isLies(hp,Point(proj.p.x, proj.p.y))) {
				return proj;
			}
			else {
				proj.p = ip;
				return proj;
			}
		}
		else {
			printf("unknown case in RAY occurred\n");
			assert(false);
		}
		
	}
	else {
		printf("unknown projection type occurred\n");
		assert(false);
	}

}

double c1,c2;

void InsertHalfplanes(vhp & halfplanes) {
	if(c1 >=0 && c2 >= 0) {
		halfplanes.insert(halfplanes.begin(), Halfplane(MIN,1,MAX));	// y = max
		halfplanes.insert(halfplanes.begin(), Halfplane(1,MIN,MAX));	// x = max
	}
	else if(c1 <= 0 && c2 >= 0) {
		halfplanes.insert(halfplanes.begin(), Halfplane(MIN,1,MAX));	// y = max
		halfplanes.insert(halfplanes.begin(), Halfplane(-1,MIN,MAX));	// x = -max
	}
	else if(c1 <= 0 && c2 <= 0) {
		halfplanes.insert(halfplanes.begin(), Halfplane(MIN,-1,MAX));	// y = -max
		halfplanes.insert(halfplanes.begin(), Halfplane(-1,MIN,MAX));	// x = -max
	}
	else if(c1 >= 0 && c2 <= 0) {
		halfplanes.insert(halfplanes.begin(), Halfplane(MIN,-1,MAX));	// y = -max
		halfplanes.insert(halfplanes.begin(), Halfplane(1,MIN,MAX));	// x = max
	}
}





int main(int argc, char const *argv[])
{
	int n;
	
	vhp halfplanes;
	printf("halfplane is of the form ax + by <= c\n");
	printf("Enter number of half planes : ");
	scanf("%d",&n);
	for (int i = 0; i < n; ++i)
	{
		double a,b,c;
		scanf("%lf %lf %lf",&a, &b, &c);
		halfplanes.push_back(Halfplane(a,b,c));
	}

	printf("Enter c1,c2 of the objective function c1x + c2y : \n");
	scanf("%lf %lf",&c1, &c2);

	/*for (int i = 0; i < n; ++i)
	{
		cout << halfplanes[i] << endl;
	}*/

	randomize(halfplanes,n);

	printf("After randomizing\n");

	for (int i = 0; i < n; ++i)
	{
		cout << halfplanes[i] << endl;
	}

	InsertHalfplanes(halfplanes);
	n += 2;

	assert(n >= 2);
	Point v = intersect(halfplanes[0], halfplanes[1]);

	cout << "Initial intersection point = " << v << endl;

	for (int i = 2; i < n; ++i)
	{
		double a,b,c;
		a = halfplanes[i].a;
		b = halfplanes[i].b;
		c = halfplanes[i].c; 
		if (!(a*v.x + b*v.y <= c)) {
			cout << v << " does not lie in " << halfplanes[i] << endl;

			Projection proj;
			proj.type = LINE;
			proj.a = a;
			proj.b = b;
			proj.c = c;
			proj.m = halfplanes[i].m;

			for (int j = 0; j < i; ++j)
			{
				proj = Project_halfplane(halfplanes[j], proj);

				cout << "proj type = " << proj.type;
				if(proj.type == LINESEGMENT) {
					cout << " : " << proj.p << " , " << proj.q ;
				}
				cout << endl;
			}

			if(proj.type == EMPTY) {
				printf("INFEASIBLE\n");
				return 0;
			}
			else {
				cout << "else projtype = " << proj.type << endl;
				// assuming no parallel halfplanes, proj is a LINESEGMENT or a RAY
				if(proj.type == RAY) {
					v = proj.p;
				}
				else {
					assert(proj.type != LINE);
					// its a LINESEGMENT
					cout << "The 2 points of the LINESEGMENT : " << proj.p << " , " << proj.q << endl;
					if(c1*proj.p.x + c2*proj.p.y >= c1*proj.q.x + c2*proj.q.y) {
						v = proj.p;
					}
					else v = proj.q;
				}
			}

			cout << "NEW v = " << v << endl;
		}
		else {
			cout << v << " does LIES in " << halfplanes[i] << endl;
		}
	}

	printf("Feasible point P that maximizes the dot product between objective function and the point (C.P) :\n");

	printf("x = %0.2lf, y = %0.2lf\n", v.x,v.y);
	printf("%0.2lfx + %0.2lfy = %0.2lf\n",c1,c2, c1*v.x+c2*v.y);


	return 0;
}

