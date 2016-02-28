/*
BANDI SUMANTH
11CS30006
Computational Geometry 2015
Assignment 1 - Implementation of Chan's algorithm for computing convex hull
January 27th, 2014
IIT Kharagpur



Assignment 1 - Implementation of Chan's algorithm for computing convex hull
---------------------------------------------------------------------------
*/

#include <vector>
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

using namespace std;


#define s(n) scanf("%d",&n)
 
#define ll long long
#define all(c) c.begin(), c.end()
#define pb push_back
#define ff first
#define ss second
 
#define rep(i, n) for(int i = 0; i < (n); ++i)

// typedef vector<int> vi;
// typedef vector< vector<int> > vvi;
typedef pair<int,int> pii;
typedef std::vector<pii> vpii;

#define sz(a) int((a).size()) 
#define szar(ar) int(sizeof(ar)/sizeof(ar[0]))

//-------------------------------------------------------------------------------------

const int COLLINEAR = 0;
const int CLOCKWISE = -1;
const int COUNTER_CLOCKWISE = 1;

struct Point
{
    int x;
    int y;

    Point() {}

    Point(int _x, int _y) : x(_x), y(_y) {}
};

bool operator < (Point const& lhs, Point const& rhs)
{
  if(lhs.x != rhs.x) return lhs.x < rhs.x;
  else if(lhs.y != rhs.y) return lhs.y < rhs.y;
  else {
    assert(lhs.x==rhs.x && lhs.y==rhs.y);
    // printf("POINT (%d,%d) APPEARS TWICE !!!\n",lhs.x, lhs.y);
    return false;
  }
}

typedef vector<Point> vpt;
typedef vector<vpt > vvpt;

// A globle point needed for  sorting points with reference to the first point
// Used in compare function for sorting in Graham scan algo
Point p0;
 
// find next to top in a stack
Point nextToTop(stack<Point> &S)
{
    Point p = S.top();
    S.pop();
    Point res = S.top();
    S.push(p);
    return res;
}
 
// returns square of distance between p1 and p2
int dist(Point p1, Point p2)
{
    return (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y);
}
 
// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// -1 --> Clockwise
// 1 --> Counterclockwise
int orientation(Point p, Point q, Point r)
{
    int val = (q.y - p.y) * (r.x - q.x) -
              (q.x - p.x) * (r.y - q.y);
 
    if (val == 0) return 0;  // colinear
    return (val > 0)? CLOCKWISE: COUNTER_CLOCKWISE; // clock or counterclock wise
}

// A custom compare function used by sort() to sort a vector of
// Points with respect to the first point p0
bool compare_with_p0(const Point& p1, const Point& p2) {
  int o = orientation(p0,p1,p2);
  if (o == COLLINEAR)  // the nearest point should come before the farther point in the sorted list
    return dist(p0, p2) >= dist(p0, p1);
  
  return (o == COUNTER_CLOCKWISE);
}

void reorderhull(vpt & hull) {
  // after reordering , leftmost point will be the first point

  int minidx = min_element(all(hull)) - hull.begin();

  vpt temp(hull.begin()+minidx, hull.end());

  for(int i=0; i<minidx; ++i) {
    temp.push_back(hull[i]);
  }

  hull.swap(temp);
}


vpt GrahamScan(vpt points)
{
  int n = sz(points);

  vpt hull;
  // Find the bottommost point
  int ymin = points[0].y, min = 0;
  for (int i = 1; i < n; i++)
  {
    int y = points[i].y; 
    // Pick the bottom-most or chose the left most point in case of tie
    if ((y < ymin) || (ymin == y && points[i].x < points[min].x))
      ymin = points[i].y, min = i;
  } 
  // Place the bottom-most point at first position
  swap(points[0], points[min]); 
  // Sort n-1 points with respect to the first point.  A point p1 comes
  // before p2 in sorted ouput if p2 has larger polar angle (in 
  // counterclockwise direction) than p1
  p0 = points[0];
  sort(points.begin()+1,points.end(),compare_with_p0);

  // EDGE ERROR CASE (n == 2)
  if(n==2) return points;

  // Create an empty stack and push first three points to it.
  stack<Point> S;
  S.push(points[0]);
  S.push(points[1]);
  S.push(points[2]); 

  // Process remaining n-3 points
  for (int i = 3; i < n; i++)
  {
    // Keep removing top while the angle formed by points next-to-top, 
    // top, and points[i] makes a non-left turn
    while (sz(S)>=2 && orientation(nextToTop(S), S.top(), points[i]) != COUNTER_CLOCKWISE)
        S.pop();
    S.push(points[i]);
  } 
  // Now stack has the output points
  while (!S.empty())
  {
      Point p = S.top();
      hull.push_back(Point(p.x,p.y));
      S.pop();
  }

  reverse(all(hull));  // after reversing, hull has Anti-clockwise list of points

  // remove the 2nd point if first 3 points are collinear points
  if(sz(hull)>=3 && orientation(hull[0],hull[1],hull[2]) == COLLINEAR) hull.erase(hull.begin()+1);

  reorderhull(hull);
  return hull;
}

pii min_hull_pt_pair(vvpt hulls) {
  // Returns the (hull index, point index) pair of the leftmost (and bottom most under a tie) point of all points
  int h=0, p=0;
  rep(i,sz(hulls)) {
    int j = min_element(all(hulls[i])) - hulls[i].begin();
    if(hulls[i][j] < hulls[h][p]) {
      h = i;
      p = j;
    }
  }
  return pii(h,p);
}

int rtangent(vpt& hull, Point p) {
  // Return the index of the point in hull that the right tangent line from p to hull touches.
  int size = sz(hull);
  int l = 0, r = sz(hull);
  int l_prev = orientation(p, hull[0], hull.back());
  int l_next = orientation(p, hull[0], hull[(l+1)%r]);

  while(l < r) {
    int c = (l + r)/2;
    int c_prev = orientation(p, hull[c], hull[(c-1+size) % size]);
    int c_next = orientation(p, hull[c], hull[(c+1) % size]);
    int c_side = orientation(p, hull[l], hull[c]);

    if(c_prev != CLOCKWISE && c_next != CLOCKWISE)
      return c;
    else if((c_side == COUNTER_CLOCKWISE && (l_next == CLOCKWISE || l_prev == l_next)) || 
              (c_side == CLOCKWISE && c_prev == CLOCKWISE))
      r = c;
    else {
      l = c+1;
      l_prev = (-1)*c_next;
      l_next = orientation(p, hull[l], hull[(l+1)%size]);
    }
  }
  return l;
}

pii next_hull_pt_pair(vvpt& hulls, pii& pair) {
  // Returns the (hull index, point index) pair of the next point in the convex hull.

  Point p = hulls[pair.ff][pair.ss];

  // Finding the next point for current hull is a little easier.
  pii next = pii(pair.ff, (pair.ss + 1) % sz(hulls[pair.ff]));

  for (int h = 0; h < sz(hulls); ++h)
  {
    if(h == pair.ff) continue;

    int s = rtangent(hulls[h], p);

    Point q,r;
    q = hulls[next.ff][next.ss];
    r = hulls[h][s];

    int t = orientation(p,q,r);
    if(t == CLOCKWISE || (t == COLLINEAR && dist(p,r) > dist(p,q)))
      next = pii(h,s);
  }
  return next;
}

bool ChanHull(vpt & P, int m, vpii& hull, vvpt & hulls) {
  printf("\n\nTaking m = H = %d\n",m );

  hull.clear();
  hulls.clear();
  int H = m;
  int n = sz(P);

  for (int i = 0; i < n; i+=m)
  { 
    vpt subhull = GrahamScan(vpt(P.begin()+i, P.begin()+min(i+m,n)));
    hulls.push_back(subhull);

    printf("Intermediate hull of Points[%d,..,%d] :  ",i,min(i+m,n)-1);
    rep(j,sz(subhull)) {
      if(j!=sz(subhull)-1) printf("(%d,%d), ",subhull[j].x, subhull[j].y );
      else printf("(%d,%d)  ",subhull[j].x, subhull[j].y );
    }
    printf("\n");
  }

  // Here we find the extreme point (leftmost) and initialize our hull with it.
  hull.push_back(min_hull_pt_pair(hulls));

  // We must ensure we stop after no more than H iterations.
  rep(iteration, H) {
    pii p = next_hull_pt_pair(hulls,hull.back());
    if(p == hull[0]) return true;
    hull.push_back(p);
  }
  return false;
}

vpt ChanHull(vpt & P) {
  vpii hull;
  vvpt hulls;

  int n = sz(P);
  for (int t = 1; t < n; ++t)
  {
    int m = min((int)pow(2, pow(2,t)), n);
    assert(m>=1 && m<=LLONG_MAX);

    if(ChanHull(P,m,hull,hulls)) {
      printf("\nCOMPLETE : Convex hull construction is Complete\n");
      vpt L;
      rep(i,sz(hull)) L.push_back(hulls[hull[i].ff][hull[i].ss]);
      return L;
    }

  }
}


 
int main()
{  
  vpt pts, hull;

  int n;
  scanf("%d",&n);

  for(int i = 0; i < n; ++i) {
    int x,y;
    scanf("%d",&x);
    scanf("%d",&y);
    pts.push_back(Point(x,y));
  }

  printf("Chan Convex Hull Algorithm :\n\nThe Intermediate hulls are found using Graham's Scan technique.\n");
  hull = ChanHull(pts);

  printf("\nResult:\n");
  printf("------\n");
  printf("Number of points in Convex Hull = %d\n",sz(hull));
  printf("\nThe list of Convex Hull points in counterclockwise order : \n{ ");
  for(int i = 0; i < sz(hull); ++i) {
    if(i!=sz(hull)-1) printf("(%d,%d), ",hull[i].x, hull[i].y);
    else printf("(%d,%d) }",hull[i].x, hull[i].y);
  }
  printf("\n");

  /*printf("\n\nThrough GrahamScan Algo :\n");
  hull = GrahamScan(pts);
  printf("Size of Convex Hull = %d\n{ ",sz(hull));
  rep(i,sz(hull)) {
    printf("(%d,%d), ",hull[i].x, hull[i].y);
  }
  printf(" }\n");*/

  return 0;
}
