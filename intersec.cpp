/*
BANDI SUMANTH
11CS30006
Computational Geometry 2015
Assignment 2 - Line Segment Intersection
February 7th, 2014
IIT Kharagpur



Assignment 2 - Line Segment Intersection
----------------------------------------
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

// -----------CHAN HULL Algorithm code  START--------------------------------------------------------------



const int COLLINEAR = 0;
const int CLOCKWISE = -1;
const int COUNTER_CLOCKWISE = 1;

struct myPoint
{
    double x;
    double y;

    myPoint() {}

    myPoint(double _x, double _y) : x(_x), y(_y) {}

    bool operator < (const myPoint & p) const {
      if(x != p.x) return x < p.x;
      else return y < p.y; 
    }

};

/*bool operator < (myPoint & lhs, myPoint & rhs)
{
  if(lhs.x != rhs.x) return lhs.x < rhs.x;
  else return lhs.y < rhs.y;
}*/

typedef vector<myPoint> vmypt;
typedef vector<vmypt > vvmypt;

// A globle point needed for  sorting points with reference to the first point
// Used in compare function for sorting in Graham scan algo
myPoint p0;
 
// find next to top in a stack
myPoint nextToTop(stack<myPoint> &S)
{
    myPoint p = S.top();
    S.pop();
    myPoint res = S.top();
    S.push(p);
    return res;
}
 
// returns square of distance between p1 and p2
double dist(myPoint p1, myPoint p2)
{
    return (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y);
}
 
// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// -1 --> Clockwise
// 1 --> Counterclockwise
int orientation(myPoint p, myPoint q, myPoint r)
{
    double val = (q.y - p.y) * (r.x - q.x) -
              (q.x - p.x) * (r.y - q.y);
 
    if (val == 0) return 0;  // colinear
    return (val > 0)? CLOCKWISE: COUNTER_CLOCKWISE; // clock or counterclock wise
}

// A custom compare function used by sort() to sort a vector of
// Points with respect to the first point p0
bool compare_with_p0(const myPoint& p1, const myPoint& p2) {
  int o = orientation(p0,p1,p2);
  if (o == COLLINEAR)  // the nearest point should come before the farther point in the sorted list
    return dist(p0, p2) >= dist(p0, p1);
  
  return (o == COUNTER_CLOCKWISE);
}

void reorderhull(vmypt & hull) {
  // after reordering , leftmost point will be the first point

  int minidx = min_element(all(hull)) - hull.begin();

  vmypt temp(hull.begin()+minidx, hull.end());

  for(int i=0; i<minidx; ++i) {
    temp.push_back(hull[i]);
  }

  hull.swap(temp);
}


vmypt GrahamScan(vmypt points)
{
  int n = sz(points);

  vmypt hull;
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
  stack<myPoint> S;
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
      myPoint p = S.top();
      hull.push_back(myPoint(p.x,p.y));
      S.pop();
  }

  reverse(all(hull));  // after reversing, hull has Anti-clockwise list of points

  // remove the 2nd point if first 3 points are collinear points
  if(sz(hull)>=3 && orientation(hull[0],hull[1],hull[2]) == COLLINEAR) hull.erase(hull.begin()+1);

  reorderhull(hull);
  return hull;
}

pii min_hull_pt_pair(vvmypt hulls) {
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

int rtangent(vmypt& hull, myPoint p) {
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

pii next_hull_pt_pair(vvmypt& hulls, pii& pair) {
  // Returns the (hull index, point index) pair of the next point in the convex hull.

  myPoint p = hulls[pair.ff][pair.ss];

  // Finding the next point for current hull is a little easier.
  pii next = pii(pair.ff, (pair.ss + 1) % sz(hulls[pair.ff]));

  for (int h = 0; h < sz(hulls); ++h)
  {
    if(h == pair.ff) continue;

    int s = rtangent(hulls[h], p);

    myPoint q,r;
    q = hulls[next.ff][next.ss];
    r = hulls[h][s];

    int t = orientation(p,q,r);
    if(t == CLOCKWISE || (t == COLLINEAR && dist(p,r) > dist(p,q)))
      next = pii(h,s);
  }
  return next;
}

bool ChanHull(vmypt & P, int m, vpii& hull, vvmypt & hulls) {
  // printf("\n\nTaking m = H = %d\n",m );



  hull.clear();
  hulls.clear();
  int H = m;
  int n = sz(P);

  for (int i = 0; i < n; i+=m)
  { 
    vmypt subhull = GrahamScan(vmypt(P.begin()+i, P.begin()+min(i+m,n)));
    hulls.push_back(subhull);

    /*printf("Intermediate hull of Points[%d,..,%d] :  ",i,min(i+m,n)-1);
    rep(j,sz(subhull)) {
      if(j!=sz(subhull)-1) printf("(%0.2lf,%0.2lf), ",subhull[j].x, subhull[j].y );
      else printf("(%0.2lf,%0.2lf)  ",subhull[j].x, subhull[j].y );
    }
    printf("\n");*/
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

vmypt ChanHull(vmypt & P) {
  if(sz(P) == 0) return vmypt();

  vpii hull;
  vvmypt hulls;

  int n = sz(P);
  for (int t = 1; t < n; ++t)
  {
    int m = min((int)pow(2, pow(2,t)), n);
    assert(m>=1 && m<=LLONG_MAX);

    if(ChanHull(P,m,hull,hulls)) {
      // printf("\nCOMPLETE : Convex hull construction is Complete\n");
      vmypt L;
      rep(i,sz(hull)) L.push_back(hulls[hull[i].ff][hull[i].ss]);
      return L;
    }

  }
}

// -----------CHAN HULL Algorithm code  END--------------------------------------------------------------



//-------------------------------------------------------------------------------------
const int LEFT_END = 0;
const int INTERSECTION = 1;
const int RIGHT_END = 2;
const double h = 1e-5;


struct LineSegment
{
    double x1,y1,x2,y2,m,c;
    int idx;
    LineSegment() {}
    // LineSegment(double _x1, double _y1, double _x2, double _y2) : x1(_x1), y1(_y1), x2(_x2), y2(_y2) {}
    LineSegment(double _x1, double _y1, double _x2, double _y2, double _m, double _c, int _idx) : x1(_x1), y1(_y1), x2(_x2), y2(_y2), m(_m), c(_c), idx(_idx) {
        
    }

    friend ostream &operator<<( ostream &output, const LineSegment &p )
    {
        /*char temp[500];
        sprintf(temp,"{(%0.2lf,%0.2lf), (%0.2lf,%0.2lf)}",p.x1,p.y1,p.x2,p.y2);
        string str(temp);
        output << temp;*/

        output << "{(" << p.x1 << "," << p.y1 << "), (" <<  p.x2 << "," << p.y2 << ")}";

        return output;
    }  
};

typedef vector<LineSegment> vls;

vls lineseg;
double currx;   // THIS SHOULD BE SET

struct Label
{
    int idx;

    Label() {}

    Label(const Label & label) {
        this->idx = label.idx;
    }

    Label(int _idx) : idx(_idx) {
        // pricntf("called with idx %d\n", idx);
    }

    bool operator < (const Label & label) const
    {
        // DEGENERACY NOT HANDLED
        // WHAT IS BOTH HAVE SAME Y-CORDINATE
        // pricntf("calling Lebel <\n");
        double y = (lineseg[idx].m)*(currx) + lineseg[idx].c;
        double labely = (lineseg[label.idx].m)*(currx) + lineseg[label.idx].c;

        // pricntf("thisy = %lf, othery = %lf\n",y,labely );

        /*if(y != labely) return y > labely;
        else {
            double thisSlope = lineseg[idx].m, otherSlope = lineseg[label.idx].m; 

            if(thisSlope != otherSlope) {
                if(thisSlope*otherSlope < 0) return thisSlope < 0;
                else {
                    // both have same sign
                    if(thisSlope > 0) return thisSlope < otherSlope;
                    else thisSlope > otherSlope;
                }
            }
            else return false;
        }*/
        return y>labely;
    }

    bool operator == (const Label & label) const
    {
        // DEGENERACY NOT HANDLED
        // WHAT IS BOTH HAVE SAME Y-CORDINATE
        double y = (lineseg[idx].m)*(currx) + lineseg[idx].c;
        double labely = (lineseg[label.idx].m)*(currx) + lineseg[label.idx].c;
        return y == labely;
    }

    const Label & operator = (const Label & label)
    {
        this->idx = label.idx;
        return *this;
    }


    bool operator > (const Label & label) const
    {
        return (!(*(this) == label)) && (!(*(this) < label));
    }


};


struct Point
{
    double x;
    double y;

    Label label;
    int type;

    Label l1,l2;  // labels of intersecting LineSegments

    Point() {}

    Point(const Point & p) {
        this->x = p.x;
        this->y = p.y;
        this->label = p.label;
        this->type = p.type;
        this->l1 = p.l1;
        this->l2 = p.l2;
    }

    // Point(double _x, double _y, int _idx) : x(_x), y(_y), idx(_idx) {}
    Point(double _x, double _y, Label l, int _type) : x(_x), y(_y), type(_type) {
        label.idx = l.idx;
    }
    Point(double _x, double _y) : x(_x), y(_y) {}

    const Point & operator = (const Point & p)
    {
      this->x = p.x;
      this->y = p.y;
      this->label = p.label;
      this->type = p.type;
      this->l1 = p.l1;
      this->l2 = p.l2;
      return *this;
    }

    friend ostream &operator<<( ostream &output, const Point &p )
    {
      if(p.type == INTERSECTION) {
        LineSegment a = lineseg[p.l1.idx], b = lineseg[p.l2.idx];
        char temp[500];
        sprintf(temp,"(%0.2lf,%0.2lf) of line segments {(%0.1lf,%0.1lf), (%0.1lf,%0.1lf)} , {(%0.1lf,%0.1lf), (%0.1lf,%0.1lf)} with line numbers %d, %d",
                      p.x,p.y,a.x1,a.y1,a.x2,a.y2,b.x1,b.y1,b.x2,b.y2,p.l1.idx+1,p.l2.idx+1);
        string str(temp);
        output << temp;
        
        /*output << setiosflags(ios::fixed) << setprecision(2) << "(" << p.x << "," << p.y << ")" ;
        output << resetiosflags(ios::fixed) << " of line segments " << a << " , " << b << " with the line numbers " << p.l1.idx+1 << ", " << p.l2.idx+1;*/

        return output;
      }
      else {
        output << "(" << p.x << "," << p.y << ")";
        return output;
      }

      
    }   

    bool operator < (const Point & p) const {
      if(x != p.x) return x < p.x;
      else return y > p.y; 
    }                                 
};


typedef set<Point> spt;
typedef vector<Point> vpt;
typedef vector<vpt > vvpt;

void cprint(spt& c) {
    // printf("Event queue\n");
    tr(c,it) {
        printf("%lf, %lf\n",(*it).x, (*it).y );
    }
}

// returns square of distance between p1 and p2
double dist(Point p1, Point p2)
{
    return (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y);
}




// An AVL tree node
struct node
{
    Label key;
    struct node *left;
    struct node *right;
    int height;
};
 
// get maximum of two integers
int max(int a, int b);
 
// get height of the tree
int height(struct node *N)
{
    if (N == NULL)
        return 0;
    return N->height;
}



 
// get maximum of two integers
int max(int a, int b)
{
    return (a > b)? a : b;
}
 
struct node* newNode(Label key)
{
    struct node* node = (struct node*)
                        malloc(sizeof(struct node));
    node->key   = key;
    node->left   = NULL;
    node->right  = NULL;
    node->height = 1;  // new node is initially added at leaf
    return(node);
}
 
// right rotate subtree rooted with y
struct node *rightRotate(struct node *y)
{
    struct node *x = y->left;
    struct node *T2 = x->right;
 
    // Perform rotation
    x->right = y;
    y->left = T2;
 
    // Update heights
    y->height = max(height(y->left), height(y->right))+1;
    x->height = max(height(x->left), height(x->right))+1;
 
    // Return new root
    return x;
}
 
// left rotate subtree rooted with x
struct node *leftRotate(struct node *x)
{
    struct node *y = x->right;
    struct node *T2 = y->left;
 
    // Perform rotation
    y->left = x;
    x->right = T2;
 
    //  Update heights
    x->height = max(height(x->left), height(x->right))+1;
    y->height = max(height(y->left), height(y->right))+1;
 
    // Return new root
    return y;
}
 
// get height difference of node N
int getHeightDifference(struct node *N)
{
    if (N == NULL)
        return 0;
    return height(N->left) - height(N->right);
}
 
struct node* insert(struct node* node, Label key)
{
    // BST rotation 
    if (node == NULL)
        return(newNode(key));
 
    if (key < node->key)
        node->left  = insert(node->left, key);
    else
        node->right = insert(node->right, key);
 
    // update height 
    node->height = max(height(node->left), height(node->right)) + 1;
 
    int heightDiff = getHeightDifference(node);
 
    // If this node becomes unbalanced, then there are 4 cases
 
    // Left Left Case
    if (heightDiff > 1 && key < node->left->key)
        return rightRotate(node);
 
    // Right Right Case
    if (heightDiff < -1 && key > node->right->key)
        return leftRotate(node);
 
    // Left Right Case
    if (heightDiff > 1 && key > node->left->key)
    {
        node->left =  leftRotate(node->left);
        return rightRotate(node);
    }
 
    // Right Left Case
    if (heightDiff < -1 && key < node->right->key)
    {
        node->right = rightRotate(node->right);
        return leftRotate(node);
    }
 
    return node;
}
 
/* returns the node with minimum
   key value in the tree */
struct node * minValueNode(struct node* node)
{
    struct node* current = node;
 
    /* loop down to the leftmost leaf */
    while (current->left != NULL)
        current = current->left;
 
    return current;
}
 
struct node* deleteNode(struct node* root, Label key)
{
    if (root == NULL)
        return root;
 
    if ( key < root->key ) {
        // printf("less than root %d %d \n",key.idx, root->key.idx);
        root->left = deleteNode(root->left, key);
    }
        
    else if( key > root->key )
        root->right = deleteNode(root->right, key);
 
    // if key is same as root's key
    else
    {
        // node with only one child or no child
        if( (root->left == NULL) || (root->right == NULL) )
        {
            struct node *temp = root->left ? root->left : root->right;
 
            // No child case
            if(temp == NULL)
            {
                temp = root;
                root = NULL;
            }
            else // One child case
             *root = *temp; 
 
            free(temp);
        }
        else
        {
            // node with two children: Get the inorder successor (smallest in the right subtree)
            struct node* temp = minValueNode(root->right);
 
            // Copy the inorder successor's data to this node
            root->key = temp->key;
 
            // delete the inorder successor
            root->right = deleteNode(root->right, temp->key);
        }
    }
 
    if (root == NULL)
      return root;
 
    // update height
    root->height = max(height(root->left), height(root->right)) + 1;
 
    int balance = getHeightDifference(root);
 
    // If this node becomes unbalanced, then there are 4 cases
 
    // Left Left Case
    if (balance > 1 && getHeightDifference(root->left) >= 0)
        return rightRotate(root);
 
    // Left Right Case
    if (balance > 1 && getHeightDifference(root->left) < 0)
    {
        root->left =  leftRotate(root->left);
        return rightRotate(root);
    }
 
    // Right Right Case
    if (balance < -1 && getHeightDifference(root->right) <= 0)
        return leftRotate(root);
 
    // Right Left Case
    if (balance < -1 && getHeightDifference(root->right) > 0)
    {
        root->right = rightRotate(root->right);
        return leftRotate(root);
    }
 
    return root;
}
 
// print inorder traversal of the tree.
void inorder(struct node *root)
{
    if(root != NULL)
    {
        
        inorder(root->left);
        printf("%d ", (root->key).idx);
        inorder(root->right);
    }
}

struct node* search(struct node * root, Label key) {
    assert(root != NULL);
    // if(root == NULL) return NULL;

    if((root->key).idx == key.idx) return root;
    else if(key > (root->key)) {
        
        assert(root->right != NULL);
        // printf("going Right | root = %d\n",root->right->key.idx);
        return search(root->right, key);
    }
    else {
        
        assert(root->left != NULL);
        // printf("going Left | root = %d\n",root->left->key.idx);
        return search(root->left, key);
    }
    // assert (false);
}

void getLineEqn(double x1, double y1, double x2, double y2, double & m, double & c) {
    assert(x2-x1 != 0);
    if(x2-x1 !=0) m = (y2-y1)/(x2-x1);
    else m = DBL_MAX;
    c = y1 - m*x1;
}

/* Given a non-empty binary search tree, return the minimum data  
    value found in that tree. Note that the entire tree does not need
    to be searched. */
struct node * minValue(struct node* node) {
  struct node* current = node;
  
  /* loop down to find the leftmost leaf */
  while (current->left != NULL) {
    current = current->left;
  }
  return current;
}

// This function finds predecessor and successor of key in BST.
// It sets pre and suc as predecessor and successor respectively
void findPreSuc(struct node* root, struct node*& pre, struct node*& suc, Label key)
{
    // Base case
    if (root == NULL)  return ;
 
    // If key is present at root
    if (root->key == key)
    {
        // the maximum value in left subtree is predecessor
        if (root->left != NULL)
        {
            node* tmp = root->left;
            while (tmp->right)
                tmp = tmp->right;
            pre = tmp ;
        }
 
        // the minimum value in right subtree is successor
        if (root->right != NULL)
        {
            node* tmp = root->right ;
            while (tmp->left)
                tmp = tmp->left ;
            suc = tmp ;
        }
        return ;
    }
 
    // If key is smaller than root's key, go to left subtree
    if (root->key > key)
    {
        suc = root ;
        findPreSuc(root->left, pre, suc, key) ;
    }
    else // go to right subtree
    {
        pre = root ;
        findPreSuc(root->right, pre, suc, key) ;
    }
}

//---------LINE SEGMENT INTERSECTION---------------------------------------------------

// Given three colinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
bool onSegment(Point p, Point q, Point r)
{
    if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
        q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
       return true;
 
    return false;
}
 
// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(Point p, Point q, Point r)
{
    // See 10th slides from following link for derivation of the formula
    // http://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf
    int val = (q.y - p.y) * (r.x - q.x) -
              (q.x - p.x) * (r.y - q.y);
 
    if (val == 0) return 0;  // colinear
 
    return (val > 0)? 1: 2; // clock or counterclock wise
}
 
// The main function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
bool doIntersect(Point p1, Point q1, Point p2, Point q2)
{
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);
 
    // General case
    if (o1 != o2 && o3 != o4)
        return true;
 
    // Special Cases
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;
 
    // p1, q1 and p2 are colinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;
 
    // p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;
 
     // p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;
 
    return false; // Doesn't fall in any of the above cases
}

bool doIntersect(LineSegment a, LineSegment b) {
    return doIntersect(Point(a.x1, a.y1), Point(a.x2, a.y2), Point(b.x1, b.y1), Point(b.x2, b.y2));
}

double det(double a,double b,double c,double d)
{
    return a*d-b*c;
}

// returns the intersection point if the pair of line segments intersect
Point intersectp(Point p1,Point q1,Point p2,Point q2)
{
    // pricntf("Received (%0.2lf, %0.2lf) , (%0.2lf, %0.2lf) , (%0.2lf, %0.2lf) , (%0.2lf, %0.2lf) ,\n", p1.x,p1.y,q1.x,q1.y,p2.x,p2.y,q2.x,q2.y);
    Point p;
    double x1=p1.x;
    double x2=q1.x;
    double x3=p2.x;
    double x4=q2.x;
    double y1=p1.y;
    double y2=q1.y;
    double y3=p2.y;
    double y4=q2.y;
    p.x=det(det(x1,y1,x2,y2),x1-x2,det(x3,y3,x4,y4)\
        ,x3-x4)/det(x1-x2,y1-y2,x3-x4,y3-y4);
    p.y=det(det(x1,y1,x2,y2),y1-y2,det(x3,y3,x4,y4)\
        ,y3-y4)/det(x1-x2,y1-y2,x3-x4,y3-y4);
    // p.label=-1;
    return p; 
}


// returns the intersection point if the pair of line segments intersect
Point intersectp(const LineSegment & a, const LineSegment & b)
{
    /*Point p1(a.x1, a.y1);
    Point q1(a.x2, a.y2); 
    Point p2(b.x1, b.y1);
    Point q2(b.x2, b.y2);*/

    // pricntf("694 : Received (%0.2lf, %0.2lf) , (%0.2lf, %0.2lf) , (%0.2lf, %0.2lf) , (%0.2lf, %0.2lf) ,\n", p1.x,p1.y,q1.x,q1.y,p2.x,p2.y,q2.x,q2.y);
    // pricntf("694 : Received (%0.2lf, %0.2lf) , (%0.2lf, %0.2lf) , (%0.2lf, %0.2lf) , (%0.2lf, %0.2lf) ,\n", p1.x,p1.y,q1.x,q1.y,b.x1,b.y1,b.x2,b.y2);

    Point ip = intersectp(Point(a.x1, a.y1), Point(a.x2, a.y2), Point(b.x1, b.y1), Point(b.x2, b.y2));
    // Point ip = intersectp(p1,q1,p2,q2);
    // pricntf("____________________________________________________________________________________________________________Intersection of %d, %d at (%0.2lf, %0.2lf)\n", a.idx, b.idx, ip.x, ip.y);
    return ip;
}

//---------LINE SEGMENT INTERSECTION---------------------------------------------------

void swapLabels(const Label & a, const Label & b, struct node *status) {
    struct node* aptr = search(status,a);
    struct node* bptr = search(status,b);

    Label temp = aptr->key;
    aptr->key = bptr->key;
    bptr->key = temp;
}

vpt getBruteForceIntersections(vls & lines) {

    vpt bfintersections;
    for (int i = 0; i < sz(lines); ++i)
    {
        for (int j = i+1; j < sz(lines); ++j)
        {
            if(doIntersect(lines[i],lines[j])){
                Point ip = intersectp(lines[i],lines[j]);  
                ip.l1 = Label(i);
                ip.l2 = Label(j);
                bfintersections.push_back(ip);
            } 
        }
    }
    return bfintersections;
}

bool equalIntersectionPoint(Point &a, Point & b) {
    return (a.x == b.x && a.y == b.y && ((a.l1 == b.l1 && a.l2 == b.l2) || (a.l1 == b.l2 && a.l2 == b.l1) ));
}

bool isAccurate(vls & lines, vpt& intersections) {
    vpt bfintersections = getBruteForceIntersections(lines);
    sort(all(bfintersections));
    sort(all(intersections));

    /*printf(" bfintersections : \n");
    rep(i,sz(bfintersections)) cout << bfintersections[i] << " , ";
    printf("\n");

    printf(" intersections : \n");
    rep(i,sz(intersections)) cout << intersections[i] << " , ";
    printf("\n");*/

    if(sz(bfintersections) != sz(intersections)) return false;

    for (int i = 0; i < sz(intersections); ++i)
    {
        if(!equalIntersectionPoint(bfintersections[i],intersections[i])) return false;
    }
    return true;
}

bool endsInInteriors(vls & lines) {
    rep(i,sz(lines)) {
        double m = lines[i].m, c = lines[i].c;
        rep(j,sz(lines)) {
            if(i==j) continue;

            // left end
            double ml, cl;
            getLineEqn(lines[i].x1,lines[i].y1,lines[j].x1,lines[j].y1,ml,cl);
            if(ml==m && cl==c) return true;

            // right end
            double mr, cr;
            getLineEqn(lines[i].x1,lines[i].y1,lines[j].x2,lines[j].y2,mr,cr);
            if(mr==m && cr==c) return true;
        }
    }
    return false;
}


 
/* Drier program to test above function*/
int main()
{
    struct node *status = NULL;
    vpt intersections;
    spt eventQ;
    printf("Enter the number of line segments - n : ");
    int n;
    s(n);

    printf("Enter the four co-ordinates x1 y1 x2 y2 of each of the 'n' line segment in separate lines\n");

    for (int i = 0; i < n; ++i)
    {
        double x1,y1,x2,y2, m,c;
        scanf("%lf %lf",&x1, &y1);
        scanf("%lf %lf",&x2, &y2);

        getLineEqn(x1,y1,x2,y2,m,c);

        lineseg.push_back(LineSegment(x1,y1,x2,y2,m,c,i));

        if(Point(x1,y1) < Point(x2,y2)) {
            eventQ.insert(Point(x1,y1,Label(i),LEFT_END));
            eventQ.insert(Point(x2,y2,Label(i),RIGHT_END));
        }
        else {
            eventQ.insert(Point(x1,y1,Label(i),RIGHT_END));
            eventQ.insert(Point(x2,y2,Label(i),LEFT_END));
        }

        
    }

    /*int i=0;
    double x1,y1,x2,y2, m,c;
    while(scanf("%lf %lf %lf %lf",&x1, &y1, &x2, &y2) != -1) {
        getLineEqn(x1,y1,x2,y2,m,c);

        lineseg.push_back(LineSegment(x1,y1,x2,y2,m,c,i));

        if(Point(x1,y1) < Point(x2,y2)) {
            eventQ.insert(Point(x1,y1,Label(i),LEFT_END));
            eventQ.insert(Point(x2,y2,Label(i),RIGHT_END));
        }
        else {
            eventQ.insert(Point(x1,y1,Label(i),RIGHT_END));
            eventQ.insert(Point(x2,y2,Label(i),LEFT_END));
        }
        i++;
    }*/

    /*if(endsInInteriors(lineseg)) {
        printf("ENDS in INTERIORS\n");
        exit(0);
    }
    else printf("DATA IS FINE\n");*/

    
    while(!eventQ.empty()) {

        Point ev = *(eventQ.begin());
        eventQ.erase(eventQ.begin());

        // printf("-----------------------GOT FROM EVENTQ - (%lf,%lf) - Type %d\n",ev.x,ev.y,ev.type);

        

        if(ev.type == LEFT_END) {
            currx = ev.x;

            // pricntf("inserting lineseg into status - %d\n",ev.label.idx);
            // printf("BEFORE inserting %d : ",ev.label.idx );
            // inorder(status); printf("\n");
            status = insert(status, ev.label);
            // printf("AFTER inserting %d : ",ev.label.idx );
            // inorder(status); printf("\n");
            
            // pricntf("insertion COMPLETE lineseg into status - %d\n",ev.label.idx);

            
            struct node* temp = search(status,ev.label);
            // printf("INSERTION SUCCESSFULL\n");

            /* FINDING SUCCESSOR , PREDECESSOR and checking if there is an intersection */
            assert(status != NULL);
            struct node* currptr = search(status,ev.label);
            assert(currptr != NULL);

            // struct node* successor = inOrderSuccessor(status,currptr);
            struct node *pre = NULL, *suc = NULL;
            findPreSuc(status, pre, suc, currptr->key);

            // printf("Label-%d, Type-%d, ", (currptr->key).idx, ev.type);
            // if(pre != NULL) printf("Pre-%d, ", (pre->key).idx);
            // if(suc != NULL) printf("Suc-%d, ", (suc->key).idx);
            // printf("\n");
            

            // pricntf("%d %d\n",(currptr->key).idx, ev.label.idx );
            // pricntf("%d\n", (successor->key).idx );
            LineSegment a,b;
            a = lineseg[(currptr->key).idx];
            if(suc != NULL) {
                b = lineseg[(suc->key).idx];

                
                if( doIntersect(a, b) ) {
                    
                    Point intersection_point = intersectp(a,b);
                    intersection_point.type = INTERSECTION;
                    intersection_point.l1 = Label((currptr->key).idx);
                    intersection_point.l2 = Label((suc->key).idx);                    
                    if(intersection_point.x > ev.x) {
                        if(!cpresent(eventQ,intersection_point)) {
                            // printf("DETECTED AT suc LEFT_END-%d : ",ev.label.idx);
                            intersections.push_back(intersection_point);
                        }
                        eventQ.insert(intersection_point);
                    }
                }
            }

            if(pre != NULL) {
                b = lineseg[(pre->key).idx];

                
                if( doIntersect(b, a) ) {
                    
                    Point intersection_point = intersectp(b,a);
                    intersection_point.type = INTERSECTION;
                    intersection_point.l1 = Label((pre->key).idx);
                    intersection_point.l2 = Label((currptr->key).idx);
                    if(intersection_point.x > ev.x) {
                        if(!cpresent(eventQ,intersection_point)) {
                            // printf("DETECTED AT pre LEFT_END-%d : ",ev.label.idx);
                            intersections.push_back(intersection_point);
                        }
                        eventQ.insert(intersection_point);
                    }
                }
            }

            

        }
        else if(ev.type == RIGHT_END) {
            currx = ev.x;

            // delete point
            // before deleting get pre and suc and check for intersection
            // printf("searching for RIGHT_END - label %d\n",ev.label.idx );
            struct node* currptr = search(status,ev.label);
            // printf("FINISHED searching for RIGHT_END - label %d, %d\n",ev.label.idx, (currptr->key).idx );
            assert(currptr != NULL);

            struct node *pre = NULL, *suc = NULL;
            findPreSuc(status, pre, suc, currptr->key);

            LineSegment a,b;
            if(suc != NULL && pre != NULL) {
                a = lineseg[(pre->key).idx];
                b = lineseg[(suc->key).idx];

                

                if( doIntersect(a, b) ) {
                    
                    Point intersection_point = intersectp(a,b);
                    intersection_point.type = INTERSECTION;
                    intersection_point.l1 = Label((pre->key).idx);
                    intersection_point.l2 = Label((suc->key).idx);
                    if(intersection_point.x > ev.x) {
                        if(!cpresent(eventQ,intersection_point)) {
                            // printf("DETECTED AT RIGHT_END-%d : ",ev.label.idx);
                            intersections.push_back(intersection_point);
                        }
                        eventQ.insert(intersection_point);
                    }
                }
            }

            // printf("Label-%d, Type-%d, ", (currptr->key).idx, ev.type);
            // if(pre != NULL) printf("Pre-%d, ", (pre->key).idx);
            // if(suc != NULL) printf("Suc-%d, ", (suc->key).idx);
            // printf("\n");

            // printf("Trying to delete linesgment %d from status\n", ev.label.idx);
            // printf("Status Before deleting %d : ", ev.label.idx);
            // inorder(status);
            // printf("\n");
            status = deleteNode(status, ev.label);
            // printf("Status After deleting %d : ", ev.label.idx);
            // inorder(status);
            // printf("\n");

        }
        else if(ev.type == INTERSECTION) {
            currx = ev.x - h;

            // printf("Intersection event - Status : ");
            // inorder(status); printf("\n");


            // printf("searching for INTERSECTION - l1 %d , (%d,%d) | root = %d\n",ev.l1.idx, ev.l1.idx, ev.l2.idx, status->key.idx );
            struct node* l1ptr = search(status,ev.l1);
            // printf("FINISHED searching for INTERSECTION - l1 %d, %d\n",ev.l1.idx, (l1ptr->key).idx );
            assert(l1ptr != NULL);

            // printf("searching for INTERSECTION - l2 %d\n",ev.l2.idx );
            struct node* l2ptr = search(status,ev.l2);
            // printf("FINISHED searching for INTERSECTION - l2 %d, %d\n",ev.l2.idx, (l2ptr->key).idx );
            assert(l2ptr != NULL);

            
            struct node *pre_l1 = NULL, *suc_l2 = NULL, *temp = NULL;
            findPreSuc(status, pre_l1, temp, ev.l1);  // suc_l2 is dummy here
            findPreSuc(status, temp, suc_l2, ev.l2);  // pre_l1 is dummy here

            // if(pre_l1 != NULL) printf("pre_l1-%d, ", (pre_l1->key).idx);
            // if(suc_l2 != NULL) printf("suc_l2-%d, ", (suc_l2->key).idx);
            // printf("\n");

             LineSegment a,b;
            // check intersection of l1 and 'successor of l2'
            if(suc_l2 != NULL) {
                a = lineseg[ev.l1.idx];
                b = lineseg[(suc_l2->key).idx];

                


                if( doIntersect(a, b) ) {
                    
                    Point intersection_point = intersectp(a,b);
                    intersection_point.type = INTERSECTION;
                    intersection_point.l1 = ev.l1;
                    intersection_point.l2 = suc_l2->key;
                    if(intersection_point.x > ev.x) {
                        if(!cpresent(eventQ,intersection_point)) {
                            // printf("DETECTED AT intersection(%d,%d) : ",ev.l1.idx, ev.l2.idx);
                            intersections.push_back(intersection_point);
                        }
                        eventQ.insert(intersection_point);
                    }
                }
            }

            // check intersection of l1 and 'successor of l2'
            if(pre_l1 != NULL) {

                a = lineseg[(pre_l1->key).idx];
                b = lineseg[ev.l2.idx];

                if( doIntersect(a, b) ) {

                    Point intersection_point = intersectp(a,b);
                    intersection_point.type = INTERSECTION;
                    intersection_point.l1 = pre_l1->key;
                    intersection_point.l2 = ev.l2;
                    if(intersection_point.x > ev.x) {
                        if(!cpresent(eventQ,intersection_point)) {
                            // printf("DETECTED AT intersection(%d,%d) : ",ev.l1.idx, ev.l2.idx);
                            intersections.push_back(intersection_point);
                        }
                        eventQ.insert(intersection_point);
                    }
                }
            }

            // printf("Before swapping %d,%d : ",ev.l1.idx, ev.l2.idx );
            // inorder(status); printf("\n");

            swapLabels(ev.l1, ev.l2,status);

            currx = ev.x + h;
        }
    }

    sort(all(intersections));
    printf("\nResult:\n");
    printf("------\n");
    printf("Number of Intersection points = %d\n\n",sz(intersections));
    for (int i = 0; i < sz(intersections); ++i)
    {   
      cout << intersections[i] << endl;
    }


    // --------------------------------------------------------------------------------------------
    // TO output the convex hull formed out of the intersection points from the code of Assignment-1

    vmypt pts, hull;

    for (int i = 0; i < sz(intersections); ++i)
    {
        pts.push_back(myPoint(intersections[i].x, intersections[i].y));
    }

    hull = ChanHull(pts);

    printf("\nThe Convex hull formed out of the intersection points:\n");
    printf("-----------------------------------------------------\n");
    printf("Number of points in Convex Hull = %d\n",sz(hull));
    printf("\nThe list of Convex Hull points in counterclockwise order : \n{ ");
    for(int i = 0; i < sz(hull); ++i) {
      if(i!=sz(hull)-1) printf("(%0.2lf,%0.2lf), ",hull[i].x, hull[i].y);
      else printf("(%0.2lf,%0.2lf) ",hull[i].x, hull[i].y);
    }
    printf("}\n\n");

    return 0;
}