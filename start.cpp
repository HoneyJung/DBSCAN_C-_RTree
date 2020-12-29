#include <iostream>
#include <string.h>
#include <fstream>
//#include "./disc/DBSCAN.cpp"
#include "./DBSCAN.cpp"
#include "./RTree.h"
#include <list>
#include <vector>

using namespace std;
typedef double ValueType;
struct Rect
{
  Rect()  {}

  Rect(double a_minX, double a_minY, double a_maxX, double a_maxY)
  {
    min[0] = a_minX;
    min[1] = a_minY;

    max[0] = a_maxX;
    max[1] = a_maxY;
  }


  double min[2];
  double max[2];
};

Rect search_rect(6, 4, 10, 6); // search will find above rects that this one overlaps



int main(int argc, char** argv) {    
    string v1,v2,v3,v4;         
    fstream fs; 
    int Size = 10000;
    list<double> x;
    vector<Point> dataset;

    typedef RTree<ValueType, double, 2, double> MyTree;
    MyTree tree;
    double exampletarget[2] = {42,41};
    double eps = 5;
    double dim = 2;
    int minPts = 5;
    int i;
    vector<Point*> nei;
    try
    {   fs.open("./dataset/Maze.csv",ios::in);
        i = 0;
        while(i<Size){
            getline(fs,v1,',');
            getline(fs,v2,',');
            getline(fs,v3,',');
            getline(fs,v4,'\n');
            struct Rect TempRect = Rect(stod(v1),stod(v2),stod(v1),stod(v2));
            tree.Insert(TempRect.min,TempRect.max,i);
            i++;
        }
        fs.close();
    }
    catch(...)
    {
        cout << "file read error" << endl;
    }
    //int nhits = tree.Search(search_rect.min, search_rect.max);
    nei = tree.Search_neighbors(exampletarget, eps, dim);
    //cout << "nei size:  " << nei.size() << endl;
    //cout << "eihrie  " << nei[0]->label << endl;
    //cout << "tree.count " << tree.Count() << endl;

    DBSCAN dbscan = DBSCAN(&tree,i,eps,dim, minPts);
    dbscan.Cluster();
    MyTree::Iterator it;
    for((tree).GetFirst(it); !(tree).IsNull(it);(tree).GetNext(it)){
      if ((*it).point.label == 1){
        cout << (*it).point.xy[0] << ' ' << (*it).point.xy[0] << ' ' << (*it).point.label << endl;
      }  
    }
}
