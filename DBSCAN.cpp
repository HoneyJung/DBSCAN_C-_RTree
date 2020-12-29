#include <iostream>
#include<list>
#include<algorithm>
#include<vector>
#include"./RTree.h"
using namespace std;
typedef double ValueType;
typedef RTree<ValueType, double, 2, double> MyTree;
class DBSCAN{
    private:
        int cluster_count = 0;
        MyTree* tree;
        int size;
        double eps;
        double dim;
        int minPts;
        vector<Point*> seed;
        vector<Point*> seed2;

    public:
        DBSCAN(MyTree* treee,int i, double eps2, double dim2,int minPts2){
            tree = treee;
            size = i;
            eps = eps2;
            dim = dim2;
            minPts = minPts2;
        }
        vector<Point*> Set_Union(vector<Point*> seed, vector<Point*> seed2){
                vector<Point*> new_seed;
                new_seed.resize(seed.size()+seed2.size());
                auto itr = set_union(seed.begin(),seed.end(),seed2.begin(),seed2.end(),new_seed.begin());
                new_seed.erase(itr,new_seed.end());
                //cout << new_seed.size() << endl;
                //cout << seed.size() << endl;
                //cout << seed2.size() << endl;
                return new_seed;
        }
        void Cluster(){
            double boundsMin[2] = {0,0};
            double boundsMax[2] = {0,0};
            MyTree::Iterator it;
            cout << "Data size : " << size << endl;
            for((*tree).GetFirst(it); !(*tree).IsNull(it);(*tree).GetNext(it))
            {
                if((*it).point.label != 0){
                    continue;
                }
                it.GetBounds(boundsMin,boundsMax);
                seed = (*tree).Search_neighbors(boundsMin, eps, dim);
                //cout << "label " << (*it).point.label << endl;
                //cout << boundsMin[0] << ' ' << boundsMin[1] << endl;
                //cout << eps << ' ' << dim << endl;
                //cout << "neighbor size " << seed.size() << endl;
                if (seed.size() < minPts){
                     (*it).point.label = -1; // noise
                     continue; 
                }
                cluster_count++;
                (*it).point.label = cluster_count; // noise

                //TODO set
                (*tree).GetNext(it);
                //cout << (*it).point.DataId << seed[0]->DataId << endl;
                vector<Point*>::iterator iter;;
                //for(iter =seed.begin(); iter != seed.end();++iter){
                for(int i=0; i < seed.size(); i++){ 
                    if( seed[i]->DataId == (*it).point.DataId){
                        continue;
                    }
                    if(seed[i]->label == -1){
                        seed[i]->label = cluster_count;
                    }
                    if(seed[i]->label != 0){
                        continue;
                    }
                    seed[i]->label = cluster_count;
                    seed2 = (*tree).Search_neighbors(seed[i]->xy, eps, dim);
                    if(seed2.size() >= minPts){
                        seed = Set_Union(seed,seed2);
                    }
                }
                
            }
            cout << cluster_count << endl;

        }

        

};


