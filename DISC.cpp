#include<vector>
#include<algorithm>
#include<set>
#include"./RTree.h"

using namespace std;


typedef double ValueType;
typedef RTree<ValueType, double, 2, double> MyTree;
class DISC{
    private:
        int cluster_ID = 0;
        int cluster_check;
        MyTree* tree;
        int size;
        double eps;
        double dim;
        int minPts;
        int tick=1;
        set<Point*> result;
        set<Point> C_out;
        set<Point> excore;
        set<Point> neocore;
        set<Point> Potential_Noises;

    public:
        Union_Find uf;
        DISC(MyTree* treee,int i, double eps2, double dim2,int minPts2){
            tree = treee;
            size = i;
            eps = eps2;
            dim = dim2;
            minPts = minPts2;
            
            // MyTree::Iterator it;
            // for((*tree).GetFirst(it); !(*tree).IsNull(it);(*tree).GetNext(it)){
            //     (*it).point.label = 0;
            // } 
            
        }
        void Collect(vector<Point> in,vector<Point> out){
            set<Point*> InNeighbors;
            set<Point*> OutNeighbors;
            set<Point*>* empty;

            for(int i=0; i<out.size(); i++){
                //set<Point> temp;
                //temp.insert(out[i]);
                if(out[i].iscore == 1){
                    C_out.insert(out[i]);
                }
                else{
                    (*tree).Remove(out[i].xy, out[i].xy, out[i].DataId);
                }
                out[i].label = DELETED;
            }

            for(int i=0; i<out.size(); i++){
                set<Point*> nei = (*tree).Search_neighbors(out[i].xy, eps, dim);
                for(set<Point*>::iterator it = nei.begin();it!=nei.end();++it){
                    if((*it)->label !=DELETED){
                        (*it)->nei_count--;
                    }
                }
                out[i].nei_count = 1; //not understood
                /////////// union
                for(set<Point*>::iterator it = nei.begin();it!=nei.end();++it){
                    OutNeighbors.insert(*it);
                }
            }            

            for(int i=0; i<in.size(); i++){
                struct Rect TempRect = Rect(in[i].xy[0],in[i].xy[1],in[i].xy[0],in[i].xy[1]);
                (*tree).Insert(TempRect.min,TempRect.max,i);
            }
            
            for(int i=0; i<in.size(); i++){
                set<Point*> nei = (*tree).Search_neighbors(in[i].xy, eps, dim);
                int nc = 1;
                for(set<Point*>::iterator it = nei.begin();it!=nei.end();++it){
                    if((*it)->label !=UNCLASSIFIED && (*it)->label != DELETED ){
                        (*it)->nei_count++;
                    }
                    if((*it)->label != DELETED){
                        nc++;
                    }
                }
                out[i].nei_count = nc; //not understood
                /////////// union
                for(set<Point*>::iterator it = nei.begin();it!=nei.end();++it){
                    InNeighbors.insert(*it);
                }
            }

            for(set<Point*>::iterator it = OutNeighbors.begin();it!=OutNeighbors.end();++it){
                if((*it)->iscore && (*it)->nei_count < minPts){
                    excore.insert((**it));
                    if((*it)->label != DELETED) (*it)->label = UNCLASSIFIED;
                }
            }           
            for(set<Point*>::iterator it = InNeighbors.begin();it!=InNeighbors.end();++it){
                //cout << (*it)->iscore << endl;
                if( !(*it)->iscore && (*it)->nei_count >= minPts){
                    neocore.insert((**it));
                    (*it)->label = UNCLASSIFIED;
                }
                else if( ((*it)->iscore == 0) && ((*it)->nei_count < minPts)) (*it)->label = NOISE;
            }
            cout << "Out size : " << OutNeighbors.size() << endl;
            cout << "In size : " << InNeighbors.size() << endl;
            cout << "Excore size : " << excore.size() << endl;
            cout << "Neocore size : " << neocore.size() << endl;

        }




        void Cluster(int cluster_count){
            cluster_ID = cluster_count;
            cluster_check = cluster_count;
            for(set<Point>::iterator it = excore.begin();it!=excore.end();++it){
                if((*it).iscore){
                    set<Point> minimal_bonding_cores;
                    find_miniaml_bonding_core(minimal_bonding_cores, *it);
                    if(minimal_bonding_cores.size() > 1){
                        MS_BFS(minimal_bonding_cores);
                    }
                    cout << "minimal_bonding_cores.size : "<< minimal_bonding_cores.size() << endl;
                }
            }
            for(set<Point>::iterator it = C_out.begin();it!=C_out.end();++it){
                (*tree).Remove((*it).xy, (*it).xy, (*it).DataId);
            }

            for(set<Point>::iterator it = neocore.begin();it!=neocore.end();++it){
                if(!(*it).iscore){
                    CreateOrMergeClusters(*it);
                }
            }
            for(auto it : Potential_Noises){ //////////////////deque
                if((it).label == UNCLASSIFIED){
                    set<Point*> nei = (*tree).Search_neighbors(it.xy, eps, dim);
                    (it).label = NOISE;
                    for(auto iter : nei){
                        if(iter->iscore){
                            it.label = iter->label;
                            break;
                        }
                    }
                }
            }
        }
        void CreateOrMergeClusters(Point p){
            tick++;
            set<int> ClustersToMerge;
            deque<Point> Queue_of_neocores_for_BFS;
            deque<Point> Queue_of_neocores_for_labeling;

            Queue_of_neocores_for_BFS.push_back(p);

            set<Point> noise_nearNewCores;

            while(!Queue_of_neocores_for_BFS.empty()){
                Point x =(*Queue_of_neocores_for_BFS.begin()); //diff point
                Queue_of_neocores_for_BFS.pop_front();
                set<Point*> nei = (*tree).Search_neighbors(x.xy,eps,dim);
                for(auto iter : nei){
                    if(iter->iscore && iter->label!=UNCLASSIFIED){
                        ClustersToMerge.insert(iter->label);
                    }
                    else if(iter->nei_count>=minPts && !iter->iscore){
                        Queue_of_neocores_for_labeling.push_back(*iter);
                        Queue_of_neocores_for_BFS.push_front(*iter);
                        iter->iscore = 1;
                    }
                    else if(iter->label == NOISE || iter->label == UNCLASSIFIED){
                        noise_nearNewCores.insert(*iter);
                    }
                }
            }

            if(ClustersToMerge.size()==0){
                int NEWclusterID = cluster_ID;
                cluster_ID++;
                for(auto c:Queue_of_neocores_for_labeling){
                    c.label = NEWclusterID;
                }
                for(auto n:noise_nearNewCores){
                    n.label = NEWclusterID;
                }
                uf.create(NEWclusterID);
            }
            else if(ClustersToMerge.size()==1){
                int NEWclusterID = *(ClustersToMerge.begin());
                for(auto c:Queue_of_neocores_for_labeling){
                    c.label = NEWclusterID;
                }
                for(auto n:noise_nearNewCores){
                    n.label = NEWclusterID;
                }
            }
            else{
                int NEWclusterID = *(ClustersToMerge.begin());
                for(auto id : ClustersToMerge){
                    uf.Union(id,NEWclusterID);
                }
                for(auto c:Queue_of_neocores_for_labeling){
                    c.label = NEWclusterID;
                }
                for(auto n:noise_nearNewCores){
                    n.label = NEWclusterID;
                }
            }

        }

        void MS_BFS(set<Point> minimal_bonding_cores){
            //int OldId = minimal_bonding_cores
            int OldId = (*minimal_bonding_cores.begin()).label;
            int L = cluster_ID;
            unordered_map<int,Point> clusterID2point; // DIFF POINT
            unordered_map<int, int> point2clusterID; // DIFF POINT
            unordered_map<int,deque<Point>> id2n;

            set<Point> copied_minimal_bonding_cores;
            copy(minimal_bonding_cores.begin(),minimal_bonding_cores.end(),inserter(copied_minimal_bonding_cores, copied_minimal_bonding_cores.begin()));
            
            for(auto iter : minimal_bonding_cores){
                if(copied_minimal_bonding_cores.find(iter) != copied_minimal_bonding_cores.end()){
                    
                    uf.create(cluster_ID);                    //
                    //

                    deque<Point> temp;
                    temp.push_back(iter);
                    id2n.insert(make_pair(cluster_ID,temp));

                    //ci2p tempci2p = {cluster_ID, iter};
                    clusterID2point.insert(make_pair(cluster_ID,iter));
                    //p2ci tempp2ci = {iter,cluster_ID};
                    point2clusterID.insert(make_pair(iter.DataId,cluster_ID));
                    iter.label = cluster_ID;

                    tick++;

                    BFS(cluster_ID, L, id2n, OldId, copied_minimal_bonding_cores, clusterID2point);
                    cluster_ID++;
                }
            }
            
            while(copied_minimal_bonding_cores.size()>0){
                for(auto iter: minimal_bonding_cores){
                    if(copied_minimal_bonding_cores.find(iter) != copied_minimal_bonding_cores.end()){
                        int cid = point2clusterID[cluster_ID];
                        BFS(cid, L, id2n, OldId, copied_minimal_bonding_cores, clusterID2point);
                    }
                }
            }
        }
        int BFS(int newID, int L, unordered_map<int,deque<Point>> set_of_Queues,int oid, set<Point> minimum_bonding_cores,unordered_map<int,Point> id2seed){
            tick++;
            deque<Point> nei = set_of_Queues[uf.find(newID)];
            deque<Point> partialQueue;
            set_of_Queues.insert(make_pair(uf.find(newID), partialQueue));

            if(minimum_bonding_cores.size()==1){
                minimum_bonding_cores.clear();
                uf.Union(oid,newID);
                return 0;
            }

            while(!nei.size()==0){
                Point pp = *(nei.begin());
                nei.pop_front();

                set<Point*> A = (*tree).Search_neighbors(pp.xy,eps,dim);
                for(auto ppp : A){
                    if(ppp->nei_count>=minPts && ppp->iscore && ppp->label != UNCLASSIFIED && ppp->label != newID){
                        // visit a core point
                        if(L > ppp->label){
                            ppp->label = newID;
                            partialQueue.push_back(*ppp);

                            if(minimum_bonding_cores.erase(*ppp)){
                                if(minimum_bonding_cores.size()==1){
                                    minimum_bonding_cores.clear();
                                    uf.Union(newID,oid);
                                    return -1;
                                }
                            }
                        }
                        else if(uf.find(newID != uf.find(ppp->label))){
                            partialQueue.push_back(*ppp);
                        
                            for(auto mbc : minimum_bonding_cores){
                                if(uf.find(ppp->label) == uf.find(mbc.label)){
                                    if(minimum_bonding_cores.erase(mbc)){
                                        break;
                                    }
                                }
                            }

                            for(auto addall:set_of_Queues[uf.find(ppp->label)]){ // diff point
                                partialQueue.push_back(addall);
                            }
                            uf.Union(newID,ppp->label);
                            set_of_Queues[uf.find(newID)] = partialQueue;

                            if(minimum_bonding_cores.size()==1){
                                minimum_bonding_cores.clear();
                                uf.Union(oid,newID);
                                return 0;
                            }
                        }
                    }
                    else if(ppp->nei_count < minPts && ppp->iscore && ppp->label != newID){
                        partialQueue.push_back(*ppp);
                        ppp->label = newID;
                    }
                    else if(ppp->nei_count < minPts && ppp->label < cluster_check){
                        if(ppp->label!=DELETED){
                            Potential_Noises.erase(*ppp);
                            ppp->label = newID;
                        }
                    }
                }
            }
            if(partialQueue.size()==0){
                minimum_bonding_cores.erase(id2seed[newID]);
                return 0;
            }
            return 0;
            
        }
        void find_miniaml_bonding_core(set<Point> minimal_bonding_cores, Point p){
            tick++;
            set<Point*> nei = (*tree).Search_neighbors(p.xy, eps, dim);
            vector<Point*> nei_v(nei.size());
            copy(nei.begin(), nei.end(), nei_v.begin());

            for(auto iter:nei_v){ // diff point
                if(iter->nei_count >= minPts && iter->iscore ){ // if core
                    minimal_bonding_cores.insert(*iter);
                }
                else if(iter->nei_count < minPts && iter->iscore){ // if excore 
                    set<Point*> nei2 = (*tree).Search_neighbors(iter->xy, eps, dim);
                    vector<Point*> nei_v2(nei.size());
                    copy(nei2.begin(), nei2.end(), nei_v2.begin());
                    
                    iter->iscore = 0;
                    if(iter->label != DELETED){
                        Potential_Noises.insert(*iter);
                        iter->label = UNCLASSIFIED;
                    }

                    for(auto iter2 : nei_v2){
                        if(iter2->nei_count >= minPts && iter2->iscore){minimal_bonding_cores.insert(*iter2);}
                        else if(iter2->nei_count < minPts && iter2->iscore){nei_v.push_back(iter2);} // infinite loop ?
                        else if(iter2->nei_count < minPts && !iter2->iscore && iter2->label!=UNCLASSIFIED && iter2->label!=DELETED && iter2->label < cluster_check){ // diff point
                            minimal_bonding_cores.insert(*iter2);
                            iter2->label= UNCLASSIFIED;
                        }
                        
                    }

                }
                else if(iter->nei_count<minPts && !iter->iscore && iter->label!=UNCLASSIFIED&&iter->label!=DELETED&&iter->label < cluster_check){ // diff point
                    Potential_Noises.insert(*iter);
                    iter->label = UNCLASSIFIED;
                }
            }
        }

};


