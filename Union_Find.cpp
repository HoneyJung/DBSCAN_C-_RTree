#include"./RTree.h"
using namespace std;
class Union_Find{
    public:
        //map<int,int> id;
        unordered_map<int,int> id;
        unordered_map<int,int> sz;

        Union_Find(){
        }
        void create(int i){
            id.insert(make_pair(i,i));
            sz.insert(make_pair(i,i));
        }
        void Union(int p ,int q){
            int pRoot = find(p);
            int qRoot = find(q);

            if(pRoot == qRoot) return;
            if(sz[pRoot] <= sz[qRoot]){
                id.insert(make_pair(pRoot,qRoot));
                sz.insert(make_pair(qRoot,sz[qRoot] + sz[pRoot]));
            }
            else{
                id.insert(make_pair(qRoot,pRoot));
                sz.insert(make_pair(pRoot,sz[pRoot] + sz[qRoot]));
            }
        }
        int find(int i){
            if(id.find(i) != id.end()){
                if(id[i] == i){
                    return i;
                }
                int parent = find(id[i]);
                id.insert(make_pair(i,parent));
                return parent;
            }
            else{
                return -1;
            }
        }
        bool connected(int p,int q){
            return find(p) == find(q);
        }
        int size(){
            return id.size();
        }
};


