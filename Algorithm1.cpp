#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <fstream>
#include <sstream>
#include <utility>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <set>
#include <numeric>

using namespace std;

using ll = long long;
struct Edge { int to, rev; double cap; };
struct Dinic {
    int n; vector<vector<Edge>> g; vector<int> lvl, ptr;
    Dinic(int N):n(N),g(N),lvl(N),ptr(N){}
    void addEdge(int u,int v,double c){
        g[u].push_back({v,(int)g[v].size(),c});
        g[v].push_back({u,(int)g[u].size()-1,0.0});
    }
    bool bfs(int s,int t){
        fill(lvl.begin(),lvl.end(),-1);
        queue<int>q; lvl[s]=0; q.push(s);
        while(!q.empty()){
            int u=q.front(); q.pop();
            for(auto &e:g[u]){
                if(lvl[e.to]<0 && e.cap>1e-12){
                    lvl[e.to]=lvl[u]+1; q.push(e.to);
                }
            }
        }
        return lvl[t]>=0;
    }
    double dfs(int u,int t,double f){
        if(u==t||f<1e-12) return f;
        for(int &i=ptr[u];i<(int)g[u].size();++i){
            auto &e=g[u][i];
            if(lvl[e.to]!=lvl[u]+1) continue;
            double p=dfs(e.to,t,min(f,e.cap));
            if(p>1e-12){
                e.cap-=p; g[e.to][e.rev].cap+=p;
                return p;
            }
        }
        return 0;
    }
    double maxFlow(int s,int t){
        double flow=0, pushed;
        while(bfs(s,t)){
            fill(ptr.begin(),ptr.end(),0);
            while((pushed=dfs(s,t,1e18))>1e-12) flow+=pushed;
        }
        return flow;
    }
    vector<char> minCut(int s){
        vector<char> inS(n);
        queue<int>q; inS[s]=1; q.push(s);
        while(!q.empty()){
            int u=q.front(); q.pop();
            for(auto &e:g[u]){
                if(e.cap>1e-12 && !inS[e.to]){
                    inS[e.to]=1; q.push(e.to);
                }
            }
        }
        return inS;
    }
};

int main(int argc,char**argv){
    if(argc!=3){
        cerr<<"Usage: "<<argv[0]<<" <edge_list> <k>\n";
        return 1;
    }
    auto t0=chrono::high_resolution_clock::now();
    int h=stoi(argv[2]);
    if(h<1){ cerr<<"k must be >=1\n"; return 1; }

    ifstream fin(argv[1]);
    string line;
    vector<pair<int,int>> edges;
    int u,v, maxID=-1;
    while(getline(fin,line)){
        if(line.empty()||line[0]=='#') continue;
        stringstream ss(line);
        if(!(ss>>u>>v)) continue;
        if(u==v) continue;
        edges.emplace_back(u,v);
        maxID=max(maxID,max(u,v));
    }
    int n=maxID+1, m=edges.size();

    vector<vector<int>> adj(n);
    for(auto &e:edges){
        int a=e.first, b=e.second;
        adj[a].push_back(b);
        adj[b].push_back(a);
    }
    for(int i=0;i<n;++i){
        auto &v=adj[i];
        sort(v.begin(),v.end());
        v.erase(unique(v.begin(),v.end()),v.end());
    }

    vector<vector<vector<int>>> cl(h+1);
    for(int i=0;i<n;++i) cl[1].push_back({i});
    if(h>=2){
        for(int i=0;i<n;++i)
            for(int j:adj[i])
                if(i<j) cl[2].push_back({i,j});
    }
    for(int k=3;k<=h;++k){
        for(auto &sm:cl[k-1]){
            int last=sm.back();
            for(int w:adj[last]){
                if(w>last){
                    bool ok=true;
                    for(int x:sm){
                        if(x==last) continue;
                        if(!binary_search(adj[x].begin(),adj[x].end(),w)){
                            ok=false; break;
                        }
                    }
                    if(ok){
                        auto ext=sm;
                        ext.push_back(w);
                        cl[k].push_back(move(ext));
                    }
                }
            }
        }
    }

    cout<<"k="<<h<<"\n";
    cout<<"Number of vertices: "<<n<<"\n";
    cout<<"Number of edges: "<<m<<"\n";
    if(h>=2)
        cout<<"Looking for cliques of size "<<h<<" and "<<h-1<<"\n";
    else
        cout<<"Looking for cliques of size 1\n";
    if(h>=2){
        cout<<"Number of cliques of size "<<h<<": "<<cl[h].size()<<"\n";
        cout<<"Number of cliques of size "<<h-1<<": "<<cl[h-1].size()<<"\n";
    } else {
        cout<<"Number of cliques of size 1: "<<cl[1].size()<<"\n";
    }

    vector<double> cdeg(n,0.0);
    if(h==1){
        for(int i=0;i<n;++i) cdeg[i]=1.0;
    } else if(h==2){
        for(int i=0;i<n;++i) cdeg[i]=(double)adj[i].size();
    } else {
        for(auto &cm:cl[h-1])
            for(int x:cm) cdeg[x]+=1.0;
    }

    double lo=0, hi=*max_element(cdeg.begin(),cdeg.end());
    double eps=1.0/(n*max(1,n-1));
    vector<int> bestD;
    while(hi-lo>eps){
        double alpha=(lo+hi)/2;
        int S=0, T=1;
        int Voff=2;
        int Loff=2+n;
        int N=2+n+(h>1?cl[h-1].size():0);
        Dinic mf(N);
        if(h==1){
            bestD.resize(n);
            iota(bestD.begin(),bestD.end(),0);
            break;
        }
        else if(h==2){
            for(int i=0;i<n;++i){
                if(cdeg[i]>0) mf.addEdge(S,Voff+i,cdeg[i]);
                mf.addEdge(Voff+i,T,alpha*2);
            }
            for(auto &e:cl[2]){
                int a=e[0], b=e[1];
                mf.addEdge(Voff+a,Voff+b,1);
                mf.addEdge(Voff+b,Voff+a,1);
            }
        }
        else {
            for(int i=0;i<n;++i)
                if(cdeg[i]>0) mf.addEdge(S,Voff+i,cdeg[i]);
            for(int i=0;i<n;++i)
                mf.addEdge(Voff+i,T,alpha*h);
            const double INF=1e15;
            for(int i=0;i<cl[h-1].size();++i){
                for(int x:cl[h-1][i])
                    mf.addEdge(Loff+i,Voff+x,INF);
            }
            for(int v=0;v<n;++v){
                for(int i=0;i<cl[h-1].size();++i){
                    auto &cm=cl[h-1][i];
                    if(find(cm.begin(),cm.end(),v)!=cm.end()) continue;
                    bool ok=true;
                    for(int x:cm){
                        if(!binary_search(adj[x].begin(),adj[x].end(),v)){
                            ok=false; break;
                        }
                    }
                    if(ok) mf.addEdge(Voff+v,Loff+i,1);
                }
            }
        }
        mf.maxFlow(S,T);
        auto inS=mf.minCut(S);
        vector<int> curr;
        for(int i=0;i<n;++i)
            if(inS[Voff+i]) curr.push_back(i);
        if(curr.empty()) hi=alpha;
        else{ bestD=move(curr); lo=alpha; }
    }

    vector<char> inD(n,0);
    for(int x:bestD) inD[x]=1;
    ll Th=0;
    if(h==1) Th=bestD.size();
    else if(h==2){
        for(auto &e:cl[2])
            if(inD[e[0]]&&inD[e[1]]) Th++;
    } else {
        for(auto &cm:cl[h])
            if(all_of(cm.begin(),cm.end(),[&](int x){return inD[x];}))
                Th++;
    }
    double density = bestD.empty()?0.0:double(Th)/bestD.size();

    cout<<"Densest subgraph vertices:";
    for(int x:bestD) cout<<" "<<x;
    cout<<"\nNumber of vertices in densest subgraph: "<<bestD.size()<<"\n";
    cout<<"Number of h-cliques in densest subgraph: "<<Th<<"\n";
    cout<<fixed<<setprecision(4)
        <<"Density of densest subgraph: "<<density<<"\n\n";

    auto t1=chrono::high_resolution_clock::now();
    ll ms=chrono::duration_cast<chrono::milliseconds>(t1-t0).count();
    cout<<"Execution time: "<<ms<<" milliseconds\n";

    return 0;
}
