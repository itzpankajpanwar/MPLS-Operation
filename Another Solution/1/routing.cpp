//#include<bits/stdc++.h>
#include<cstdio>
#include<iostream>
#include<fstream>
#include<vector>
#include<set>
#include<queue>
#include<cstring>

using namespace std;

typedef pair<int, int> ii;
typedef vector<int> vi;
typedef pair<ii, ii> i4;

#define INF 1000000000
#define pb push_back
#define mp make_pair

#define x first
#define y second

#define MAXPATHSIZE 10


struct conr{
	int u, v, id;
	float equi;
	bool admit;
	int path_n;
	conr(int _u, int _v, int _id, float _equi) : u(_u), v(_v), id(_id), equi(_equi), admit(false), path_n(-1){}
};

struct node{
	int v, delay;
	float speed, metric, used;
	node(int _v, int _delay, float _speed, float _metric) : v(_v), delay(_delay), speed(_speed), metric(_metric), used(0) {}
};

struct rout_node{
	int dest, path_cost, path_delay;
	vi path;
	rout_node(){ }
};

struct forward_node{
	int rout_id, inc_port, inc_vcid, out_port, out_vcid;
	forward_node(int a, int b, int c, int d, int e) : rout_id(a), inc_port(b), inc_vcid(c), out_port(d), out_vcid(e) { }
};

struct path_node{
	int id, src, dest, path_cost;
	bool admitted;
	vi path, vcid_list;
	path_node(int _id, int _src, int _dest, int _path_cost) : id(_id), src(_src), dest(_dest), path_cost(_path_cost), admitted(false){}
};

vector<vector<node>> nw;
vi par;
vi dist;
vector<conr> conns;					//conns vector store all connection requests.
set<ii> edgeign;						//Edges to ignore while considering second shortest path.

vector<vector<vector<i4>>> path[2];			//2 shortest paths between all pairs.
								//Path are stored in ((u, iv), (v, iu)) format
								//Where iv is index of v node in the adjacency list of u.
								//iu is defined similarly.

//Output vectors.
vector<vector<pair<rout_node, rout_node>>> rout_table;			//For two shortest paths.
vector<vector<forward_node>> forward_table;					//Contains all the information for forwarding table
vector<path_node> path_table;								//Contains all the information for path table

char *topfile, *confile, *routfile, *forwdfile, *pathfile;
int hop=0, pess=0;
int n, m, commcount;

void parse(int argc, char *argv[]);
void create_graph();
void print_graph();
void con_req();
void dijkstras(int src, int dest);
void _2spath();
void print_path(int src, int dest, int pn);
void pathsd(int src, int v, int dest, int pi);	//pathsd saves edges to ignore, pi==0 shortest, 1--> 2nd shortest
								//also saves adds path to indices.
void process_reqs();
void route();
void forward();
void save_out();
void out();


int main(int argc, char *argv[])
{		
	if(argc < 2)
		return 0;
	topfile = argv[1];
	confile = NULL;
	commcount = 0;
	
	parse(argc, argv);
	con_req();
	create_graph();
	
	
	path[0].assign(n, vector<vector<i4>>(n, vector<i4>()));
	path[1].assign(n, vector<vector<i4>>(n, vector<i4>()));
	
	_2spath();
	process_reqs();
	route();
	forward();
	save_out();
	out();
}

//Printing output on the screen.
void out()
{
	ifstream ff("pf.txt");
	string str;
	while(getline(ff, str))
		cout << str << "\n";
}

//This function stores two shortest paths and related info between every node pair in the rout_table vector
void route()
{
	int i, j, k;
	int cost, delay;
	rout_table.assign(n, vector<pair<rout_node, rout_node>>());
	for(i=0; i<n; i++)
	{
		for(j=0; j<n; j++)
		{
			if(i==j)
				continue;
			rout_node rn[2];
			for(k=0; k<2; k++)
			{
				rn[k].dest = j;
				cost = delay = 0;					//Currently cost is hop
				cost = path[k][i][j].size();
				rn[k].path_cost = cost;
				rn[k].path_delay = delay;
				if(cost==0)
					continue;
				for(auto p: path[k][i][j])
				{
					delay += nw[p.x.x][p.x.y].delay;
					rn[k].path.pb(p.x.x);
				}
				rn[k].path_delay = delay;
				rn[k].path.pb(j);
			}
			rout_table[i].pb(mp(rn[0], rn[1]));
			
		}
	}
}

//This function assigns vcid for each accepted conns
//It stores forwading table for each node
//It also stores Path table for stdout.
void forward()
{
	vi invcid(n, 0), outvcid(n, 0);
	forward_table.assign(n, vector<forward_node>());
	int i, j, cur, u, v;
	
	for(int l=0; l<(int)conns.size(); l++)
	{
		auto p = conns[l];
		
		if(!p.admit || p.path_n==-1)
			continue;
		
		auto &patht = path_table[l];
		auto &paath = path[p.path_n][p.u][p.v];
		patht.admitted = true;
		vector<ii> pairs;					//Incoming port vci
		for(i=0; i<(int)paath.size(); i++)
		{
			auto k = paath[i].x;
			u = k.x;
			v = nw[u][k.y].v;
			invcid[v] = outvcid[u] = cur = max(invcid[v], outvcid[u]) + 1;
			if(i==0)
				pairs.pb(mp(-1, -1));
			pairs.pb(mp(k.y+1, cur));
			pairs.pb(mp(paath[i].y.y+1, cur));
			patht.path.pb(k.x);
			if(i==(int)paath.size()-1)
			{
				patht.path.pb(paath[i].y.x);
				pairs.pb(mp(-1, -1));
			}
			patht.vcid_list.pb(cur);
		}
		patht.path_cost = paath.size();
		j=0;
		
		for(i=0; i<(int)paath.size(); i++)
		{
			u = paath[i].x.x;
			forward_table[u].pb(forward_node(u, pairs[j].x, pairs[j].y, pairs[j+1].x, pairs[j+1].y));
			j+=2;
			
			if(i==(int)paath.size()-1)
			{
				
				u = paath[i].y.x;
				forward_table[u].pb(forward_node(u, pairs[j].x, pairs[j].y, pairs[j+1].x, pairs[j+1].y));
			}
		}
	}
}

//This function takes requests one by one
//If request can be admitted it assigns particular bandwidth to it.
void process_reqs()
{
	int i;
	int count=0;
	for(auto &p:conns)						//Each requests
	{
		int u=p.u;
		int v = p.v;
		for(i=0; i<2; i++)				//2 shortest paths
		{
			float mn = INF;
			if(path[i][u][v].size()==0)
				continue;
			for(auto k:path[i][u][v])
			{
				auto &nd = nw[k.x.x][k.x.y];
				mn = min(mn, nd.speed - nd.used);
			}
			if(mn < p.equi)
				continue;
			p.admit = true;
			//Subtract values.
			for(auto k:path[i][u][v])
			{
				auto &nd = nw[k.x.x][k.x.y];
				auto &nd2 = nw[k.y.x][k.y.y];
				nd.used += p.equi;
				nd2.used += p.equi;
			}
			p.path_n = i;
			break;
		}
		if(p.admit){
			count++;
			//printf("Request (%d, %d) id %d is admitted.\n", u, v, p.id);
		}
	}
}

//This function calls dijkstras twice
//second time edgeign contains some edges that are path of first shortest path
//Hence link disjoint.
// O(V^3. E)	//Worst case
void _2spath()
{
	int i, j;
	
	for(i=0; i<n; i++)
	{
		for(j=i+1; j<n; j++)
		{
			dijkstras(i, j);
			
			if(dist[j] != INF)
				pathsd(i, j, j, 0);
			dijkstras(i, j);
			
			if(dist[j] != INF)
				pathsd(i, j, j, 1);
			edgeign.clear();
		}
	}
	//Since the above loop calculates path for (i, j) where i<j
	//Remaining paths are just reverse of their respective entry.
	for(i=0; i<n; i++)
		for(j=0; j<i; j++)
			for(int pn=0; pn<2; pn++)
				for(auto k=path[pn][j][i].rbegin(); k!= path[pn][j][i].rend(); k++)
					path[pn][i][j].pb(mp((*k).y, (*k).x));
}

//This function saves output to respective files
void save_out()
{
	//Saving Routing Table.
	ofstream fr(routfile);
	if(fr)
	{
		for(int i=0; i<n; i++)
		{
			fr << "\n\nRouting Table for Node " << i << "\n\n";
			fr << "Destination\tPath";
			for(int kk=0; kk<MAXPATHSIZE; kk++)
				fr << "\t";
			fr << "Path Delay\tPath Cost\n";
			for(auto rn : rout_table[i])
			{
				for(int k=0; k<2; k++)
				{
					auto rnk = rn.x;
					if(k==1)	rnk = rn.y;
					fr << rnk.dest << "\t\t\t";
					for(auto p:rnk.path)
						fr << p << "\t";
					for(int kk=0; kk<MAXPATHSIZE-(int)rnk.path.size(); kk++)
						fr << "\t";
					fr << rnk.path_delay << "\t\t\t" << rnk.path_cost << "\n";
				}
			}
		}
		fr.close();
	}
	//Saving Forwarding Table
	ofstream ff(forwdfile);
	if(ff)
	{
		for(int i=0; i<n; i++)
		{
			ff << "\n\nForwarding Table for Node " << i << "\n";
			ff << "Router's ID\tID of Incoming Port\t\tVCID\t\tID of Outgoing Port\t\tVCID\n";
			for(auto p:forward_table[i])
			{
				ff << p.rout_id << "\t\t\t" << p.inc_port << "\t\t\t\t\t" << p.inc_vcid << "\t\t" << p.out_port << "\t\t\t\t\t" << p.out_vcid << "\n";
			}
		}
		ff.close();
	}
	//Save Path file
	ofstream fpp(pathfile);
	if(fpp)
	{
		int count = 0;
		fpp << conns.size() << " ";
		for(auto p:conns)
			if(p.admit)
				count++;
		fpp << count << "\n";
		fpp.close();
	}
	ofstream fp("pf.txt");
	if(fp)
	{
		fp << "Conn ID\tSource\tDestination\tPath";
		for(int kk=0; kk<MAXPATHSIZE-2; kk++)
			fp << "\t";
		fp << "VCID List";
		for(int kk=0; kk<MAXPATHSIZE-2; kk++)
			fp << "\t";
		fp << "Path Cost\n\n";
		for(auto p:path_table)
		{
			if(!p.admitted)
				continue;
			fp << p.id << "\t" << p.src << "\t" << p.dest << "\t\t";
			for(auto k:p.path)
				fp << k << "\t";
			for(int kk=0; kk<MAXPATHSIZE-(int)p.path.size()-2; kk++)
				fp << "\t";
			for(auto k:p.vcid_list)
				fp << k << "\t";
			for(int kk=0; kk<MAXPATHSIZE-(int)p.vcid_list.size()-2; kk++)
				fp << "\t";
			fp << p.path_cost << "\n";
		}
		fp.close();
	}
}

/* It is a recursive function that is called for every node pair twice.
 * After the recursive part, u, v appears in order for the required path.
 * And hence we're abe to store the path in a correct sequence.
 * O(V.E)
 */
 
void pathsd(int src, int v, int dest, int pi)
{
	if(par[v] == -1)	return;
	edgeign.insert(mp(v, par[v]));
	edgeign.insert(mp(par[v], v));
	
	if(par[v] != src)
	{
		pathsd(src, par[v], dest, pi);
	}
	
	int u = par[v], tov=-1, tou=-1;
	//Finding v in u.
	for(int i=0; i<(int)nw[u].size(); i++)
		if(v==nw[u][i].v)
		{
			tov = i;
			break;
		}
	//Finding u in v.
	for(int i=0; i<(int)nw[v].size(); i++)
		if(u==nw[v][i].v)
		{
			tou = i;
			break;
		}
	path[pi][src][dest].pb(mp(mp(u, tov), mp(v, tou)));
	
}

// Used to calculate shortest paths.
// Called twice for each node pair.
// Second time it ignores some edges by checking edgeign.
// O((V+E)(logv)^2)
void dijkstras(int src, int dest)
{
	par.clear();
	par.assign(n, -1);
	dist.clear();
	dist.assign(n, INF);
	dist[src] = 0;
	priority_queue<ii, vector<ii>, greater<ii>> pq;
	pq.push(mp(0, src));
	int d, u;
	while(!pq.empty())
	{
		ii front = pq.top();
		pq.pop();
		d = front.x;
		u = front.y;
		if(d > dist[u])
			continue;
		for(auto to:nw[u])
		{
			if(edgeign.find(mp(u, to.v)) != edgeign.end())
				continue;
			if(dist[u] + to.metric < dist[to.v])
			{
				par[to.v] = u;
				dist[to.v] = dist[u] + to.metric;
				pq.push(mp(dist[to.v], to.v));
			}
		}
	}
	edgeign.clear();
}

// This function parses the cmd argument
void parse2(char *comma, char *arga)
{
	string comm(comma), arg(arga);
	if(comm=="-top")
	{
		topfile = arga;
		commcount |= (1<<0);
		return;
	}
	else if(comm=="-conn")
	{
		confile = arga;
		commcount |= (1<<1);
		return;
	}
	else if(comm=="-rt")
	{
		routfile = arga;
		commcount |= (1<<2);
	}
	else if(comm=="-ft")
	{
		forwdfile = arga;
		commcount |= (1<<3);
	}
	else if(comm=="-path")
	{
		pathfile = arga;
		commcount |= (1<<4);
	}
	else if(comm=="-flag")
	{
		if(arg=="dist")
		{
			hop = 0;
			commcount |= (1<<5);
		}
		else if(arg=="hop")
		{
			hop = 1;
			commcount |= (1<<5);
		}
	}
	else if(comm=="-p")
	{
		if(arg=="0")
		{
			pess = 0;
			commcount |=(1<<6);
		}
		else if(arg=="1")
		{
			pess = 1;
			commcount |= (1<<6);
		}
	}
	else
		cout << comm << arg << endl;
}

void parse(int argc, char *argv[])
{
	if(argc != 15)
	{
		printf("input error.\n");
		exit(0);
	}
	for(int i=1; i<14; i++, i++)
	{
		parse2(argv[i], argv[i+1]);
	}
	if(commcount != (1<<7) - 1)
	{
		printf("Wrong arguments.\n");
		printf("comm count is %d.\n", commcount);
		exit(0);
	}
}

// Function for printing the graph
void print_graph()
{
	int i;
	printf("\nGraph is\n");
	for(i=0; i<n; i++)
	{
		printf("%d:\t", i);
		for(auto v:nw[i])
			printf("(%d, (%d, %f), %f)\t\t", v.v, v.delay, v.speed, v.used);
		printf("\n\n");
	}
}

// Function which takes input of conn requests.
void con_req()
{
	int con, u, v;
	float cost, mn, mx, ave;
	ofstream log("log");
	ifstream conin(confile);
	conin >> con;
	log << "conn size is " << con << endl;
	for(int i=1; i<=con; i++)
	{
		conin >> u >> v >> mn >> ave >> mx;
		((pess) ? (cost = mx) : (cost = (min(mx, (float)(ave + 0.25*(mx-mn))))));
		
		conns.pb(conr(u, v, i, cost));
		path_table.pb(path_node(i, u, v, 0));
	}
	conin.close();
	log.close();
}

//Function which takes top and saves it as graph.
void create_graph()
{
	int u, v, delay;
	float rel, speed, metric;
	ifstream topin(topfile);
	topin >> n >> m;
	
	metric = 1;
	nw.assign(n, vector<node>());
	while(m--)
	{
		topin >> u >> v >> delay >> speed >> rel;
		if(!hop)
			metric = delay;
		nw[u].pb(node(v, delay, speed, metric));
		nw[v].pb(node(u, delay, speed, metric));
		
	}
	topin.close();
}

//Function to print the path.
void print_path(int src, int dest, int pn)
{
	printf("\n%d -- %d\n", src, dest);
	
	int u, v;
	if(pn!=0 && pn!=1)
		return;
	u = src;
	v = dest;
	
	for(auto p:path[pn][u][v])
	{
		printf("[(%d, %d), (%d, %d)]\t", p.x.x, nw[p.x.x][p.x.y].v, p.y.x, nw[p.y.x][p.y.y].v);
	}
	
	
	printf("\n");
}
