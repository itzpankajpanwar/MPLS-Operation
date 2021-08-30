#include<bits/stdc++.h>
using namespace std;
#define population 100000.0

double exponential(double x)
{
	double z;                     
	do {
		z = ((double) rand() / RAND_MAX);
	}
	while ((z == 0) || (z == 1));
	return(-x * log(z));
}

int main(){
    double lambda, mu, hours;
    cout<<"Enter arrival rate in passengers per hour"<<endl;
    cin>>lambda;
    cout<<"Enter service rate in passengers per hour"<<endl;
    cin>>mu;
    /*if(mu <= lambda){
        cout<<"System will not be in stable state"<<endl;
        return 0;
    }*/
    double pm,em;
    pm = 3600.0/lambda;
    em = 3600.0/mu;
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generatorp (seed), generatore (seed);
    poisson_distribution<long int> poisson (pm);
    queue<pair<double, double>> Queue[3]; //arrival_time , service_time
    double wait_time[3] = {0.0,0.0,0.0};
    double last_arrival[3] = {0.0,0.0,0.0};
    double last_finished[3] = {0.0,0.0,0.0};
    double service_time[3] = {0.0,0.0,0.0};
    long int passengers_in_queue[3] = {0, 0, 0};
    long int being_inspected[3] = {0, 0, 0};
    long int count[3] = {0, 0, 0};
    vector<long int> arrival[3];
    vector<long int>::iterator itr;
    for(long int m = 0; m < 3; m++){
        double number = exponential(em);
        Queue[m].push(make_pair(0, number));
        arrival[m].push_back(0);
        service_time[m] += number;
    }
    long int k = 0, p = 0;
    long int y = 3*population;
    while(y){
        double wt, at, st;
        for(int j = 0; j < 3; j++){  
            for(int i = 0; i < 3; i++){
                if(!Queue[i].empty() && Queue[i].front().first <= last_finished[i]){
                    being_inspected[i]++;
                }
            }
            if(!Queue[j].empty()){
                at = Queue[j].front().first;
                st = Queue[j].front().second; 
                service_time[j] += st;
                if(at <= last_finished[j]){
                    wt = last_finished[j] - at;
                }else{
                    wt = 0;
                }
                wait_time[j] += wt;
                last_finished[j] = max(last_finished[j], at) + st;
                Queue[j].pop();
                arrival[j].erase(arrival[j].begin());
                for(int i = 0; i < 3; i++){
                    if(!Queue[i].empty() && Queue[i].front().first <= last_finished[i]){
                        being_inspected[i]++;
                    }
                }
            }
        }
        for(int n = 0; n < 3; n++){
            count[n] = 0;
            itr = arrival[n].begin();
            while(itr != arrival[n].end() && *itr < last_finished[n]){
                count[n]++;
                itr++;
            }
            passengers_in_queue[n] += count[n];
        }
        while(p <= 3*population){
            double a, s;
            int x = 0;
            for(long int j = 0; j < 3; j++){
                a = poisson(generatorp);
                s = exponential(em);
                a += last_arrival[j]; 
                service_time[j] += s;
                if(a <= last_finished[j]){
                    if(Queue[j].size() < 10){
                        Queue[j].push(make_pair(a, s));
                        arrival[j].push_back(a);
                        last_arrival[j] = a;
                        p++;
                    }
                    else{
                        p++;
                    }
                }else{
                    x++;
                }
            }
            if(x == 3){
                break;
            }
        }
        y-=3;
        if(y <= 0){
            break;
        }
    }
    double avg_psgr, avg_wt, avg_resp, avg_psgr_in_queue;
    double bi = 0, wt = 0, pas = 0, st = 0;
    for(long int j = 0; j < 3; j++){
        wt += wait_time[j];
        pas += passengers_in_queue[j];
        bi += being_inspected[j];
        st += service_time[j];
    } 
    avg_wt = wt/(population);
    avg_resp = avg_wt + st/(3*population);
    avg_psgr_in_queue = pas/(population);
    avg_psgr = bi/(population*3);
    cout<<"Average number of passengers getting inspected: "<<avg_psgr<<endl;
    cout<<"Average response time for passengers in getting inspected: "<<avg_resp<<" seconds."<<endl;
    cout<<"Average time for which a passenger has to wait until getting inspected: "<<avg_wt<<" seconds."<<endl;
    cout<<"Average number of passengers waiting in queue before each officer(except the one getting inspected): "<<avg_psgr_in_queue<<endl;
    return 0;
}

