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
    if(mu <= lambda){
        cout<<"System will not be in stable state"<<endl;
        return 0;
    }
    double pm,em;
    pm = 3600.0/lambda;
    em = 3600.0/mu;
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generatorp (seed), generatore (seed);
    poisson_distribution<long int> poisson (pm);
    queue<pair<double, double>> Queue; //arrival_time , service_time
    double wait_time = 0.0, last_arrival = 0.0, service_time = 0.0;
    double last_finished[3] = {0.0,0.0,0.0};
    double last_started[3] = {0.0,0.0,0.0};
    long int waiting_in_queue = 0;
    long int being_inspected = 0;
    long int waiting;
    vector<double> arrival;
    vector<double>::iterator itr;
    Queue.push(make_pair(0, exponential(em)));
    arrival.push_back(0);
    for(long int i = 1; i < population; i++){
        double at, st;
        at = poisson(generatorp);
        st = exponential(em);
        at += last_arrival; 
        Queue.push(make_pair(at, st));
        arrival.push_back(at);
        last_arrival = at;
    } 
    int index = 0;
    while(true){
        waiting = 0;
        itr = arrival.begin();
        double avg;
        avg = (last_finished[0] + last_finished[1] + last_finished[2])/3;
        while(itr != arrival.end() && *itr < avg){
            waiting++;
            itr++;
        }
        waiting_in_queue += waiting;
        double wt, at, st;
        if(!Queue.empty()){
            int count = 0;
            for(int i = 0; i < 2; i++){
                if(at <= last_finished[i]){
                    count++;
                }
            }
            being_inspected += count; 
            at = Queue.front().first;
            st = Queue.front().second; 
            Queue.pop();
            service_time += st;
            count = 0;
            if(last_finished[index] - at > 0){
                wt = last_finished[index] - at;
            }else{
                wt = 0;
            }
            wait_time += wt;
            last_started[index] = max(last_finished[index], at);
            last_finished[index] = max(last_finished[index], at) + st;
            arrival.erase(arrival.begin());     
            for(int i = 0; i < 3; i++){
                if(at <= last_finished[i]){
                    count++;
                }
            }
            being_inspected += count; 
            index++;
            index %= 3;
        }else{
            break;
        }
    }
    double avg_psgr, avg_wt, avg_resp, avg_psgr_in_queue;
    avg_wt = wait_time/population;
    avg_resp = avg_wt + service_time/population;
    avg_psgr_in_queue = waiting_in_queue/(3*population);
    avg_psgr = being_inspected/(3*population);
    cout<<"Average number of passengers getting inspected: "<<avg_psgr<<endl;
    cout<<"Average response time for passengers in getting inspected: "<<avg_resp<<" seconds."<<endl;
    cout<<"Average time for which a passenger has to wait until getting inspected: "<<avg_wt<<" seconds."<<endl;
    cout<<"Average number of passengers waiting in queue before each officer(except the one getting inspected): "<<avg_psgr_in_queue<<endl;
    return 0;
}

