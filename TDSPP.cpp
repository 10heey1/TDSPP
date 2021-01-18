// TDSPP.cpp : This file contains the 'main' function. Program execution begins and ends there.
// TO DO: add more breakpoints associated with minimums (not just THE minimum).


#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <fstream>
#include <chrono> 
#include <sstream>
#include <algorithm>
#include <queue>
#include <set>
#include <stack>
#include <assert.h>
using namespace std;
typedef pair<int,int> type_arc;
typedef pair<int, int> type_bp;

struct Output {
    bool md = 0;
    bool ddd = 0;
    int bpexplored = 0;
    double runtime = 0;
    double addruntime = 0;
    double spruntime = 0;
    int arctotal = 0;
    int subpathtotal = 1;
    double optval = 0;
    int iter = 0;
    void writeOutputCSV(string& filename) {
        ofstream myFile(filename, ofstream::app);
        myFile << this->md << ',';
        myFile << this->ddd << ',';
        myFile << this->bpexplored << ',';
        myFile << this->runtime << ',';
        myFile << this->addruntime << ',';
        myFile << this->spruntime << ',';
        myFile << this->arctotal << ',';
        myFile << this->subpathtotal << ',';
        myFile << this->iter << ',';
        myFile << this->optval << endl;
        myFile.close();
        return;
    }
};
class SummaryOutput {
public:
    bool md;
    int n;
    int T;
    int gtype;
    int ttype;
    int numseed;
    int bpexplored = 0;
    double runtime = 0;
    double enumtime = 0;
    int arctotal = 0;
    int subpathtotal = 0;
    double avgbpexplored = 0;
    int bptotal = 0;
    double percentbp = 0;
    double avgruntime = 0;
    double avgenumtime = 0;
    double percenttime = 0;
    double gavgpercenttime = 1;
    double avgarc = 0;
    double avgsubpath = 0;
    double totaladdruntime = 0;
    double totalspruntime = 0;
    double avgaddruntime = 0;
    double avgspruntime = 0;
    int iter = 0;
    string filename;
    SummaryOutput(string filename, int n, int T, int gtype, int ttype, int numseed, bool md){
        this->filename = filename;
        this->n = n;
        this->T = T;
        this->gtype = gtype;
        this->ttype = ttype;
        this->numseed = numseed;
        this->md = md;
    }
    void updateSummaryOutput(Output& output) {
        if (output.ddd) {
            bpexplored += output.bpexplored;
            runtime += output.runtime;
            totaladdruntime += output.addruntime;
            totalspruntime += output.spruntime;
            arctotal += output.arctotal;
            subpathtotal += output.subpathtotal;
            gavgpercenttime *= output.runtime;
            iter += output.iter;
        }
        else {
            bptotal = output.bpexplored;
            enumtime += output.runtime;
            gavgpercenttime /= output.runtime;
        }
        return;
    }
    void calcSummaryOutput() {
        avgarc = (double)arctotal / numseed;
        avgbpexplored = (double)bpexplored / numseed;
        avgenumtime = enumtime / numseed;
        avgruntime = runtime / numseed;
        avgaddruntime = totaladdruntime / numseed;
        avgspruntime = totalspruntime / numseed;
        avgsubpath = (double)subpathtotal / numseed;
        percenttime = 100 * avgruntime / avgenumtime;
        percentbp = 100 * (double)avgbpexplored / bptotal;
        return;
    }
    void writeSOutputCSV() {
        ofstream myFile(filename, ofstream::app);
        myFile << n << ',';
        myFile << T << ',';
        myFile << gtype << ',';
        myFile << ttype << ',';
        myFile << md << ',';
        myFile << avgbpexplored << ',';
        myFile << bptotal << ',';
        myFile << percentbp << ',';
        myFile << pow(gavgpercenttime,1.0/numseed) << ',';
        myFile << avgaddruntime << ',';
        myFile << avgspruntime << ',';
        myFile << avgarc << ',';
        myFile << avgsubpath << ',';
        myFile << avgruntime << ',';
        myFile << avgenumtime << ',';
        myFile << percenttime << ',';
        myFile << iter << endl;
        myFile.close();
        return;
    }
};
void writeOutputCSVHeader(const string& filename) {
    ofstream myFile(filename, ofstream::app);
    myFile << "Is MD?" << ',';
    myFile << "Is DDD?" << ',';
    myFile << "BP Explored" << ',';
    myFile << "Run-Time" << ',';
    myFile << "Add Run-Time" << ',';
    myFile << "SP Run-Time" << ',';
    myFile << "#Arcs in Soln" << ',';
    myFile << "#Subpaths in Soln" << ',';
    myFile << "Iterations" << ',';
    myFile << "Opt Obj" << endl;
    myFile.close();
    return;
}
void writeSOutputCSVHeader(const string& filename) {
    ofstream myFile(filename, ofstream::app);
    myFile << "n" << ',';
    myFile << "T" << ',';
    myFile << "gtype" << ',';
    myFile << "ttype" << ',';
    myFile << "md" << ',';
    myFile << "Avg BP Explored" << ',';
    myFile << "Total BP" << ',';
    myFile << "BP Explored (%)" << ',';
    myFile << "Geom Run-Time (%)" << ',';
    myFile << "Add Run-Time (ms)" << ',';
    myFile << "SP Run-Time (ms)" << ',';
    myFile << "Avg #Arcs in Soln" << ',';
    myFile << "Avg #Subpaths in Soln" << ',';
    myFile << "Avg DDD Run-Time (ms)" << ',';
    myFile << "Avg Enum Run-Time (ms)" << ',';
    myFile << "Run-Time (%)" << ',';
    myFile << "Iterations" << endl;
    myFile.close();
    return;
}

class Graph {
public:
    //Node
    struct Node {
        int nodeID = 0;
        string nodeName = "";
        Node(int ID, string name = "") {
            nodeID = ID;
            if (name == "") name = to_string(ID);
            nodeName = name;
        }
    };

    //Initialization
    Graph(const int n, const int eT, const int sT = 0, const int startn = 0, const int endn = -1) {
        startT = sT;
        endT = eT;
        this->n = n;
        startN = startn;
        if (endn != -1) endN = endn;
        else endN = n - 1;
        addNodes(n);
    }

    void addNode(const int i) {
        Node newNode(i);
        nodes.push_back(newNode);
        //indMap[i] = &newNode;
        return;
    }

    void addNodes(const int n) {
        for (int i = 0; i < n; i++) {
            addNode(i);
        }
        return;
    }

    void addArc(const int startind, const int endind, vector<double> traveltimes) {
        /*
        adds arc information to maps
        to get ifloortraveltimes:
            for each integer going forwards, see where it lands, if it overshoots, assign the previous to ifloor
        to get iceiltraveltimes:
            for each integer going backwards, see where it lands, if it undershoots, assign the next to iceil
        to get minttMap:
            

        */
        //Simple Maps
        type_arc arc = { startind,endind };
        inMap[endind].push_back(startind);
        outMap[startind].push_back(endind);
        ttMap[arc] = traveltimes;
        //ifloorMap
        int T = endT - startT + 1;
        int curr = 0;
        vector<int> ifloortraveltimes(T);
        for (int i = 0; i < T; i++) {
            while (curr <= T - 1 && curr + traveltimes[curr] <= i) {
                curr++;
            }
            ifloortraveltimes[i] = curr - 1;
            //cout << "floor:";
            //cout << i << ',' << traveltimes[i] << ',';
            //cout << i << ',' << ifloortraveltimes[i] << endl;
        }
        ittfloorMap[arc] = ifloortraveltimes;
        //iceilMap
        curr = T - 1;
        vector<int> iceiltraveltimes(T);
        for (int i = T - 1; i >= 0; i--) {
            while (curr >=0 && curr + traveltimes[curr] >= i) {
                curr--;
            }
            iceiltraveltimes[i] = curr + 1;
            //cout << "ceil:";
            //cout << i << ',' << traveltimes[i] << ',';
            //cout << i << ',' << iceiltraveltimes[i] << endl;
        }
        ittceilMap[arc] = iceiltraveltimes;
        //minttMap and minindttMap
        vector<vector<double>> mintt(T, vector<double>(T, INT_MAX));
        vector<vector<int>> minindtt(T, vector<int>(T, -1));
        for (int left = 0; left < T; left++) {
            double currmin = traveltimes[left];
            int currind = left;
            for (int right = left; right < T; right++) {
                if (currmin > traveltimes[right]) {
                    currmin = traveltimes[right];
                    currind = right;
                }
                mintt[left][right] = currmin;
                minindtt[left][right] = currind;
            }
        }
        minttMap[arc] = mintt;
        minindttMap[arc] = minindtt;
        return;
    }

    //Properties
    int startT = 0;
    int endT = 0;
    int n = 0;
    int startN = 0;
    int endN = 0;
    vector<Node> nodes;
    //map<int, Node*> indMap;
    map<int, vector<int>> inMap;
    map<int, vector<int>> outMap;
    map<type_arc, vector<double>> ttMap;
    map<type_arc, vector<int>> ittfloorMap; //returns floor guess for each integer endt
    map<type_arc, vector<int>> ittceilMap; //returns ceil guess for each integer endt
    map<type_arc, vector<vector<double>>> minttMap; //returns table of min travel time for arc in interval
    map<type_arc, vector<vector<int>>> minindttMap; //returns table of departure time of min travel times for arc in interval
    typedef pair<int, double> type_timednode;
    struct timednode_less {
        bool operator() (const type_timednode& lhs, const type_timednode& rhs) const {
            return lhs.second < rhs.second;
        }
    };
    struct timednode_greater {
        bool operator() (const type_timednode& lhs, const type_timednode& rhs) const {
            return lhs.second > rhs.second;
        }
    };

    //Functions
    double TT(const int startind, const int endind, const double startt) {
        /*
        returns travel time of arc from startind to endind with departure time startt
        procedure:
        if startt<startT use travel time at startT
        if startt>endT use travel time at endT
        otherwise interpolate travel time between floor(startt) and ceil(startt)
        */
        pair<int, int> arc = { startind, endind };
        if (startt < startT) return ttMap[arc][0];
        else if (startt > endT) return ttMap[arc][endT - startT];
        else if (startt == (int)startt) return ttMap[arc][(int)startt];
        else {
            double fractpart, intpart;
            fractpart = modf(startt - startT, &intpart);
            return ttMap[arc][(int)intpart] + fractpart * (ttMap[arc][ceil(startt) - startT] - ttMap[arc][(int)intpart]);
        }
    }

    double iTT(const int startind, const int endind, const double endt) {
        /*
        returns travel time of arc from startind to endind with arrival time endt
        procedure:
        if endt<startT+ttMap[arc][0] use travel time at startT
        if endt>endT+ttMap[arc][endT-startT] use travel time at endT
        otherwise interpolate travel time between floor*(endt) and ceil*(endt)
        (separate case for endt integer)
        where * indicates rounding to nearest endpoint of breakpoints, see below.

        get guess of left and right breakpoint: ittMap[arc][floor(endt)],ittMap[arc][floor(endt)]
        refine guesses until we get correct breakpoints leftprev and rightprev
        get leftnext and rightnext

        find a such that ln+a(rn-ln)=endt
        a=(endt-ln)/(rn-ln)
        return ttMap[arc][lp+a(rp-lp)]=ttMap[arc][lp+(endt-ln)/(rn-ln)] since rp-lp=1
        

        */
        pair<int, int> arc = { startind, endind };
        if (endt < startT + ttMap[arc][0]) return ttMap[arc][0];
        if (endt > endT + ttMap[arc][endT - startT]) return ttMap[arc][endT - startT];
        int leftprev, rightprev;
        if (endt > endT) {
            leftprev = ittfloorMap[arc][endT - startT];
            rightprev = ittceilMap[arc][endT - startT];
        }
        else {
            leftprev = ittfloorMap[arc][floor(endt)];
            rightprev = ittceilMap[arc][ceil(endt)];
        }
        if (leftprev == rightprev) return ttMap[arc][leftprev];
        //cout << "before" << leftprev << ',' << rightprev << endl;
        while (rightprev - leftprev > 1) {
            if (leftprev + 1 + ttMap[arc][leftprev + 1] <= endt) leftprev++;
            if (rightprev - 1 + ttMap[arc][rightprev - 1] > endt) rightprev--;
        }
        //cout << "after" << leftprev << ',' << rightprev << endl;
        double leftnext = leftprev + ttMap[arc][leftprev];
        double rightnext = rightprev + ttMap[arc][rightprev];
        //cout << "next" << leftnext << ',' << rightnext << endl;
        return TT(startind, endind, leftprev + (endt - leftnext) / (rightnext - leftnext));
    }

    double FSP(const vector<double>& startts, const vector<double>& endts) {
        vector<double> fspt(n, endT + 1);
        priority_queue<type_timednode, vector<type_timednode>, timednode_greater> pq;
        pq.push({ startN,startts[startN] });
        fspt[startN] = startts[startN];
        while (!pq.empty()) {
            int i = pq.top().first;
            if (i == endN) break;
            //Care double comparison
            if (pq.top().second != fspt[i]) {
                //cout << "(not) popping:" << i << ',' << pq.top().first << endl;
                pq.pop();
                continue;
            }
            //cout << "popping:" << i << ',' << fspt[i] << endl;
            pq.pop();
            for (int j : outMap[i]) {
                //cout << i << j << startts[i] << endts[i] << endl;
                double weight = minTT(i, j, startts[i], endts[i]);
                if (fspt[j] > fspt[i] + weight) {
                    fspt[j] = fspt[i] + weight;
                    pq.push({ j,fspt[j] });
                    //cout << "pushing:" << j << ',' << fspt[j] << endl;
                }
            }
            //printf("Times of Nodes in FSPT\n");
            //for (int i = 0; i < n; i++)
            //    printf("%d \t\t %f\n", i, fspt[i]);
        }
        return fspt[endN] - fspt[startN];

    }

    vector<double> FSPT(const int startind, const double startt, vector<type_arc>& arcs) {
        vector<double> fspt(n, endT + 1);
        vector<int> pred(n, -1);
        priority_queue<type_timednode, vector<type_timednode>, timednode_greater> pq;
        pq.push({ startind,startt });
        fspt[startind] = startt;
        while (!pq.empty()) {
            int i = pq.top().first;
            //Care double comparison
            if (pq.top().second != fspt[i]) {
                //cout << "(not) popping:" << i << ',' << pq.top().first << endl;
                pq.pop();
                continue;
            }
            //cout << "popping:" << i << ',' << fspt[i] << endl;
            pq.pop();
            for (int j : outMap[i]) {
                double weight = TT(i, j, fspt[i]);
                if (fspt[j] > fspt[i] + weight) {
                    fspt[j] = fspt[i] + weight;
                    pred[j] = i;
                    pq.push({ j,fspt[j] });
                    //cout << "pushing:" << j << ',' << fspt[j] << endl;
                }
            }
            //printf("Times of Nodes in FSPT\n");
            //for (int i = 0; i < n; i++)
            //    printf("%d \t\t %f\n", i, fspt[i]);
        }
        for (int j = 0; j < n; j++) {
            if (pred[j] == -1) continue;
            else {
                arcs.push_back({ pred[j],j });
            }
        }
        //printf("Times of Nodes in FSPT\n");
        //for (int i = 0; i < n; i++)
        //    printf("%d \t\t %f\n", i, fspt[i]);
        //printf("Arcs in FSPT\n");
        //for (int ind = 0; ind < arcs.size(); ind++)
        //    printf("%d \t\t %d\n", arcs[ind].first, arcs[ind].second);
        return fspt;
    }
    vector<double> BSPT(const int endind, const double endt, vector<type_arc>& arcs) {
        vector<double> bspt(n, startT - 1);
        vector<int> succ(n, -1);
        priority_queue<type_timednode, vector<type_timednode>, timednode_less> pq;
        pq.push({ endind,endt });
        bspt[endind] = endt;
        while (!pq.empty()) {
            int j = pq.top().first;
            //Care double comparison
            if (pq.top().second != bspt[j]) {
                //cout << "(not) popping:" << j << ',' << pq.top().first << endl;
                pq.pop();
                continue;
            }
            //cout << "popping:" << j << ',' << bspt[j] << endl;
            pq.pop();
            for (int i : inMap[j]) {
                double weight = iTT(i, j, bspt[j]);
                if (bspt[i] < bspt[j] - weight) {
                    bspt[i] = bspt[j] - weight;
                    succ[i] = j;
                    pq.push({ i,bspt[i] });
                    //cout << "pushing:" << i << ',' << bspt[i] << endl;
                }
            }
            //printf("Times of Nodes in BSPT\n");
            //for (int j = 0; j < n; j++)
            //    printf("%d \t\t %f\n", j, bspt[j]);
        }
        for (int i = 0; i < n; i++) {
            if (succ[i] == -1) continue;
            else {
                arcs.push_back({ i,succ[i] });
            }
        }
        //printf("Times of Nodes in BSPT\n");
        //for (int j = 0; j < n; j++)
        //    printf("%d \t\t %f\n", j, bspt[j]);
        //printf("Arcs in BSPT\n");
        //for (int ind = 0; ind < arcs.size(); ind++)
        //    printf("%d \t\t %d\n", arcs[ind].first, arcs[ind].second);
        return bspt;
    }

    double minTT(const int startind, const int endind, double leftt, double rightt) {
        /* Finds departure time of minimum travel time for arc (startind,endind) in interval [leftt,rightt]*/
        leftt = min((double)endT - startT, max(0.0, leftt - startT));
        rightt = min((double)endT - startT, max(0.0, rightt - startT));
        double ans = min(TT(startind, endind, leftt), TT(startind, endind, rightt));
        return min(ans, minttMap[{startind, endind}][ceil(leftt)][floor(rightt)]);
    }

    int minindTT(const int startind, const int endind, const int leftt, const int rightt) {
        return minindttMap[{startind, endind}][leftt][rightt];
    }

    bool checkFIFO() {
        for (int i = 0; i < n; i++) {
            for (int j : outMap[i]) {
                for (int t = 1; t <= endT-startT; t++) {
                    if (ttMap[{i, j}][t] + 1 < ttMap[{i, j}][t - 1]) {
                        cout << i << ',' << j << endl;
                        cout << t - 1 << ',' << ttMap[{i, j}][t - 1] << endl;
                        cout << t << ',' << ttMap[{i, j}][t] << endl;
                        return 0;
                    }
                }
            }
        }
        return 1;
    }
};

class TEN {
public:
    //Common
    //Minimum Duration Structures
    struct MDTimedNode {
        double time = 0;
        MDTimedNode(double t) {
            time = t;
        }
    };
    struct MDTimednode_compare {
        bool operator() (const MDTimedNode& lhs, const MDTimedNode& rhs) const {
            return lhs.time < rhs.time;
        }
    };
    vector<set<MDTimedNode, MDTimednode_compare>> timednodes;
    struct Abspt {
        vector<const MDTimedNode*> nodes;
        vector<double> times;
        double lb = 0;
        double ub = INT_MAX;
        int bpnode = -1;
        int bptime = -1;
        Abspt() :nodes(0), times(0) {
            bpnode = 0;
            bptime = 0;
        }
        Abspt(int i, int t, int n) :nodes(n), times(n) {
            bpnode = i;
            bptime = t;
        }
    };
    struct Abspt_compare {
        bool operator() (const Abspt& lhs, const Abspt& rhs) const {
            return lhs.times[0] < rhs.times[0];
        }
    };
    struct Abspt_comparelb {
        bool operator() (const set<Abspt, Abspt_compare>::iterator& lhs, const set<Abspt, Abspt_compare>::iterator& rhs) const {
            return lhs->lb < rhs->lb;
        }
    };
    struct Abspt_compareub {
        bool operator() (const set<Abspt, Abspt_compare>::iterator& lhs, const set<Abspt, Abspt_compare>::iterator& rhs) const {
            return lhs->ub < rhs->ub;
        }
    };
    set<Abspt,Abspt_compare> abspts;
    typedef set<Abspt, Abspt_compare>::iterator Absptit;
    set<Absptit, Abspt_comparelb> absptlbs;
    set<Absptit, Abspt_compareub> absptubs;

    //Minimum Duration Functions
    void addABSPT(const int bpnode, const int bptime) {
        //cout << "Adding ABSPT:(" << bpnode << ',' << bptime << ')' << endl;
        Abspt newabspt(bpnode, bptime, G.n);
        vector<type_arc> arcs = {};
        double endTime = G.FSPT(bpnode, bptime, arcs)[G.endN];
        arcs = {};
        newabspt.times = G.BSPT(G.endN, endTime, arcs);
        newabspt.times[bpnode] = bptime; //for floating point errors.
        newabspt.ub = newabspt.times[G.endN] - newabspt.times[G.startN];
        for (int i = 0; i < G.n; i++) {
            MDTimedNode timednode(newabspt.times[i]);
            auto it = timednodes[i].insert(timednode);
            newabspt.nodes[i] = &(*it.first);
        }
        //Find LB for new ABSPT
        set<Abspt, Abspt_compare>::iterator it = abspts.upper_bound(newabspt);
        if (it == abspts.end()) {
            newabspt.lb = newabspt.ub;
        }
        else {
            newabspt.lb = G.FSP(newabspt.times, it->times);
        }
        //Update LB for prev ABSPT
        if (it == abspts.begin()) {
            cout << "no previous ABSPT, this should only occur during first or second ABSPT" << endl;
        }
        else {
            it--;
            Abspt prevabspt = *it;
            absptlbs.erase(it);
            absptubs.erase(it);
            abspts.erase(*it);
            prevabspt.lb = G.FSP(prevabspt.times, newabspt.times);
            auto itbool = abspts.insert(prevabspt);
            it = itbool.first;
            //cout << "Updating ABSPT:(" << it->bpnode << ',' << it->bptime << ')' << endl;
            absptlbs.insert(it);
            absptubs.insert(it);
        }
        auto itbool = abspts.insert(newabspt);
        it = itbool.first;
        //cout << "Inserting ABSPT:(" << it->bpnode << ',' << it->bptime << ')' << endl;
        absptlbs.insert(it);
        absptubs.insert(it);
        return;
    }
    void resolveABSPT(Abspt& curr) {
        Abspt currabspt = curr;
        Absptit it = abspts.find(curr);
        absptlbs.erase(it);
        absptubs.erase(it);
        abspts.erase(*it);
        currabspt.lb = currabspt.ub;
        auto itbool = abspts.insert(currabspt);
        it = itbool.first;
        //cout << "Resolving ABSPT:(" << it->bpnode << ',' << it->bptime << ')' << endl;
        absptlbs.insert(it);
        absptubs.insert(it);
    }
    type_bp findBP(Abspt& curr, int option = 1) {
        auto it = abspts.find(curr);
        it++;
        int bpnode = -1;
        int bptime = -1;
        if (it == abspts.end()) return { bpnode,bptime };
        for (int i = 0; i < G.n; i++) {
            if (G.outMap[i].empty()) continue;
            double dleftt = curr.times[i];
            double drightt = it->times[i];
            int leftt, rightt;
            if (dleftt == (int)dleftt) leftt = dleftt + 1;
            else leftt = ceil(dleftt);
            if (drightt == (int)drightt) rightt = drightt - 1;
            else rightt = floor(drightt);
            if (leftt > rightt) continue;
            else {
                bpnode = i;
                if (option == 0) bptime = G.minindTT(i, G.outMap[i][0], leftt, rightt); //min
                else if (option == 1) bptime = (leftt + rightt) / 2; //med
                else if (option == 2) bptime = leftt + rand() % (rightt - leftt + 1); //rand
                return { bpnode, bptime };
            }
        }
        return { bpnode,bptime };
    }
    Output findMD() {
        int bpexplored = 0;
        int iter = 0;
        auto start = chrono::high_resolution_clock::now();
        addABSPT(G.endN, G.endT);
        addABSPT(G.startN, G.startT);
        iter += 2;
        bpexplored += 2;
        auto lbit = absptlbs.begin();
        auto ubit = absptubs.begin();
        //cout << "current lower bound=" << (*lbit)->lb << endl;
        //cout << "current upper bound=" << (*ubit)->ub << endl;
        //printCurrentABSPTs();
        while ((*lbit)->lb != (*ubit)->ub) {
            auto myabspt = **lbit;
            type_bp nextbp = findBP(myabspt);
            if (nextbp.first == -1) {
                resolveABSPT(myabspt);
            }
            else {
                addABSPT(nextbp.first, nextbp.second);
                bpexplored++;
            }
            lbit = absptlbs.begin();
            ubit = absptubs.begin();
            iter++;
            //cout << "current lower bound=" << (*lbit)->lb << endl;
            //cout << "current upper bound=" << (*ubit)->ub << endl;
            //printCurrentABSPTs();
        }
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        cout << "Time taken by function: " << duration.count() << " milliseconds" << endl;
        Output output;
        output.md = 1;
        output.ddd = 1;
        output.bpexplored = bpexplored;
        output.runtime = duration.count();
        output.iter = iter;
        printOptPath(*lbit, output);
        //printOptMD(*lbit);
        output.subpathtotal = 1;
        output.optval = (*lbit)->lb;
        return output;
    }

    //Minimum Travel Time Structures
    struct Mangrove;
    struct TimedNode {
        double time = 0;
        int nodeID = -1;
        int bpnode = -1;
        int bptime = -1;
        bool forward = 0;
        TimedNode(int n, double t, int i, int k, bool forward) {
            nodeID = n;
            time = t;
            bpnode = i;
            bptime = k;
            this->forward = forward;
        }
    };
    struct Timednode_compare {
        bool operator() (const TimedNode& lhs, const TimedNode& rhs) const {
            return lhs.time < rhs.time;
        }
    };
    const TimedNode* origin = nullptr;
    const TimedNode* destination = nullptr;
    typedef set<TimedNode, Timednode_compare> set_timednode;
    vector<vector<set_timednode>> ftimednodes; //takes ownership of forward nodes, used to add waiting arcs, indexed by [node][bpnode];
    vector<vector<set_timednode>> btimednodes; //takes ownership of backward nodes, used to add waiting arcs, indexed by [node][bpnode];
    struct Mangrove {
        mutable vector<const TimedNode*> fnodes;
        mutable vector<const TimedNode*> bnodes;
        mutable vector<double> ftimes;
        mutable vector<double> btimes;
        mutable int bpnode = -1;
        int bptime = -1;
        mutable bool resolved = 0;
        Mangrove(int i, int t, int n) :fnodes(n), bnodes(n), ftimes(n), btimes(n) {
            bpnode = i;
            bptime = t;
        }
    };
    struct Mangrove_compare {
        bool operator() (const Mangrove& lhs, const Mangrove& rhs) const {
            return lhs.bptime < rhs.bptime;
        }
    };
    vector<set<Mangrove, Mangrove_compare>> mangroves; //takes ownership of mangroves, indexed by [bpnode];
    vector<map<int, const Mangrove*>> mangroveMap; //finds mangrove, indexed by [bpnode][bptime];
    typedef pair<const TimedNode*, const TimedNode*> TimedArcptr;
    map<TimedArcptr, double> ttMapLB;
    map<TimedArcptr, double> ttMapUB;
    map<const TimedNode*, vector<const TimedNode*>> outMapLB; //internal arcs for LB
    map<const TimedNode*, vector<const TimedNode*>> outMapUB; //internal arcs for UB
    typedef pair<const TimedNode*, double> TnTT;
    struct TnTT_compare {
        bool operator() (const TnTT& lhs, const TnTT& rhs) const {
            return lhs.second > rhs.second;
        }
    };
    typedef vector<const TimedNode*> type_path;

    //Minimum Travel Time Functions
    const Mangrove* addMangrove(const int bpnode, const int bptime) {
        //cout << "Adding Mangrove:(" << bpnode << ',' << bptime << ')' << endl;
        Mangrove newmangrove(bpnode, bptime, G.n);
        auto nextit = mangroves[bpnode].upper_bound(newmangrove);
        const Mangrove* nextmanptr = (nextit == mangroves[bpnode].end() ? nullptr: &(*nextit));
        const Mangrove* prevmanptr = (nextit == mangroves[bpnode].begin() ? nullptr : &(*prev(nextit, 1)));
        ////Forward
        vector<type_arc> farcs = {};
        newmangrove.ftimes = G.FSPT(bpnode, bptime, farcs);
        //Add nodes to mangrove and ftimednodes
        for (int i = 0; i < G.n; i++) {
            if (newmangrove.ftimes[i] <= G.endT) {
                TimedNode ftimednode(i, newmangrove.ftimes[i], bpnode, bptime, 1);
                auto itbool = ftimednodes[i][bpnode].insert(ftimednode);
                auto it = itbool.first;
                auto currtimednode = &(*it);
                newmangrove.fnodes[i] = currtimednode;
            }
            else {
                newmangrove.fnodes[i] = nullptr;
            }
        }
        ////Backward
        vector<type_arc> barcs = {};
        newmangrove.btimes = G.BSPT(bpnode, bptime, barcs);
        //Add nodes to mangrove and btimednodes
        for (int i = 0; i < G.n; i++) {
            if (newmangrove.btimes[i] >= G.startT) {
                TimedNode btimednode(i, newmangrove.btimes[i], bpnode, bptime, 0);
                auto itbool = btimednodes[i][bpnode].insert(btimednode);
                auto it = itbool.first;
                auto currtimednode = &(*it);
                newmangrove.bnodes[i] = currtimednode;
            }
            else {
                newmangrove.bnodes[i] = nullptr;
            }
        }

        //Split into two cases depending on whether added mangrove is resolved(no next mangrove or next mangrove is one unit away)
        //Add entries for outMapLB, outMapUB, ttMapLB, ttMapUB
        if (nextmanptr==nullptr || nextmanptr->bptime == bptime + 1) {
            //Do everything resolved
            newmangrove.resolved = 1;
            //Add arcs in FSPT to outMapLB and outMapUB with travel times in ttMapLB and ttMapUB
            for (type_arc arc : farcs) {
                int i = arc.first, j = arc.second;
                if (newmangrove.ftimes[i] <= G.endT && newmangrove.ftimes[j] <= G.endT) {
                    const TimedNode* timednodei = newmangrove.fnodes[i];
                    const TimedNode* timednodej = newmangrove.fnodes[j];
                    outMapLB[timednodei].push_back(timednodej);
                    outMapUB[timednodei].push_back(timednodej);
                    ttMapUB[{timednodei, timednodej}] = G.TT(i, j, newmangrove.ftimes[i]);
                    ttMapLB[{timednodei, timednodej}] = ttMapUB[{timednodei, timednodej}];
                }
            }
            //Add arcs in BSPT to outMapLB and outMapUB with travel times in ttMapLB and ttMapUB
            for (type_arc arc : barcs) {
                int i = arc.first, j = arc.second;
                if (newmangrove.btimes[i] >= G.startT && newmangrove.btimes[j] >= G.startT) {
                    const TimedNode* timednodei = newmangrove.bnodes[i];
                    const TimedNode* timednodej = newmangrove.bnodes[j];
                    outMapLB[timednodei].push_back(timednodej);
                    outMapUB[timednodei].push_back(timednodej);
                    ttMapUB[{timednodei, timednodej}] = G.TT(i, j, newmangrove.btimes[i]);
                    ttMapLB[{timednodei, timednodej}] = ttMapUB[{timednodei, timednodej}];
                }
            }
        }
        else {
            //Do everything normally
            //Add arcs in FSPT to outMapLB and outMapUB
            for (int i = 0; i < G.n; i++) {
                const TimedNode* ftimednode = newmangrove.fnodes[i];
                if (ftimednode == nullptr) continue;
                for (int j : G.outMap[i]) {
                    if (newmangrove.ftimes[j] <= G.endT) {
                        const TimedNode* timednodei = newmangrove.fnodes[i];
                        const TimedNode* timednodej = newmangrove.fnodes[j];
                        outMapLB[timednodei].push_back(timednodej);
                        ttMapLB[{timednodei, timednodej}] = G.minTT(i, j, newmangrove.ftimes[i], nextmanptr->ftimes[i]);
                    }
                }
            }
            for (type_arc arc : farcs) {
                int i = arc.first, j = arc.second;
                if (newmangrove.ftimes[i] <= G.endT && newmangrove.ftimes[j] <= G.endT) {
                    const TimedNode* timednodei = newmangrove.fnodes[i];
                    const TimedNode* timednodej = newmangrove.fnodes[j];
                    outMapUB[timednodei].push_back(timednodej);
                    ttMapUB[{timednodei, timednodej}] = G.TT(i, j, newmangrove.ftimes[i]);
                }
            }
            //Add arcs in BSPT to outMapLB and outMapUB
            for (int i = 0; i < G.n; i++) {
                const TimedNode* btimednode = newmangrove.bnodes[i];
                if (btimednode == nullptr) continue;
                for (int j : G.outMap[i]) {
                    if (newmangrove.btimes[j] >= G.startT) {
                        const TimedNode* timednodei = newmangrove.bnodes[i];
                        const TimedNode* timednodej = newmangrove.bnodes[j];
                        outMapLB[timednodei].push_back(timednodej);
                        ttMapLB[{timednodei, timednodej}] = G.minTT(i, j, newmangrove.btimes[i], nextmanptr->btimes[i]);
                    }
                }
            }
            for (type_arc arc : barcs) {
                int i = arc.first, j = arc.second;
                if (newmangrove.btimes[i] >= G.startT && newmangrove.btimes[j] >= G.startT) {
                    const TimedNode* timednodei = newmangrove.bnodes[i];
                    const TimedNode* timednodej = newmangrove.bnodes[j];
                    outMapUB[timednodei].push_back(timednodej);
                    ttMapUB[{timednodei, timednodej}] = G.TT(i, j, newmangrove.btimes[i]);
                }
            }
        }
        //Update costs for prev mangrove ttMapLB and ttMapUB depending on whether prev mangrove needs to be resolved or not
        if (prevmanptr!=nullptr) {
            if (prevmanptr->bptime == bptime - 1) {
                //Resolve
                prevmanptr->resolved = 1;
                for (const TimedNode* ftimednode : prevmanptr->fnodes) {
                    if (ftimednode == nullptr) continue;
                    outMapLB[ftimednode] = outMapUB[ftimednode];
                    for (const TimedNode* nextftimednode : outMapUB[ftimednode]) {
                        TimedArcptr timedarcptr = { ftimednode,nextftimednode };
                        ttMapLB[timedarcptr] = ttMapUB[timedarcptr];
                    }
                }
                for (const TimedNode* btimednode : prevmanptr->bnodes) {
                    if (btimednode == nullptr) continue;
                    outMapLB[btimednode] = outMapUB[btimednode];
                    for (const TimedNode* nextbtimednode : outMapUB[btimednode]) {
                        TimedArcptr timedarcptr = { btimednode,nextbtimednode };
                        ttMapLB[timedarcptr] = ttMapUB[timedarcptr];
                    }
                }
            }
            else {
                //Update ttMapLB
                for (const TimedNode* ftimednode : prevmanptr->fnodes) {
                    if (ftimednode == nullptr) continue;
                    int i = ftimednode->nodeID;
                    double leftt = prevmanptr->ftimes[i];
                    for (const TimedNode* nextftimednode : outMapLB[ftimednode]) {
                        int j = nextftimednode->nodeID;
                        double rightt = newmangrove.ftimes[j];
                        TimedArcptr timedarcptr = { ftimednode,nextftimednode };
                        ttMapLB[timedarcptr] = G.minTT(i, j, leftt, rightt);
                    }
                }
                for (const TimedNode* btimednode : prevmanptr->bnodes) {
                    if (btimednode == nullptr) continue;
                    int i = btimednode->nodeID;
                    double leftt = prevmanptr->btimes[i];
                    for (const TimedNode* nextbtimednode : outMapLB[btimednode]) {
                        int j = nextbtimednode->nodeID;
                        double rightt = newmangrove.btimes[j];
                        TimedArcptr timedarcptr = { btimednode,nextbtimednode };
                        ttMapLB[timedarcptr] = G.minTT(i, j, leftt, rightt);
                    }
                }
            }
        }
        auto itbool = mangroves[bpnode].insert(newmangrove);
        const Mangrove* manptr = &(*itbool.first);
        mangroveMap[bpnode][bptime] = manptr;
        return manptr;
    }
    pair<type_path,double> findLB() {
        map<const TimedNode*,double> dp;
        map<const TimedNode*, const TimedNode*> pred;
        priority_queue<TnTT, vector<TnTT>, TnTT_compare> pq;
        pq.push({ origin,0 });
        dp[origin] = 0;
        while (!pq.empty()) {
            const TimedNode* ptr = pq.top().first;
            //cout << "Popping: (" << ptr->nodeID << ',' << ptr->nodeID << ") from: (" << ptr->bpnode << ',' << ptr->bptime << ")" << endl;
            //cout << "DP value: " << pq.top().second << endl;
            if (ptr == destination) break;
            //Care double comparison
            if (pq.top().second > dp[ptr]) {
                pq.pop();
                continue;
            }
            pq.pop();
            for (const TimedNode* nextptr : outMapLB[ptr]) {
                double weight = ttMapLB[{ptr, nextptr}];
                auto it = dp.find(nextptr);
                if (it == dp.end()) dp[nextptr] = INT_MAX;
                if (dp[nextptr] > dp[ptr] + weight) {
                    dp[nextptr] = dp[ptr] + weight;
                    pq.push({ nextptr,dp[nextptr] });
                    //cout << "Pushing: (" << nextptr->nodeID << ',' << nextptr->nodeID << ") from: (" << nextptr->bpnode << ',' << nextptr->bptime << ")" << endl;
                    //cout << "DP value: " << dp[nextptr] << endl;
                    pred[nextptr] = ptr;
                }
            }
            int node = ptr->nodeID;
            int i = ptr->bpnode;
            if (ptr->forward) {
                for (int j = 0; j < G.n; j++) {
                    if (btimednodes[node][j].empty()) continue;
                    auto nextit = btimednodes[node][j].upper_bound(*ptr);
                    const TimedNode* nextptr;
                    if (j == i) {
                        //Waiting arc to next copy
                        if (nextit == btimednodes[node][j].end()) continue;
                        else nextptr = &(*nextit);
                    }
                    else {
                        //Waiting arc from i-FSPT to j-BSPT
                        if (nextit == btimednodes[node][j].begin()) nextptr = &(*nextit);
                        else if (nextit == btimednodes[node][j].end()) continue;
                        else {
                            auto previt = prev(nextit, 1);
                            if (mangroveMap[j][previt->bptime]->resolved) nextptr = &(*nextit);
                            else nextptr = &(*previt);
                        }
                    }
                    if (nextptr == nullptr) continue;
                    auto it = dp.find(nextptr);
                    if (it == dp.end()) dp[nextptr] = INT_MAX;
                    if (dp[nextptr] > dp[ptr]) {
                        dp[nextptr] = dp[ptr];
                        pq.push({ nextptr,dp[nextptr] });
                        //cout << "Pushing: (" << nextptr->nodeID << ',' << nextptr->nodeID << ") from: (" << nextptr->bpnode << ',' << nextptr->bptime << ")" << endl;
                        //cout << "DP value: " << dp[nextptr] << endl;
                        pred[nextptr] = ptr;
                    }
                }
            }
            else {
                if (i == node) {
                    //Waiting arc from BSPT to FSPT through node
                    auto nextit = ftimednodes[node][i].lower_bound(*ptr);
                    const TimedNode* nextptr;
                    assert(nextit != ftimednodes[node][i].end());
                    nextptr = &(*nextit);
                    auto it = dp.find(nextptr);
                    if (it == dp.end()) dp[nextptr] = INT_MAX;
                    if (dp[nextptr] > dp[ptr]) {
                        dp[nextptr] = dp[ptr];
                        pq.push({ nextptr,dp[nextptr] });
                        //cout << "Pushing: (" << nextptr->nodeID << ',' << nextptr->nodeID << ") from: (" << nextptr->bpnode << ',' << nextptr->bptime << ")" << endl;
                        //cout << "DP value: " << dp[nextptr] << endl;
                        pred[nextptr] = ptr;
                    }
                }
            }
        }
        type_path path;
        const TimedNode* curr = destination;
        while (curr != origin) {
            path.push_back(curr);
            curr = pred[curr];
        }
        path.push_back(curr);
        reverse(path.begin(), path.end());
        return { path,dp[destination] - dp[origin] };
    }
    pair<type_path, double> findUB() {
        map<const TimedNode*, double> dp;
        map<const TimedNode*, const TimedNode*> pred;
        priority_queue<TnTT, vector<TnTT>, TnTT_compare> pq;
        pq.push({ origin,0 });
        dp[origin] = 0;
        while (!pq.empty()) {
            const TimedNode* ptr = pq.top().first;
            //cout << "Popping: (" << ptr->nodeID << ',' << ptr->nodeID << ") from: (" << ptr->bpnode << ',' << ptr->bptime << ")" << endl;
            //cout << "DP value: " << pq.top().second << endl;
            if (ptr == destination) break;
            //Care double comparison
            if (pq.top().second > dp[ptr]) {
                pq.pop();
                continue;
            }
            pq.pop();
            for (const TimedNode* nextptr : outMapUB[ptr]) {
                double weight = ttMapUB[{ptr, nextptr}];
                auto it = dp.find(nextptr);
                if (it == dp.end()) dp[nextptr] = INT_MAX;
                if (dp[nextptr] > dp[ptr] + weight) {
                    dp[nextptr] = dp[ptr] + weight;
                    pq.push({ nextptr,dp[nextptr] });
                    //cout << "Pushing: (" << nextptr->nodeID << ',' << nextptr->nodeID << ") from: (" << nextptr->bpnode << ',' << nextptr->bptime << ")" << endl;
                    //cout << "DP value: " << dp[nextptr] << endl;
                    pred[nextptr] = ptr;
                }
            }
            int node = ptr->nodeID;
            int i = ptr->bpnode;
            if (ptr->forward) {
                for (int j = 0; j < G.n; j++) {
                    if (btimednodes[node][j].empty()) continue;
                    auto nextit = btimednodes[node][j].upper_bound(*ptr);
                    const TimedNode* nextptr;
                    if (j == i) {
                        //Waiting arc to next copy
                        if (nextit == btimednodes[node][j].end()) continue;
                        else nextptr = &(*nextit);
                    }
                    else {
                        //Waiting arc from i-FSPT to j-BSPT
                        if (nextit == btimednodes[node][j].begin()) nextptr = &(*nextit);
                        else if (nextit == btimednodes[node][j].end()) continue;
                        nextptr = &(*nextit);
                    }
                    if (nextptr == nullptr) continue;
                    auto it = dp.find(nextptr);
                    if (it == dp.end()) dp[nextptr] = INT_MAX;
                    if (dp[nextptr] > dp[ptr]) {
                        dp[nextptr] = dp[ptr];
                        pq.push({ nextptr,dp[nextptr] });
                        //cout << "Pushing: (" << nextptr->nodeID << ',' << nextptr->nodeID << ") from: (" << nextptr->bpnode << ',' << nextptr->bptime << ")" << endl;
                        //cout << "DP value: " << dp[nextptr] << endl;
                        pred[nextptr] = ptr;
                    }
                }
            }
            else {
                if (i == node) {
                    //Waiting arc from BSPT to FSPT through node
                    auto nextit = ftimednodes[node][i].lower_bound(*ptr);
                    const TimedNode* nextptr;
                    assert(nextit != ftimednodes[node][i].end());
                    nextptr = &(*nextit);
                    auto it = dp.find(nextptr);
                    if (it == dp.end()) dp[nextptr] = INT_MAX;
                    if (dp[nextptr] > dp[ptr]) {
                        dp[nextptr] = dp[ptr];
                        pq.push({ nextptr,dp[nextptr] });
                        //cout << "Pushing: (" << nextptr->nodeID << ',' << nextptr->nodeID << ") from: (" << nextptr->bpnode << ',' << nextptr->bptime << ")" << endl;
                        //cout << "DP value: " << dp[nextptr] << endl;
                        pred[nextptr] = ptr;
                    }
                }
            }
        }
        type_path path;
        const TimedNode* curr = destination;
        while (curr != origin) {
            path.push_back(curr);
            curr = pred[curr];
        }
        path.push_back(curr);
        reverse(path.begin(), path.end());
        return { path,dp[destination] - dp[origin] };
    }
    vector<type_bp> findBP(const Mangrove* curr, int addmult = 4, int option = 0) {
        vector<type_bp> multbps;
        int bpnode = curr->bpnode;
        int bptime = -1;
        auto nextit = mangroves[curr->bpnode].upper_bound(*curr);
        if (nextit == mangroves[curr->bpnode].end()) return multbps;
        //Check if any breakpoints in between
        int leftt = curr->bptime;
        int rightt = nextit->bptime;
        leftt++, rightt--;
        if (leftt <= rightt){
            if (addmult == 4) {
                stack<pair<int,int>> st;
                st.push({ leftt,rightt });
                while (!st.empty()) {
                    pair<int, int> interval = st.top();
                    st.pop();
                    leftt = interval.first;
                    rightt = interval.second;
                    if (G.outMap[bpnode].empty()) {
                        if (option == 0) bptime = G.minindTT(G.inMap[bpnode][0], bpnode, leftt, rightt); //min
                        else if (option == 1) bptime = (leftt + rightt) / 2; //med
                        else if (option == 2) bptime = leftt + rand() % (rightt - leftt + 1); //rand
                    }
                    else {
                        if (option == 0) bptime = G.minindTT(bpnode, G.outMap[bpnode][0], leftt, rightt); //min
                        else if (option == 1) bptime = (leftt + rightt) / 2; //med
                        else if (option == 2) bptime = leftt + rand() % (rightt - leftt + 1); //rand
                    }
                    if (bptime - 1 > leftt) {
                        st.push({ leftt ,bptime - 2 });
                    }
                    if (bptime + 1 < rightt) {
                        st.push({ bptime + 2,rightt });
                    }
                    multbps.push_back({ bpnode,bptime });
                    if (bptime - 1 >= leftt) multbps.push_back({ bpnode,bptime - 1 });
                    if (bptime + 1 <= rightt) multbps.push_back({ bpnode,bptime + 1 });
                }
            }
            else if (addmult == 5) {
                for (int i = leftt; i <= rightt; i++) {
                    multbps.push_back({ bpnode,i });
                }
            }
            else {
                if (G.outMap[bpnode].empty()) {
                    if (option == 0) bptime = G.minindTT(G.inMap[bpnode][0], bpnode, leftt, rightt); //min
                    else if (option == 1) bptime = (leftt + rightt) / 2; //med
                    else if (option == 2) bptime = leftt + rand() % (rightt - leftt + 1); //rand
                }
                else {
                    if (option == 0) bptime = G.minindTT(bpnode, G.outMap[bpnode][0], leftt, rightt); //min
                    else if (option == 1) bptime = (leftt + rightt) / 2; //med
                    else if (option == 2) bptime = leftt + rand() % (rightt - leftt + 1); //rand
                }
                multbps.push_back({ bpnode,bptime });
                if (addmult >= 1 && bptime - 1 >= leftt) multbps.push_back({ bpnode,bptime - 1 });
                if (addmult >= 2 && bptime + 1 <= rightt) multbps.push_back({ bpnode,bptime + 1 });
                if (addmult >= 3 && leftt < bptime - 1) multbps.push_back({ bpnode,leftt });
            }
        }
        return multbps;
    }
    set<type_bp> findBP(type_path& path) {
        set<type_bp> bps;
        for (int i = 1; i < path.size(); i++) { 
            const TimedNode* timednode = path[i];
            if (path[i]->bpnode != path[i - 1]->bpnode) {
                //cout << "Not Investigating: (" << timednode->bpnode << ',' << timednode->bptime << ')' << endl;
                continue;
            }
            //cout << "Investigating: (" << timednode->bpnode << ',' << timednode->bptime << ')' << endl;
            vector<type_bp> multbp = findBP(mangroveMap[timednode->bpnode][timednode->bptime]);
            for (type_bp bp: multbp) {
                bps.insert(bp);
            }
        }
        return bps;
    }
    bool isResolved(type_path& path) {
        bool wait = 1;
        for (int i = 1; i < path.size(); i++) {
            if (path[i]->nodeID == path[i - 1]->nodeID) {
                wait = 1;
            }
            else if (wait) {
                wait = 0;
                //cout << path[i]->bpnode << path[i]->bptime << endl;
                if (mangroveMap[path[i]->bpnode][path[i]->bptime]->resolved == false) {
                    return false;
                }
            }
        }
        return true;
    }
    Output findMTT() {
        int bpexplored = 0;
        int iter = 0;
        auto start = chrono::high_resolution_clock::now();
        auto start2 = chrono::high_resolution_clock::now();
        const Mangrove* lastmangrove = addMangrove(G.endN, G.endT);
        auto stop2 = chrono::high_resolution_clock::now();
        auto duration2 = chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
        bpexplored++;
        destination = lastmangrove->bnodes[G.endN];
        for (int i = 0; i < G.n; i++) {
            if (i == G.endN) continue;
            start2 = chrono::high_resolution_clock::now();
            addMangrove(i, floor(lastmangrove->btimes[i]));
            stop2 = chrono::high_resolution_clock::now();
            duration2 += chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
            bpexplored++;
        }
        start2 = chrono::high_resolution_clock::now();
        const Mangrove* firstmangrove = addMangrove(G.startN, G.startT);
        stop2 = chrono::high_resolution_clock::now();
        duration2 += chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
        bpexplored++;
        origin = firstmangrove->fnodes[G.startN];
        for (int i = 0; i < G.n; i++) {
            if (i == G.startN) continue;
            start2 = chrono::high_resolution_clock::now();
            addMangrove(i, ceil(firstmangrove->ftimes[i]));
            stop2 = chrono::high_resolution_clock::now();
            duration2 += chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
            bpexplored++;
        }
        //Add breakpoints associated with minimums or points near minimums
        if (false) {
            for (int i = 0; i < G.n; i++) {
                int j = G.outMap[i][0];
                int leftt = ceil(firstmangrove->ftimes[i]) + 1;
                int rightt = floor(lastmangrove->btimes[i]) - 1;
                double reftime = G.minttMap[{i, j}][leftt][rightt];
                queue<vector<int>> q;
                q.push({ G.minindttMap[{i, j}][leftt][rightt],leftt,rightt });
                while (!q.empty()) {
                    vector<int> v = q.front();
                    int t = v[0];
                    leftt = v[1];
                    rightt = v[2];
                    q.pop();
                    start2 = chrono::high_resolution_clock::now();
                    addMangrove(i, t);
                    if (t > leftt) {
                        addMangrove(i, t - 1);
                        bpexplored++;
                    }
                    if (t - 1 > leftt) {
                        addMangrove(i, t - 2);
                        bpexplored++;
                    }
                    if (t < rightt) {
                        addMangrove(i, t + 1);
                        bpexplored++;
                    }
                    if (t + 1 < rightt) {
                        addMangrove(i, t + 2);
                        bpexplored++;
                    }
                    if (t + 2 < rightt) {
                        addMangrove(i, t + 3);
                        bpexplored++;
                    }
                    stop2 = chrono::high_resolution_clock::now();
                    duration2 += chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
                    bpexplored++;
                    if (leftt < t - 2 && G.minttMap[{i, j}][leftt][t - 3] < 1.03 * reftime) {
                        q.push({ G.minindttMap[{i, j}][leftt][t - 3],leftt,t - 3 });
                    }
                    if (rightt > t + 3 && G.minttMap[{i, j}][t + 4][rightt] < 1.03 * reftime) {
                        q.push({ G.minindttMap[{i, j}][t + 4][rightt],t + 4,rightt });
                    }
                }
            }
        }
        //printCurrentMangroves();
        auto start3 = chrono::high_resolution_clock::now();
        pair<type_path, double> pathLB = findLB();
        auto stop3 = chrono::high_resolution_clock::now();
        auto duration3 = chrono::duration_cast<chrono::milliseconds>(stop3 - start3);
        iter++;
        //pair<type_path, double> pathUB = findUB();
        //printPathInfo(pathLB, pathUB);
        //printCurrentMangroves();
        while (!isResolved(pathLB.first)) {
            set<type_bp> nextbps = findBP(pathLB.first);
            //printBPtoAdd(nextbps);
            for (type_bp nextbp : nextbps) {
                start2 = chrono::high_resolution_clock::now();
                addMangrove(nextbp.first, nextbp.second);
                stop2 = chrono::high_resolution_clock::now();
                duration2 += chrono::duration_cast<chrono::milliseconds>(stop2 - start2);
                bpexplored++;
            }
            start3 = chrono::high_resolution_clock::now();
            pathLB = findLB();
            stop3 = chrono::high_resolution_clock::now();
            //cout << chrono::duration_cast<chrono::milliseconds>(stop3 - start3).count() << endl;
            duration3 += chrono::duration_cast<chrono::milliseconds>(stop3 - start3);
            //pair<type_path, double> pathUB = findUB();
            //printPathInfo(pathLB, pathUB);
            iter++;
        }
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        cout << "Time taken by function: " << duration.count() << " milliseconds" << endl;
        cout << "Time taken by adding mangrove: " << duration2.count() << " milliseconds" << endl;
        cout << "Time taken by calc: " << duration3.count() << " milliseconds" << endl;
        cout << "Optimal Solution is: " << pathLB.second << endl;
        //pair<type_path, double> pathUB = findUB();
        //printPathInfo(pathLB, pathUB);
        Output output;
        output.md = 0;
        output.ddd = 1;
        output.bpexplored = bpexplored;
        output.runtime = duration.count();
        output.spruntime = duration3.count();
        output.addruntime = duration2.count();
        output.iter = iter;
        printLBPath(pathLB.first, output);
        output.optval = pathLB.second;
        return output;
    }

    //Enumeration Functions
    Output findEnumMD() {
        auto start = chrono::high_resolution_clock::now();
        const Mangrove* optmangrove = nullptr;
        double optval;
        const Mangrove* lastmangrove = addEnumMangrove(G.endN, G.endT);
        optmangrove = lastmangrove;
        optval = lastmangrove->ftimes[G.endN] - lastmangrove->btimes[G.startN];
        const Mangrove* firstmangrove = addEnumMangrove(G.startN, G.startT);
        double temp = firstmangrove->ftimes[G.endN] - firstmangrove->btimes[G.startN];
        if (temp < optval) {
            optmangrove = firstmangrove;
            optval = temp;
        }
        for (int i = 0; i < G.n; i++) {
            int leftt, rightt;
            if ((int)lastmangrove->btimes[i] == lastmangrove->btimes[i]) {
                rightt = lastmangrove->btimes[i] - 1;
            }
            else {
                rightt = floor(lastmangrove->btimes[i]);
            }
            if ((int)firstmangrove->ftimes[i] == firstmangrove->ftimes[i]) {
                leftt = firstmangrove->ftimes[i] + 1;
            }
            else {
                leftt = ceil(firstmangrove->ftimes[i]);
            }
            for (int t = leftt; t <= rightt; t++) {
                const Mangrove* newmangrove = addEnumMangrove(i, t);
                double temp = newmangrove->ftimes[G.endN] - newmangrove->btimes[G.startN];
                if (temp < optval) {
                    optmangrove = newmangrove;
                    optval = temp;
                }
            }
        }
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        cout << "Time taken by function: " << duration.count() << " milliseconds" << endl;
        cout << "Optimal Solution is: " << optval << endl;
        Output output;
        output.md = 1;
        output.ddd = 0;
        output.bpexplored = (G.n - 1) * (G.endT - G.startT - 1) + 2;
        output.runtime = duration.count();
        printOptPath(*optmangrove, output);
        output.subpathtotal = 1;
        output.optval = optval;
        return output;
    }
    const Mangrove* addEnumMangrove(const int bpnode, const int bptime) {
        //cout << "Adding Mangrove:(" << bpnode << ',' << bptime << ')' << endl;
        Mangrove newmangrove(bpnode, bptime, G.n);
        ////Forward
        vector<type_arc> farcs = {};
        newmangrove.ftimes = G.FSPT(bpnode, bptime, farcs);
        //Add nodes to mangrove and ftimednodes
        for (int i = 0; i < G.n; i++) {
            if (newmangrove.ftimes[i] <= G.endT) {
                TimedNode ftimednode(i, newmangrove.ftimes[i], bpnode, bptime, 1);
                auto itbool = ftimednodes[i][bpnode].insert(ftimednode);
                auto it = itbool.first;
                auto currtimednode = &(*it);
                newmangrove.fnodes[i] = currtimednode;
            }
            else {
                newmangrove.fnodes[i] = nullptr;
            }
        }
        ////Backward
        vector<type_arc> barcs = {};
        newmangrove.btimes = G.BSPT(bpnode, bptime, barcs);
        //Add nodes to mangrove and btimednodes
        for (int i = 0; i < G.n; i++) {
            if (newmangrove.btimes[i] >= G.startT) {
                TimedNode btimednode(i, newmangrove.btimes[i], bpnode, bptime, 0);
                auto itbool = btimednodes[i][bpnode].insert(btimednode);
                auto it = itbool.first;
                auto currtimednode = &(*it);
                newmangrove.bnodes[i] = currtimednode;
            }
            else {
                newmangrove.bnodes[i] = nullptr;
            }
        }
        //Do everything resolved
        newmangrove.resolved = 1;
        //Add arcs in FSPT to outMapLB and outMapUB with travel times in ttMapLB and ttMapUB
        for (type_arc arc : farcs) {
            int i = arc.first, j = arc.second;
            if (newmangrove.ftimes[i] <= G.endT && newmangrove.ftimes[j] <= G.endT) {
                const TimedNode* timednodei = newmangrove.fnodes[i];
                const TimedNode* timednodej = newmangrove.fnodes[j];
                outMapUB[timednodei].push_back(timednodej);
                ttMapUB[{timednodei, timednodej}] = G.TT(i, j, newmangrove.ftimes[i]);
            }
        }
        //Add arcs in BSPT to outMapLB and outMapUB with travel times in ttMapLB and ttMapUB
        for (type_arc arc : barcs) {
            int i = arc.first, j = arc.second;
            if (newmangrove.btimes[i] >= G.startT && newmangrove.btimes[j] >= G.startT) {
                const TimedNode* timednodei = newmangrove.bnodes[i];
                const TimedNode* timednodej = newmangrove.bnodes[j];
                outMapUB[timednodei].push_back(timednodej);
                ttMapUB[{timednodei, timednodej}] = G.TT(i, j, newmangrove.btimes[i]);
            }
        }
        auto itbool = mangroves[bpnode].insert(newmangrove);
        const Mangrove* manptr = &(*itbool.first);
        mangroveMap[bpnode][bptime] = manptr;
        return manptr;
    }
    Output findEnumMTT() {
        auto start = chrono::high_resolution_clock::now();
        const Mangrove* lastmangrove = addEnumMangrove(G.endN, G.endT);
        destination = lastmangrove->bnodes[G.endN];
        const Mangrove* firstmangrove = addEnumMangrove(G.startN, G.startT);
        origin = firstmangrove->fnodes[G.startN];
        for (int i = 0; i < G.n; i++) {
            int leftt, rightt;
            if ((int)lastmangrove->btimes[i] == lastmangrove->btimes[i]) {
                rightt = lastmangrove->btimes[i] - 1;
            }
            else {
                rightt = floor(lastmangrove->btimes[i]);
            }
            if ((int)firstmangrove->ftimes[i] == firstmangrove->ftimes[i]) {
                leftt = firstmangrove->ftimes[i] + 1;
            }
            else {
                leftt = ceil(firstmangrove->ftimes[i]);
            }
            for (int t = leftt; t <= rightt; t++) {
                const Mangrove* newmangrove = addEnumMangrove(i, t);
            }
        }
        pair<type_path, double> pathOpt = findUB();
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        cout << "Time taken by function: " << duration.count() << " milliseconds" << endl;
        cout << "Optimal Solution is: " << pathOpt.second << endl;
        Output output;
        output.md = 0;
        output.ddd = 0;
        output.bpexplored = (G.n - 1) * (G.endT - G.startT - 1) + 2;
        output.runtime = duration.count();
        printUBPath(pathOpt.first, output);
        output.optval = pathOpt.second;
        return output;
    };

    //Print Functions
    void printCurrentABSPTs() {
        printBar();
        cout << "Printing Current ABSPTs:" << endl;
        for (auto abspt : abspts) {
            cout << "(" << abspt.bpnode << ',' << abspt.bptime << "), resolved = " << (bool)(abspt.lb == abspt.ub);
            cout<< ", lb = " << abspt.lb << ", ub = " << abspt.ub << endl;
        }
        printBar();
    }
    void printOptPath(const Absptit& abspt, Output& output, bool flag=0) {
        if (flag) printBar();
        if (flag) cout << "Printing Opt Path: " << endl;
        output.arctotal = 0;
        vector<type_arc> arcs = {};
        vector<double> times = G.FSPT(G.startN, abspt->times[G.startN], arcs);
        int curr = G.endN;
        stack<pair<int, double>> st;
        while (curr != G.startN) {
            for (type_arc arc : arcs) {
                if (curr == arc.second) {
                    st.push({ curr,times[curr] });
                    curr = arc.first;
                    break;
                }
            }
        }
        if (flag) cout << "(" << curr << ',' << times[curr] << ")" << endl;
        while (!st.empty()) {
            auto currtimes = st.top();
            if(flag) cout << "(" << currtimes.first << ',' << currtimes.second << ")" << endl;
            output.arctotal++;
            st.pop();
        }
        if (flag) printBar();
    }
    void printOptPath(const Mangrove& mangrove, Output& output, bool flag = 0) {
        if (flag) printBar();
        if (flag) cout << "Printing Opt Path: " << endl;
        output.arctotal = 0;
        vector<type_arc> arcs = {};
        vector<double> times = G.FSPT(G.startN, mangrove.btimes[G.startN], arcs);
        int curr = G.endN;
        stack<pair<int, double>> st;
        while (curr != G.startN) {
            for (type_arc arc : arcs) {
                if (curr == arc.second) {
                    st.push({ curr,times[curr] });
                    curr = arc.first;
                    break;
                }
            }
        }
        if (flag) cout << "(" << curr << ',' << times[curr] << ")" << endl;
        while (!st.empty()) {
            auto currtimes = st.top();
            if (flag) cout << "(" << currtimes.first << ',' << currtimes.second << ")" << endl;
            output.arctotal++;
            st.pop();
        }
        if (flag) printBar();
    }
    void printOptMD(const Absptit& abspt) {
        cout << "Optimal MD starting time:" << abspt->times[G.startN] << endl;
        cout << "Optimal MD val:" << abspt->lb << endl;
        cout << "Optimal MD breakpoint: (" << abspt->bpnode << ',' << abspt->bptime << ')' << endl;
    }
    void printPathInfo(pair<type_path, double>& pathLB, pair<type_path, double>& pathUB) {
        Output output;
        printLBPath(pathLB.first, output, 1);
        cout << "current lower bound = " << pathLB.second << endl;
        printUBPath(pathUB.first, output, 1);
        cout << "current upper bound = " << pathUB.second << endl;
    }
    void printLBPath(type_path path, Output& output, bool flag = 0) {
        if (flag) printBar();
        if (flag) cout << "Printing LB Path: " << endl;
        output.arctotal = 0;
        output.subpathtotal = 0;
        bool waiting = 1;
        for (int i = 1; i < path.size(); i++) {
            const TimedNode* node = path[i];
            if (path[i - 1]->nodeID != path[i]->nodeID) {
                TimedArcptr arc = { path[i - 1],path[i] };
                if (flag) cout << "arc cost: " << ttMapLB[arc] << endl;
                double UBcost = (path[i - 1]->nodeID == path[i]->nodeID) ? 0 : path[i]->time - path[i - 1]->time;
                if (flag) cout << "UB arc cost: " << UBcost << endl;
                output.arctotal++;
                if (flag) cout << "at: (" << node->nodeID << ',' << node->time << "), in mangrove: (" << node->bpnode << ',' << node->bptime << ")" << endl;
                if (waiting) {
                    waiting = 0;
                    output.subpathtotal++;
                }
            }
            else {
                waiting = 1;
            }
        }
        if (flag) printBar();
    }
    void printUBPath(type_path path, Output& output, bool flag = 0) {
        if (flag) printBar();
        if (flag) cout << "Printing UB Path: " << endl;
        output.arctotal = 0;
        output.subpathtotal = 0;
        bool waiting = 1;
        for (int i = 1; i < path.size(); i++) {
            const TimedNode* node = path[i];
            if (path[i - 1]->nodeID != path[i]->nodeID) {
                TimedArcptr arc = { path[i - 1],path[i] };
                if (flag) cout << "arc cost: " << ttMapUB[arc] << endl;
                output.arctotal++;
                if (flag) cout << "at: (" << node->nodeID << ',' << node->time << "), in mangrove: (" << node->bpnode << ',' << node->bptime << ")" << endl;
                if (waiting) {
                    waiting = 0;
                    output.subpathtotal++;
                }
            }
            else {
                waiting = 1;
            }
        }
        if (flag) printBar();
    }
    void printCurrentMangroves() {
        printBar();
        cout << "Printing Current Mangroves:" << endl;
        for (auto mangroveset : mangroves) {
            for (auto mangrove : mangroveset) {
                cout << "(" << mangrove.bpnode << ',' << mangrove.bptime << "), resolved = " << mangrove.resolved << endl;
            }
        }
        printBar();
    }
    void printBPtoAdd(set<type_bp>& nextbps) {
        printBar();
        cout << "Printing BP to Add:" << endl;
        for (type_bp bp : nextbps) {
            cout << "(" << bp.first << ',' << bp.second << ")" << endl;
        }
        printBar();
    }
    void printBar() {
        cout << "====================================================================" << endl;
    }

    //Initializer
    Graph G;
    TEN(const int n, const int eT, const int gtype, const int ttype, const int seed,  const int sT = 0) : G(n, eT, sT), timednodes(n), ftimednodes(n, vector<set_timednode>(n)), btimednodes(n, vector<set_timednode>(n)), mangroves(n), mangroveMap(n) {
        /*number of nodes, end time, graph type, travel time type, seed, start time (optional)*/
        string filename = "Data/n" + to_string(n) + "T" + to_string(eT) + "gt" + to_string(gtype) + "tt" + to_string(ttype) + "s" + to_string(seed) + ".csv";
        cout << "opening:" << filename << endl;
        std::ifstream nodeFile(filename.c_str());
        if (!nodeFile.is_open()) {
            cerr << "error opening file" << endl;
        }
        string headers;
        getline(nodeFile, headers);
        string line;
        string startID;
        string endID;
        string time;
        while (getline(nodeFile, startID, ',')) {
            vector<double> times;
            getline(nodeFile, endID, ',');
            for (int i = sT; i < eT; i++) {
                getline(nodeFile, time, ',');
                times.push_back(stod(time, nullptr));
            }
            getline(nodeFile, time);
            times.push_back(stod(time, nullptr));
            G.addArc(stoi(startID, nullptr), stoi(endID, nullptr), times);
            //cout << "i:" << stoi(startID, nullptr) << "j:" << stoi(endID, nullptr) << "t:" << times[eT-sT] << endl;
        }
        if (G.checkFIFO()) {
            cout << "read successful!" << endl;
        }
        else {
            cerr << "Error: FIFO violated" << endl;
            cin.get();
        }
        return;
    }
};

void runTest(vector<int> ns, vector<int> endTs, vector<int> gtypes, vector<int> ttypes, vector<int> seeds, const int commonsT=0, bool md = 1, bool mtt = 1, bool mdenum = 1, bool mttenum = 1) {
    string mdfilename = "Results/experimentsMDMED.csv";
    string mttfilename = "Results/MTTBPMMED.csv";
    if (md) writeSOutputCSVHeader(mdfilename);
    if (mtt) writeSOutputCSVHeader(mttfilename);
    for (int n : ns) {
        for (int commoneT : endTs) {
            for (int gtype : gtypes) {
                for (int ttype : ttypes) {
                    SummaryOutput myMDOutput(mdfilename, n, commoneT - commonsT, gtype, ttype, seeds.size(), 1);
                    SummaryOutput myMTTOutput(mttfilename, n, commoneT - commonsT, gtype, ttype, seeds.size(), 0);
                    string filename = "Results/experimentsMDMEDlog.csv";
                    //string filename = "Results/n" + to_string(n) + "T" + to_string(commoneT - commonsT) + "g" + to_string(gtype) + "t" + to_string(ttype) + ".csv";
                    writeOutputCSVHeader(filename);
                    for (int seed : seeds) {
                        if (md) {
                            TEN my_ten(n, commoneT, gtype, ttype, seed);
                            Output output = my_ten.findMD();
                            output.writeOutputCSV(filename);
                            myMDOutput.updateSummaryOutput(output);
                        }
                        if (mdenum) {
                            TEN my_ten(n, commoneT, gtype, ttype, seed);
                            Output output = my_ten.findEnumMD();
                            output.writeOutputCSV(filename);
                            myMDOutput.updateSummaryOutput(output);
                        }
                        if (mtt) {
                            TEN my_ten(n, commoneT, gtype, ttype, seed);
                            Output output = my_ten.findMTT();
                            output.writeOutputCSV(filename);
                            myMTTOutput.updateSummaryOutput(output);
                        }
                        if (mttenum) {
                            TEN my_ten(n, commoneT, gtype, ttype, seed);
                            Output output = my_ten.findEnumMTT();
                            output.writeOutputCSV(filename);
                            myMTTOutput.updateSummaryOutput(output);
                        }
                    }
                    if (md) myMDOutput.calcSummaryOutput();
                    if (mtt) myMTTOutput.calcSummaryOutput();
                    if (md) myMDOutput.writeSOutputCSV();
                    if (mtt) myMTTOutput.writeSOutputCSV();
                }
            }
        }
    }
}


int main()
{
    //runTest({ 30 }, 40, { 1 }, { 1 }, { 1 });
    //runTest({ 30,50 }, 20, { 1,2,3 }, { 1,2 }, { 1,2,3,4,5,6,7,8,9,10 });
    //runTest({ 1000 }, { 20,40 }, { 1,2,3 }, { 1,2 }, { 1,2,3,4,5 }, 0, 1, 0, 1, 0);
    //runTest({ 1000 }, { 40 }, { 3 }, { 2 }, { 1,2,3,4,5 }, 0, 1, 0, 1, 0);
    //runTest({ 20 }, { 1000 }, { 1,2,3 }, { 1,2 }, { 1,2,3 }, 0, 0, 1, 0, 1);
    //runTest({ 100 }, { 20 }, { 2 }, { 1 }, { 1 }, 0, 0, 1, 0, 1);
    //runTest({ 30 }, { 20,60,100 }, { 2,3 }, { 1,2 }, { 1,2,3,4,5 }, 0, 0, 1, 0, 1);
    //runTest({ 50 }, { 40 }, { 1,2,3 }, { 1,2 }, { 1,2,3 }, 0, 0, 1, 0, 0);
    //runTest({ 50 }, { 20,40,60,80,100 }, { 1,2,3 }, { 1,2 }, { 1,2,3,4,5,6,7,8,9,10 }, 0, 1, 0, 1, 0);
    //runTest({ 60 }, { 50,100,150,200,250 }, { 1,2,3 }, { 1,2 }, { 1,2,3 }, 0, 0, 1, 0, 1);
    runTest({ 10000 }, { 40 }, { 1 }, { 1 }, { 1 }, 0, 1, 0, 1, 0);
    

    //TEN my_ten(30, 40, 1, 1, 1);
    //Output output = my_ten.findMTT();
    return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
