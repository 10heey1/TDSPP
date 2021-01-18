// TDSPP.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
//Don't create duplicate nodes, keep nodes sorted to avoid nodes being too similar. Map calls timed node instead of time.
//Careful what to set for btimes and ftimes for when they don't exist.
//Optimize outMapexLB and outMapexUB using vector allocation + resize?
//addInitialWaiting assumes time horizon large enough
//Currently findMTT half is uncommented
//Updating waiting arcs procedure not written yet
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <queue>
#include <set>
#include <assert.h>
using namespace std;
typedef pair<int,int> type_arc;
typedef pair<int, int> type_bp;

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
    map<type_arc, vector<vector<double>>> minttMap; //returns table of min travel times
    map<type_arc, vector<vector<int>>> minindttMap; //returns table of min travel times
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
        printf("Times of Nodes in FSPT\n");
        for (int i = 0; i < n; i++)
            printf("%d \t\t %f\n", i, fspt[i]);
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
        printf("Times of Nodes in BSPT\n");
        for (int j = 0; j < n; j++)
            printf("%d \t\t %f\n", j, bspt[j]);
        //printf("Arcs in BSPT\n");
        //for (int ind = 0; ind < arcs.size(); ind++)
        //    printf("%d \t\t %d\n", arcs[ind].first, arcs[ind].second);
        return bspt;
    }

    double minTT(const int startind, const int endind, double leftt, double rightt) {
        leftt = min((double)endT - startT, max(0.0, leftt - startT));
        rightt = min((double)endT - startT, max(0.0, rightt - startT));
        double ans = min(TT(startind, endind, leftt), TT(startind, endind, rightt));
        return min(ans, minttMap[{startind, endind}][ceil(leftt)][floor(rightt)]);
    }

    int minindTT(const int startind, const int endind, const int leftt, const int rightt) {
        return minindttMap[{startind, endind}][leftt][rightt];
    }

};

class TEN {
public:
    //Common
    //Minimum Duration
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
    void addABSPT(const int i, const int k) {
        cout << "Adding ABSPT:(" << i << ',' << k << ')' << endl;
        Abspt newabspt(i, k, G.n);
        vector<type_arc> arcs = {};
        double endTime = G.FSPT(i, k, arcs)[G.endN];
        arcs = {};
        newabspt.times = G.BSPT(G.endN, endTime, arcs);
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
            cout << "Updating ABSPT:(" << it->bpnode << ',' << it->bptime << ')' << endl;
            absptlbs.insert(it);
            absptubs.insert(it);
        }
        auto itbool = abspts.insert(newabspt);
        it = itbool.first;
        cout << "Inserting ABSPT:(" << it->bpnode << ',' << it->bptime << ')' << endl;
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
        cout << "Resolving ABSPT:(" << it->bpnode << ',' << it->bptime << ')' << endl;
        absptlbs.insert(it);
        absptubs.insert(it);
    }
    type_bp findBP(Abspt& curr) {
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
                bptime = G.minindTT(i, G.outMap[i][0], leftt, rightt);
                return { bpnode, bptime };
            }
        }
        return { bpnode,bptime };
    }
    Abspt findMD() {
        addABSPT(3, 5);
        addABSPT(0, 0);
        auto lbit = absptlbs.begin();
        auto ubit = absptubs.begin();
        cout << "current lower bound=" << (*lbit)->lb << endl;
        cout << "current upper bound=" << (*ubit)->ub << endl;
        while ((*lbit)->lb != (*ubit)->ub) {
            auto myabspt = **lbit;
            type_bp nextbp = findBP(myabspt);
            if (nextbp.first == -1) {
                resolveABSPT(myabspt);
            }
            else {
                addABSPT(nextbp.first, nextbp.second);
            }
            lbit = absptlbs.begin();
            ubit = absptubs.begin();
            cout << "current lower bound=" << (*lbit)->lb << endl;
            cout << "current upper bound=" << (*ubit)->ub << endl;
        }

        return **lbit;
    }
    //Minimum Travel Time
    struct Mangrove;
    struct TimedNode {
        double time = 0;
        int nodeID = -1;
        int bpnode = -1;
        int bptime = -1;
        TimedNode(int n, double t, int i, int k) {
            nodeID = n;
            time = t;
            bpnode = i;
            bptime = k;
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
    map<const TimedNode*, vector<const TimedNode*>> outMapexLB; //external (waiting) arcs for LB vector<const TimedNode*> is of size G.n
    map<const TimedNode*, vector<const TimedNode*>> outMapexUB; //external (waiting) arcs for UB vector<const TimedNode*> is of size G.n
    typedef pair<const TimedNode*, double> TnTT;
    struct TnTT_compare {
        bool operator() (const TnTT& lhs, const TnTT& rhs) const {
            return lhs.second > rhs.second;
        }
    };
    typedef vector<const TimedNode*> type_path;

    //Functions
    const Mangrove* addMangrove(const int bpnode, const int bptime) {
        cout << "Adding Mangrove:(" << bpnode << ',' << bptime << ')' << endl;
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
                TimedNode ftimednode(i, newmangrove.ftimes[i], bpnode, bptime);
                auto itbool = ftimednodes[i][bpnode].insert(ftimednode);
                auto it = itbool.first;
                auto currtimednode = &(*it);
                newmangrove.fnodes[i] = currtimednode;
                vector<const TimedNode*> outVec(G.n, nullptr);
                outMapexLB[currtimednode] = outVec;
                outMapexUB[currtimednode] = outVec;
                if (i == G.startN && origin != nullptr) {
                    outMapexLB[origin].push_back(currtimednode);
                    outMapexUB[origin].push_back(currtimednode);
                }
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
                TimedNode btimednode(i, newmangrove.btimes[i], bpnode, bptime);
                auto itbool = btimednodes[i][bpnode].insert(btimednode);
                auto it = itbool.first;
                auto currtimednode = &(*it);
                newmangrove.bnodes[i] = currtimednode;
                vector<const TimedNode*> outVec(G.n + 1, nullptr);
                outMapexLB[currtimednode] = outVec;
                outMapexUB[currtimednode] = outVec;
                if (i == G.endN && destination != nullptr) {
                    outMapexLB[currtimednode][G.n] = destination;
                    outMapexUB[currtimednode][G.n] = destination;
                }
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
            const TimedNode* i = pq.top().first;
            if (i == destination) break;
            //Care double comparison
            if (pq.top().second != dp[i]) {
                //cout << "(not) popping:" << i << ',' << pq.top().first << endl;
                pq.pop();
                continue;
            }
            //cout << "popping:" << i << ',' << fspt[i] << endl;
            pq.pop();
            for (const TimedNode* j : outMapLB[i]) {
                //cout << i << j << startts[i] << endts[i] << endl;
                double weight = ttMapLB[{i, j}];
                auto it = dp.find(j);
                if (it == dp.end()) dp[j] = INT_MAX;
                if (dp[j] > dp[i] + weight) {
                    dp[j] = dp[i] + weight;
                    pq.push({ j,dp[j] });
                    pred[j] = i;
                    //cout << "pushing:" << j << ',' << fspt[j] << endl;
                }
            }
            for (const TimedNode* j : outMapexLB[i]) {
                if (j == nullptr) continue;
                //cout << i << j << startts[i] << endts[i] << endl;
                auto it = dp.find(j);
                if (it == dp.end()) dp[j] = INT_MAX;
                if (dp[j] > dp[i]) {
                    dp[j] = dp[i];
                    pq.push({ j,dp[j] });
                    pred[j] = i;
                    //cout << "pushing:" << j << ',' << fspt[j] << endl;
                }
            }
            //printf("Times of Nodes in FSPT\n");
            //for (int i = 0; i < n; i++)
            //    printf("%d \t\t %f\n", i, fspt[i]);
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
    double findUB() {
        map<const TimedNode*, double> dp;
        priority_queue<TnTT, vector<TnTT>, TnTT_compare> pq;
        pq.push({ origin,0 });
        dp[origin] = 0;
        while (!pq.empty()) {
            const TimedNode* i = pq.top().first;
            if (i == destination) break;
            //Care double comparison
            if (pq.top().second != dp[i]) {
                //cout << "(not) popping:" << i << ',' << pq.top().first << endl;
                pq.pop();
                continue;
            }
            //cout << "popping:" << i << ',' << fspt[i] << endl;
            pq.pop();
            for (const TimedNode* j : outMapUB[i]) {
                //cout << i << j << startts[i] << endts[i] << endl;
                double weight = ttMapUB[{i, j}];
                auto it = dp.find(j);
                if (it == dp.end()) dp[j] = INT_MAX;
                if (dp[j] > dp[i] + weight) {
                    dp[j] = dp[i] + weight;
                    pq.push({ j,dp[j] });
                    //cout << "pushing:" << j << ',' << fspt[j] << endl;
                }
            }
            for (const TimedNode* j : outMapexUB[i]) {
                if (j == nullptr) continue;
                //cout << i << j << startts[i] << endts[i] << endl;
                auto it = dp.find(j);
                if (it == dp.end()) dp[j] = INT_MAX;
                if (dp[j] > dp[i]) {
                    dp[j] = dp[i];
                    pq.push({ j,dp[j] });
                    //cout << "pushing:" << j << ',' << fspt[j] << endl;
                }
            }
            //printf("Times of Nodes in FSPT\n");
            //for (int i = 0; i < n; i++)
            //    printf("%d \t\t %f\n", i, fspt[i]);
        }
        return dp[destination] - dp[origin];
    }
    type_bp findBP(const Mangrove* curr) {
        int bpnode = -1;
        int bptime = -1;
        auto it = mangroves[curr->bpnode].find(*curr);
        it++;
        if (it == mangroves[curr->bpnode].end()) return { bpnode,bptime };
        for (int i = 0; i < G.n; i++) {
            if (G.outMap[i].empty()) continue;
            //backward
            double dleftt = curr->btimes[i];
            double drightt = it->btimes[i];
            int leftt, rightt;
            if (dleftt == (int)dleftt) leftt = dleftt + 1;
            else leftt = ceil(dleftt);
            if (drightt == (int)drightt) rightt = drightt - 1;
            else rightt = floor(drightt);
            if (leftt > rightt) continue;
            else {
                bpnode = i;
                bptime = G.minindTT(i, G.outMap[i][0], leftt, rightt);
                return { bpnode, bptime };
            }
            //forward
            dleftt = curr->ftimes[i];
            drightt = it->ftimes[i];
            if (dleftt == (int)dleftt) leftt = dleftt + 1;
            else leftt = ceil(dleftt);
            if (drightt == (int)drightt) rightt = drightt - 1;
            else rightt = floor(drightt);
            if (leftt > rightt) continue;
            else {
                bpnode = i;
                bptime = G.minindTT(i, G.outMap[i][0], leftt, rightt);
                return { bpnode, bptime };
            }
        }
        return { bpnode, bptime };
    }
    vector<type_bp> findBP(type_path& path) {
        vector<type_bp> bps;
        for (const TimedNode* timednode : path) {
            type_bp bp = findBP(mangroveMap[timednode->bpnode][timednode->bptime]);
            if (bp.first != -1) bps.push_back(bp);
        }
        return bps;
    }
    pair<type_path, double> findMTT() {
        const Mangrove* firstmangrove = addMangrove(G.startN, G.startT);
        origin = firstmangrove->fnodes[G.startN];
        const Mangrove* lastmangrove = addMangrove(G.endN, G.endT);
        destination = lastmangrove->bnodes[G.endN];
        for (int i = 0; i < G.n; i++) {
            if (i == G.endN) continue;
            addMangrove(i, floor(lastmangrove->btimes[i]));
        }
        for (int i = 0; i < G.n; i++) {
            if (i == G.startN) continue;
            addMangrove(i, ceil(firstmangrove->ftimes[i]));
        }
        addInitialWaiting();
        pair<type_path, double> pathLB = findLB();
        double UB = findUB();
        /*
        while (UB != pathLB.second) {
            vector<type_bp> nextbps = findBP(pathLB.first);
            for (type_bp nextbp : nextbps) {
                addABSPT(nextbp.first, nextbp.second);
            }
            pair<type_path, double> pathLB = findLB();
            double UB = findUB();
            cout << "current lower bound=" << pathLB.second << endl;
            cout << "current upper bound=" << UB << endl;
        }
        */
        return pathLB;
    }
    void addInitialWaiting() {
        /*
        Assumptions are two mangroves per bpnode, time horizon large enough
        */
        for (int i = 0; i < G.n; i++) {
            auto it = mangroveMap[i].begin();
            const Mangrove* firstmanptr = ((*it).second);
            const Mangrove* lastmanptr = ((*next(it,1)).second);
            for (int node = 0; node < G.n; node++) {
                //Waiting from i-FSPT to i-FSPT
                const TimedNode* first = firstmanptr->fnodes[node];
                const TimedNode* last = lastmanptr->fnodes[node];
                if (first == nullptr || last == nullptr) continue;
                outMapexLB[first][i] = last;
                outMapexUB[first][i] = last;
            }
            for (int j = 0; j < G.n; j++) {
                if (i == j) continue; //covered already
                auto otherit = mangroveMap[j].begin();
                const Mangrove* firstothermanptr = ((*otherit).second);
                const Mangrove* lastothermanptr = ((*next(otherit, 1)).second);
                for (int node = 0; node < G.n; node++) {
                    //Waiting from i-FSPT to j-BSPT
                    const TimedNode* ifirst = firstmanptr->fnodes[node];
                    const TimedNode* ilast = lastmanptr->fnodes[node];
                    const TimedNode* jfirst = firstothermanptr->bnodes[node];
                    const TimedNode* jlast = lastothermanptr->bnodes[node];
                    if (ifirst == nullptr) {
                        assert(ilast == nullptr);
                        continue;
                    }
                    if (jfirst == nullptr) {
                        assert(jlast == nullptr);
                        continue;
                    }
                    assert(ilast != nullptr);
                    assert(jlast != nullptr);
                    if (ilast->time < jlast->time) {
                        outMapexLB[ilast][j] = jfirst;
                        if (ilast->time < jfirst->time) {
                            outMapexUB[ilast][j] = jfirst;
                        }
                        else {
                            outMapexUB[ilast][j] = jlast;
                            if (ifirst->time < jfirst->time) {
                                outMapexUB[ifirst][j] = jfirst;
                            }
                        }
                    }
                    else {
                        outMapexLB[ilast][j] = jlast;
                        if (ifirst->time < jlast->time) {
                            outMapexLB[ifirst][j] = jfirst;
                            if (ifirst->time < jfirst->time) {
                                outMapexUB[ifirst][j] = jfirst;
                            }
                            else {
                                outMapexUB[ifirst][j] = jlast;
                            }
                        }
                    }
                }
            }
        }
    }
    void addiiWaiting(Mangrove& newmangrove, const Mangrove* prevmanptr, const Mangrove* nextmanptr) {
        int i = newmangrove.bpnode;
        for (int node = 0; node < G.n; node++) {
            if (newmangrove.fnodes[node] != nullptr) {
                assert(prevmanptr->fnodes[node] != nullptr);
                assert(nextmanptr->fnodes[node] != nullptr);
                outMapexLB[prevmanptr->fnodes[node]][i] = newmangrove.fnodes[node];
                outMapexLB[newmangrove.fnodes[node]][i] = prevmanptr->fnodes[node];
            }
        }
    }
    void addijWaiting(Mangrove& newmangrove, const Mangrove* prevmanptr, const Mangrove* nextmanptr) {
        int i = newmangrove.bpnode;
        for (int node = 0; node < G.n; node++) {
            const TimedNode* currptr = newmangrove.fnodes[node];
            if (currptr != nullptr) {
                for (int j = 0; j < G.n; j++) {
                    if (j == i) continue;
                    auto nextjit = btimednodes[node][j].upper_bound(*currptr);
                    const TimedNode* nextjptr = &(*nextjit);
                    assert(nextjptr != nullptr);
                    auto prevjit = prev(nextjit, 1);
                    const TimedNode* prevjptr = &(*prevjit);
                    assert(prevjptr != nullptr);
                    if (mangroveMap[prevjptr->bpnode][prevjptr->bptime]->resolved) {
                        outMapexLB[currptr][j] = nextjptr;
                    }
                    else {
                        outMapexLB[currptr][j] = prevjptr;
                    }
                }
            }
        }
    }
    void printAllWaitingArcs() {
        for (auto kv : outMapexLB) {
            auto i = kv.first;
            cout << "from: (" << i->nodeID << ',' << i->time << "), in mangrove: (" << i->bpnode << ',' << i->bptime << ")" << endl;
            for (auto j : kv.second) {
                if (j == nullptr) continue;
                cout << "to: (" << j->nodeID << ',' << j->time << "), in mangrove: (" << j->bpnode << ',' << j->bptime << ")" << endl;
            }
        }
    }
    void printLBPath(type_path path) {
        for (int i = 0; i < path.size(); i++) {
            const TimedNode* node = path[i];
            if (i != 0) {
                TimedArcptr arc = { path[i - 1],path[i] };
                cout << "arc cost: " << ttMapLB[arc] << endl;
                double UBcost = (path[i - 1]->nodeID == path[i]->nodeID) ? 0 : path[i]->time - path[i - 1]->time;
                cout << "UB arc cost: " << UBcost << endl;
            }
            cout << "at: (" << node->nodeID << ',' << node->time << "), in mangrove: (" << node->bpnode << ',' << node->bptime << ")" << endl;
        }
    }
    //Initializer
    Graph G;
    TEN(const int n, const int eT, const int gtype, const int seed, const int ttype, const int sT = 0) : G(n, eT, sT), timednodes(n), ftimednodes(n, vector<set_timednode>(n)), btimednodes(n, vector<set_timednode>(n)), mangroves(n), mangroveMap(n) {
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
        cout << "read successful!" << endl;
        return;
    }
};

int main()
{
    TEN my_ten(4, 7, 1, 1, 1);
    vector<pair<int, int>> arcs = {};
    vector<double> times = my_ten.G.FSPT(0, 0, arcs);
    TEN::Abspt optMD = my_ten.findMD();
    cout << "Optimal MD starting time:" << optMD.times[my_ten.G.startN] << endl;
    cout << "Optimal MD val:" << optMD.lb << endl;
    cout << "Optimal MD breakpoint: (" << optMD.bpnode << ',' << optMD.bptime << ')' << endl;
    auto pathval = my_ten.findMTT();
    cout << pathval.second << endl;
    my_ten.printAllWaitingArcs();
    my_ten.printLBPath(pathval.first);
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
