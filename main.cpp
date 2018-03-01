#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
//#include <boost/filesystem.hpp>
#include "eigen3/Eigen/Eigen"
#include "clustering.h"

// distance from a point to a line
double frag2grag(std::pair<Eigen::Vector3d, Eigen::Vector3d> frag,
               Eigen::Vector3d a, Eigen::Vector3d b)
{
    /*
    double x1 = frag.first[0];
    double x2 = frag.second[0];
    double y1 = frag.first[1];
    double y2 = frag.second[1];
    double z1 = frag.first[2];
    double z2 = frag.second[2];
    double A = z2 - z1;
    double B = x1 - x2;
    double C = x2*z1 - x1*z2;

    double d1 = (A*a[0] + B*a[2] + C)/sqrt(A*A + B*B);
    double d2 = (A*b[0] + B*b[2] + C)/sqrt(A*A + B*B);
    double d3 = (A*(a[0]+b[0])/2 + B*(a[2]+b[2])/2 + C)/sqrt(A*A + B*B);
    return fabs(d1+d2+2*d3)/sqrt(pow(a[0]-b[0],2) + pow(a[2]-b[2],2));
     */

    /*
    int sampleDense = 10;
    Eigen::Vector3d delta = (b - a)/sampleDense;
    std::vector<Eigen::Vector3d> pts;
    pts.reserve(sampleDense+1);
    Eigen::Vector3d pt(0,0,0);
    for (int i = 0; i < sampleDense; ++i)
    {
        pt = a + delta*i;
        pts.push_back(pt);
    }

    Eigen::Vector3d cd = frag.second - frag.first;
    std::vector<double> dists;
    dists.reserve(pts.size());
    Eigen::Vector3d c = frag.first;
    Eigen::Vector3d d = frag.second;

    for (int i = 0; i < pts.size(); ++i)
    {
        pt = pts[i];
        Eigen::Vector3d cp = pt - frag.first;
        double dist = cp.cross(cd).norm();
        dists.push_back(dist);
    }

    std::vector<double>::iterator iter = std::min_element(dists.begin(), dists.end());
    double candi = *iter;
    int idx = std::distance(dists.begin(), iter);
    Eigen::Vector3d np = pts[idx];
    Eigen::Vector3d tmpc = np-c;
    double tmp = cd.dot(np-c)/cd.norm();
    Eigen::Vector3d vdp = cd.dot(np - c)/cd.norm()*cd/cd.norm();
    Eigen::Vector3d dp = c + vdp;
    if ((dp - frag.first).dot(dp - frag.second) < 0)
    {
        return *iter;
    }
    else
    {
        return 10;
    }

    // add angle into account
//    return fabs(dist1 + dist2)/sqrt(pow(a[0]-b[0],2) + pow(a[2]-b[2],2));
//    return minDist;
    */


    Eigen::Vector3d c = frag.first;
    Eigen::Vector3d d = frag.second;
    Eigen::Vector3d ca = a - c;
    Eigen::Vector3d cb = b - c;
    Eigen::Vector3d cd = d - c;

    double dist1 = ca.cross(cd).norm();
    double dist2 = cb.cross(cd).norm();
    return fabs(dist1 + dist2)/sqrt(pow(a[0]-b[0],2) + pow(a[2]-b[2],2));

}

//double pt2frag(std::pair<Eigen::Vector3d, Eigen::Vector3d> frag,
//               Eigen::Vector3d a)
//{
//    Eigen::Vector3d c = frag.first;
//    Eigen::Vector3d d = frag.second;
//    double dist1 = sqrt(pow(a[0]-c[0],2) + pow(a[2]-c[2],2));
//    double dist2 = sqrt(pow(a[0]-d[0],2) + pow(a[2]-d[2],2));
//    // todo:: add angle into account
//    return dist1 < dist2 ? dist1 : dist2;
//}

double pt2frag(std::pair<Eigen::Vector3d, Eigen::Vector3d> frag,
               Eigen::Vector3d c)
{
    Eigen::Vector3d a = frag.first;
    Eigen::Vector3d b = frag.second;
    Eigen::Vector3d ac = c-a;
    Eigen::Vector3d ab = b-a;
    Eigen::Vector3d Iab = ab/ab.norm();
    return ac.cross(Iab).norm();
//    Eigen::Vector3d n = ac.cross(ab);
//    return sqrt(pow(n[0], 2) + pow(n[1], 2));
}

double pt2Group(std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> lines,
                Eigen::Vector3d a)
{
    int num_segment = lines.size();
    std::vector<double> dist(num_segment,100);

    for (int i = 0; i < lines.size(); ++i)
    {
        std::pair<Eigen::Vector3d, Eigen::Vector3d> frag = lines[i];
        dist[i] = pt2frag(frag, a);
    }

    return *(std::min_element(dist.begin(), dist.end()));
}


// judege the similarity between two line segments
bool line2line(const std::pair<Eigen::Vector3d, Eigen::Vector3d> lines1,
                 const std::pair<Eigen::Vector3d, Eigen::Vector3d> lines2,
                 const double distTrd, const double angleTrd)
{
    Eigen::Vector2d a = Eigen::Vector2d(lines1.first[0], lines1.first[2]);
    Eigen::Vector2d b = Eigen::Vector2d(lines1.second[0], lines1.second[2]);
    Eigen::Vector2d c = Eigen::Vector2d(lines2.first[0], lines2.first[2]);
    Eigen::Vector2d d = Eigen::Vector2d(lines2.second[0], lines2.second[2]);

    double x1 = a[0];double y1 = a[1];
    double x2 = b[0];double y2 = b[1];
    double x3 = c[0];double y3 = c[1];
    double x4 = d[0];double y4 = d[1];

    double D = (y1-y2)*(x4-x3) - (y3-y4)*(x2-x1);
    double D1 = (x2*y1 - x1*y2)*(x4-x3) - (x4*y3 - x3*y4)*(x2-x1);
    double D2 = (y1-y2)*(x4*y3-x3*y4) - (y3-y4)*(x2*y1-x1*y2);
    double crossX = D1/D;
    double crossZ = D2/D;
    Eigen::Vector2d p = Eigen::Vector2d(crossX, crossZ);

    Eigen::Vector2d pa = a - p;
    Eigen::Vector2d pb = b - p;
    Eigen::Vector2d pc = c - p;
    Eigen::Vector2d pd = d - p;

    bool onAB = true;
    bool onCD = true;

    if (pa.dot(pb) > 0)
    {
        onAB = false;
    }

    if (pc.dot(pd) > 0)
    {
        onCD = false;
    }

    bool isCross = false;
    double dist = 100;

    if (onAB&&onCD)       // the two line segments do cross
    {
        isCross = true;
        dist = 0;
    }
    else                 // the two line segments does not cross
    {
        // find the distance of a to cd, b to cd, c to ab and d to ab, return the smallest one
        double distA = pt2frag(lines2, lines1.first);
        double distB = pt2frag(lines2, lines1.second);
        double distC = pt2frag(lines1, lines2.first);
        double distD = pt2frag(lines1, lines2.second);
        double minDist1 = distA < distB ? distA : distB;
        double minDist2 = distC < distD ? distC : distD;
        isCross = false;
        dist = minDist1 < minDist2 ? minDist1 : minDist2;
    }

    double angle = 2;
    Eigen::Vector2d ab = b-a;
    Eigen::Vector2d cd = d-c;
    angle = acos(fabs(ab.dot(cd))/(ab.norm()*cd.norm()));

    if (dist < distTrd && angle < angleTrd)
    {
        return true;
    }

    return false;


    /*
    // note the foot point of c to ab be point e, foot point of d to ab be point f, find them
    Eigen::Vector3d a = lines1.first;
    Eigen::Vector3d b = lines1.second;
    Eigen::Vector3d c = lines2.first;
    Eigen::Vector3d d = lines2.second;
    Eigen::Vector3d ab = b - a;
    Eigen::Vector3d ac = c - a;
    Eigen::Vector3d ad = d - a;
    double lenab = ab.norm();               // length of vector ab
    Eigen::Vector3d Iab = ab/lenab;         // Identity vector of vector ab
    Eigen::Vector3d ae = ac.cross(Iab)*Iab;
    Eigen::Vector3d af = ad.cross(Iab)*Iab;
    Eigen::Vector3d e = a + ae;
    Eigen::Vector3d f = a + af;

    bool eOn = true;
    bool fOn = true;
    bool isCross = false;
    double dist = 100;

    if (ae.cross(ab) < 0)  // e is along the vector ab, not ba
    {
        eOn = false;
    }

    if ((f-b).cross(a-b) < 0)
    {
        fOn = false;
    }

    Eigen::Vector3d ec = c-e;
    Eigen::Vector3d fd = d-f;
    if (ec.cross(fd) > 0)          // c,d on the same side
    {
        isCross = false;
        double lenec = ec.norm();
        double lenfd = fd.norm();
        dist = lenec < lenfd ? lenec : lenfd;
    }

    if (eOn&&fOn)
    {

    }
    else if (eOn&&!fOn)
    {

    }
    */
}

int main()
{
    std::ifstream fin("lines.txt");
    std::string ptline;
    double x, y, z;
    std::vector<Eigen::Vector3d> points;

    while (getline(fin, ptline))
    {
        std::stringstream ss(ptline);
        ss >> x >> y >> z;
        points.push_back(Eigen::Vector3d(x,y,z));
    }

//    std::vector<Eigen::Vector3d>::iterator iter = coordinates.begin();
//    for (; iter != coordinates.end(); iter++)
//    {
//        std::cout << *iter << std::endl;
//    }

    std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>> categories;
    int pointSize = points.size()/2;
    categories.resize(points.size()/2);  // points.size()/2 categories at most
    int cate = 1;
    Eigen::Vector3d a,b;

    double avaDistTrh = 1.5;
    double dist1Trh = 3;
    double dist2Trh = 3;
    int tmp = points.size();

    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> lines;

    for (int i = 0; i < points.size(); i=i+2)
    {
        a = points[i];
        b = points[i+1];
        lines.push_back(std::pair<Eigen::Vector3d, Eigen::Vector3d>(a,b));
        Eigen::Vector3d pt(-243.354000000000, 9.36483000000000, -271.059000000000);
        Eigen::Vector3d offset = pt-a;

        if (offset.norm() < 0.01)
        {
            int tem = 1;
        }
        if (i == 0)
        {
            categories[cate-1].push_back(std::pair<Eigen::Vector3d, Eigen::Vector3d>(a, b));
        }

        if (i > 0)
        {
            bool findAffli = false;
            std::vector<int> cluster;
            cluster.reserve(cate);

            for (int j = 0; j < cate; j++)
            {
                int num_segment = categories[j].size();
                std::vector<double> avaDist(num_segment,100);
//                std::vector<double> dist1(num_segment,100);
//                std::vector<double> dist2(num_segment,100);

                for (int k = 0; k < num_segment; k++)
                {
                    std::pair<Eigen::Vector3d, Eigen::Vector3d> frag = categories[j][k];
                    avaDist[k] = frag2grag(frag, a, b);
//                    dist1[k] = pt2frag(frag, a);
//                    dist2[k] = pt2frag(frag, b);
                }

                double minAvaDist =
                *(std::min_element(avaDist.begin(), avaDist.end()));
                double minDist1 = pt2Group(categories[j], a);
                double minDist2 = pt2Group(categories[j], b);

                if (minAvaDist < avaDistTrh && minDist1 < dist1Trh && minDist2 < dist2Trh)
//                if (minAvaDist < avaDistTrh)
                {
                    cluster.push_back(j);
//                    categories[j].push_back(std::pair<Eigen::Vector3d, Eigen::Vector3d>(a,b));
                    findAffli = true;
//                    break;
                }
            }



            if (findAffli)
            {
                int clusterSize = cluster.size();
                int ID = cluster[0];

                if (clusterSize == 1)
                {

                    categories[ID].push_back(std::pair<Eigen::Vector3d, Eigen::Vector3d>(a,b));
                }

                else
                {
                    cate -= clusterSize-1;

                    for (int k = 1; k < clusterSize; k++)
                    {
                        categories[ID].insert(categories[ID].end(),
                        categories[cluster[k]].begin(), categories[cluster[k]].begin());
                        categories.erase(categories.begin()+cluster[k]-(k-1));
                    }
                }
            }
            else
            {
                cate++;
                categories[cate-1].push_back(std::pair<Eigen::Vector3d, Eigen::Vector3d>(a,b));
            }

            /*
            std::vector<double> bias(cate, 100);  // innitial to be 100 which is bigger than threshold
            std::cout << "for points " << i << " and " << i+1 << " dist are\n";

            for (int j = 0; j < cate; j++)
            {
                int num_segment = categories[j].size();
                std::vector<double> bia(num_segment, 100);
                for (int k = num_segment; k > 0 && k > num_segment-3; k--)
                {
                    std::pair<Eigen::Vector3d, Eigen::Vector3d> frag = categories[j][k-1];
                    bia.push_back(pt2line(frag, a, b));
                }
                bias[j] = *(std::min_element(bia.begin(), bia.end()));
                std::cout << bias[j] << std::endl;
//                categories[j].push_back(std::pair<Eigen::Vector3d, Eigen::Vector3d>(a,b));
//                std::cout << similarities[j] << std::endl;
            }

            std::vector<double>::iterator minBaisIter = std::min_element(bias.begin(), bias.end());
            double minBais = *minBaisIter;

            if (minBais > 1)
            {
                cate++;
                categories[cate-1].push_back(std::pair<Eigen::Vector3d, Eigen::Vector3d>(a,b));
            }
            else
            {
                int idx = std::distance(bias.begin(), minBaisIter);
                categories[idx].push_back(std::pair<Eigen::Vector3d, Eigen::Vector3d>(a,b));
            }
             */


        }

    }

    int num_lines = lines.size();
    double distTrd = 1;
    double angleTrd = 10;
    std::vector<std::pair<int, int>> A;

    for (int i = 0; i < num_lines; ++i)
    {
        std::pair<Eigen::Vector3d, Eigen::Vector3d> ab = lines[i];

        for (int j = i+1; j < num_lines; ++j)
        {
            std::pair<Eigen::Vector3d, Eigen::Vector3d> cd = lines[j];
            bool isConnect = line2line(ab, cd, distTrd, angleTrd);

            if (isConnect)
            {
                A.push_back(std::pair<int, int>(i, j));
                A.push_back(std::pair<int, int>(j, i));
            }
        }
    }

    cluGraph graph(A.size());

    std::vector<std::pair<int, int>>::iterator iter = A.begin();
    for (; iter != A.end(); iter++)
    {
        std::pair<int, int> edge = *iter;
        int a = graph.findCluID(edge.first);
        int b = graph.findCluID(edge.second);

        if (a != b)
        {
            graph.join(a, b);
        }
    }

    std::map<int, std::list<int> > cluster2lines;

    for (int i = 0; i < num_lines; ++i)
    {
        int clusterID = graph.findCluID(i);

        if (cluster2lines.find(clusterID) == cluster2lines.end())
        {
            cluster2lines[clusterID].push_back(i);
        }
    }


    std::ofstream fout("/media/psf/Home/Desktop/cluster.txt");

//    for (int i = 0; i < cate; i++)
//    {
//        for (int j = 0; j <  categories[i].size(); j++)
//        {
//            Eigen::Vector3d a = categories[i][j].first;
//            Eigen::Vector3d b = categories[i][j].second;
//            fout << a[0] << " " << a[1] << " " << a[2] << " " << i << std::endl;
//            fout << b[0] << " " << b[1] << " " << b[2] << " " << i << std::endl;
//        }
//    }

    std::map<int, std::list<int> >::iterator iter1 = cluster2lines.begin();
    for (int i = 0; iter1 != cluster2lines.end(); iter1++, i++)
    {
        std::list<int> linesID = (*iter1).second;

        std::list<int>::iterator iter2;
        for (iter2 = linesID.begin(); iter2!= linesID.end(); iter2++)
        {
            Eigen::Vector3d a = lines[*iter2].first;
            Eigen::Vector3d b = lines[*iter2].second;
            fout << a[0] << " " << a[1] << " " << a[2] << " " << i << std::endl;
            fout << b[0] << " " << b[1] << " " << b[2] << " " << i << std::endl;
        }
    }

}