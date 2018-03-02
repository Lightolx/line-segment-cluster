#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
//#include <boost/filesystem.hpp>
#include "eigen3/Eigen/Eigen"
#include "clustering.h"

// compute the min distance from a point to a line segment
double pt2frag(std::pair<Eigen::Vector3d, Eigen::Vector3d> frag,
               Eigen::Vector3d c)
{
    Eigen::Vector3d a = frag.first;
    Eigen::Vector3d b = frag.second;
    Eigen::Vector3d ac = c-a;
    Eigen::Vector3d ab = b-a;
    Eigen::Vector3d Iab = ab/ab.norm();
    Eigen::Vector3d ad = ac.dot(Iab)*Iab;
    Eigen::Vector3d d = a+ad;
    Eigen::Vector3d da = -ad;
    Eigen::Vector3d db = b-d;
    if (da.dot(db) < 0)
    {
        return ac.cross(Iab).norm();
    }

    double lenac = ac.norm();
    double lenbc = (b-c).norm();

    return lenac < lenbc ? lenac : lenbc;

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
}

int main()
{
    // ***************read in lines************************//
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

    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> lines;

    int numPoints = points.size();

    for (int i = 0; i < numPoints; i = i+2)
    {
        std::pair<Eigen::Vector3d, Eigen::Vector3d> line(points[i], points[i+1]);
        lines.push_back(line);
    }

    // get relation edge matrix A
    int num_lines = lines.size();
    double distTrd = 2;
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

    int tmp = A.size();
    int n = 0;

    // return a conjoint graph
    cluGraph graph(num_lines);

    std::vector<std::pair<int, int>>::iterator iter = A.begin();
    for (; iter != A.end(); iter++)
    {
        std::pair<int, int> edge = *iter;
        int a = graph.findCluID(edge.first);
        int b = graph.findCluID(edge.second);

        if (a != b)
        {
            graph.join(a, b);
            n++;
        }
    }

    // outout the conjoint graph
    std::map<int, std::list<int> > cluster2lines;

    for (int i = 0; i < num_lines; ++i)
    {
        int clusterID = graph.findCluID(i);
        cluster2lines[clusterID].push_back(i);
    }


    std::ofstream fout("/media/psf/Home/Desktop/cluster.txt");

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