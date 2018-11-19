using namespace std;
class SolutionPts
{
public:
    SolutionPts(){}
    ~SolutionPts(){}

    void setvals(int i, double a, double b, double c){id = i; x = a; y=b; z=c;}
    int getid(){return id;}
    double  getx(){return x;}
    double gety(){return y;}
    void display(){cout<<"#"<<id<<" at ("<<x<<", "<<y<<", "<< z<< ")"<<endl;}
private:
    int id;
    double x;
    double y;
    double z;

};
