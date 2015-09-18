#include <iostream>
using namespace std;
#include <stdlib.h>
#include <string.h>

#define AAA 100
class test
{
  private:
    int L;
    double *p;
        
  public:
    void set_L(int Lp=32)
    {
        L = Lp;
    }
    
    void print_L() { cout << "L=" << L << endl;}
    
    void set_p()
    {
        p = new double[L];
        for(int i=0; i<L; i++)
        {
            p[i] = i;
            cout << "p[i]" << p[i] << endl;
        }
    }
    
    void print_p()
    {
        int i=0;
        for(i=0; i<L; i++)
        {
            cout << " " << p[i] << endl;
        }
    }
};

int main()
{
    test a;
    a.set_L(3);
    //a.print_L();
    //a.set_p();
    //a.print_p();
    
    char name[8];
    double b = 6272.1;
    
    cout << b/AAA << endl;
    strncpy(name, "123", 8);
    //cout << name << endl;
    return(0);
}
