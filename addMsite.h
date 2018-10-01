// create a new subclase that inherits the gmx reader class
#include <string.h>
#include <complex.h>
#include <gmx_reader.h>

#define PSTR        "||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PWID        50


#ifndef HBondDistribution_H
#define HBondDistribution_H
class model: public gmx_reader
{
    public:
        // class variables
        string  outf="trajMsite.xtc";           // name for output files
        string  watermodel="tip4p2005";         // name of model
        const int OW=0, HW1=1, HW2=2, MW=3;     // integers to access water atoms
        rvec *xnew;

    
        // Default constructor and destructor
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        model( string _inpf_ );
        ~model();

        // functions
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        void  addMsite();
};
#endif
