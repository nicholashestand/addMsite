#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <xdrfile.h>
#include <xdrfile_xtc.h>
#include <gmx_reader.h>
#include <fftw3.h>
#include <omp.h>
#include <pthread.h>
#include "addMsite.h"

#define PI 3.14159265359 

using namespace std;

// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
model::model( string _inpf_ ) : gmx_reader::gmx_reader( _inpf_ )
// Default constructor
{

    // set userparams from input file
    for ( int i = 0; i < nuParams; i ++ )
    {
        if ( uParams[i] == "outf" )  outf  = uValues[i];
        if ( uParams[i] == "watermodel" ) watermodel = uValues[i];
    }

    cout << "Set outf to: "         << outf         << endl;
    cout << "Set watermodel to: "        << watermodel   << endl;

    if ( watermodel != "tip4p2005" and watermodel != "tip4p" and watermodel != "e3b2" and watermodel != "e3b3" ){
        cout << "Warning:: watermodel: " << watermodel << " unknown.\nOnly 'tip4p2005', 'tip4p', 'e3b3' and 'e3b2' are acceptable values. Aborting!" << endl;
        exit(EXIT_FAILURE);
    }

    // allocate new array for xnew
    xnew = new rvec[ nmol * 4 ];
}

model::~model()
// Default Destructor
{
    delete [] xnew;
}

void model::addMsite()
// set msite OH distance to tip4p geometry
{
    float OMlen;

    if ( watermodel == "tip4p2005" or watermodel == "e3b3" ) OMlen = 0.01546;
    else OMlen = 0.01500;
    #pragma omp parallel for
    for ( int mol = 0; mol < nmol; mol ++ ){
        int   i;
        float oh1_vec[3], oh2_vec[3], r, om_vec[3];

        // the OH unit vectors
        for ( i = 0; i < 3; i ++ ) {
            oh1_vec[i] = x[ mol*natoms_mol + HW1 ][i] - x[ mol*natoms_mol + OW ][i];
            oh2_vec[i] = x[ mol*natoms_mol + HW2 ][i] - x[ mol*natoms_mol + OW ][i];
        }
        minImage( oh1_vec );
        minImage( oh2_vec );

        // get OM unit vec
        for ( i = 0; i < 3; i ++ ) om_vec[i] = oh1_vec[i] + oh2_vec[i];
        minImage( om_vec );
        r = mag3( om_vec );
        for ( i = 0; i < 3; i ++ ) om_vec[i] /= r;

        // set the m site based on the unit vector
        for ( i = 0; i < 3; i ++ ){
            xnew[ mol*4 + OW ][i]  = x[ mol*natoms_mol + OW  ][i];
            xnew[ mol*4 + HW1 ][i] = x[ mol*natoms_mol + HW1 ][i];
            xnew[ mol*4 + HW2 ][i] = x[ mol*natoms_mol + HW2 ][i];
            xnew[ mol*4 + MW ][i]  = x[ mol*natoms_mol + OW  ][i]+ OMlen * om_vec[i];
        }
    }
}

// ************************************************************************ 
int main( int argc, char* argv[] )
{

    int     currentSample, frameno, tcfpoint;
    float   currentTime, tcfTime;

    // Check program input
    if ( argc != 2 ){
        printf("Program expects only one argument, which is the name of \n\
                an input file containing the details of the analysis.\nAborting...\n");
        exit(EXIT_FAILURE);
    }

    // get filename for parameters
    string inpf(argv[1]);

    // initialize class 
    model reader( inpf );

    // open new xtc file for output writing
    XDRFILE *xtcout = xdrfile_open( reader.outf.c_str(), "w" );
    
    // reset file pointer to frame 0
    reader.find_frame( 0 );

    // loop over entire trajectory
    for ( currentSample = 0; currentSample < reader.nframes; currentSample ++ ){
        if ( currentSample % 1000 == 0 ){
            cout << "\rNow Processing frame: " << currentSample << "/" << reader.nframes;
        }
        fflush(stdout);
        reader.read_next_frame(); // read current frame
        
        int natoms = reader.nmol*4;
        reader.addMsite(); // add msite

        // write to new xtc file
        write_xtc( xtcout, natoms , reader.step, reader.gmxtime, reader.box, \
                   reader.xnew, reader.prec );
    }
    xdrfile_close( xtcout );
    cout << endl << "DONE!" << endl;
}
