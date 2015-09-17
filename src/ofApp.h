#pragma once

#include <stack>
#include <vector>
#include "ofMain.h"
#include "ofMath.h"
#include "ofxGui.h"

#define H_MARGIN 100
#define W_MARGIN 100
#define TH_MARGIN 20
#define TW_MARGIN 20

#define FRAME_RATE 60 
#define MAX_N 1000

using namespace std;

// -----------------------------------------------------------------------------------------
// Pseudo-random number generator
// -----------------------------------------------------------------------------------------
static double ran2(long *idum) 
{
  // Long period (> 2e18) random number generator of L'Ecuyer 
  // with Bays-Durham shuffle and added safeguards. 
  // Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values). 
  // Call with idum a negative integer to initialize; 
  // thereafter, do not alter idum between successive deviates in a sequence. 
  // RNMX should approximate the largest floating value that is less than 1. 
  //
  // change into a double version

  static const long IM1 = 2147483563;    
  static const long IM2 = 2147483399;    
  static const double AM = 1.0L / IM1;    
  static const long IMM1 = (IM1-1);      
  static const long IA1 = 40014;          
  static const long IA2 = 40692;         
  static const long IQ1 = 53668;        
  static const long IQ2 = 52774;        
  static const long IR1 = 12211;         
  static const long IR2 = 3791;         
  static const long NTAB = 32;          
  static const long NDIV =(1+IMM1/NTAB); 
  //static const double EPS = 4.4e-16; //1.2e-7;    
  //static const double RNMX = (1.0-EPS); 

  int j; 
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  if (*idum <= 0) {              // Initialize. 
    if (-(*idum) < 1) *idum=1;   // Be sure to prevent idum = 0. 
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {    // Load the shuffle table (after 8 warm-ups).
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1; 
      if (*idum < 0) *idum += IM1; 
      if (j < NTAB) iv[j] = *idum;
    } 
    iy=iv[0]; 
  } 
  k=(*idum)/IQ1;                 // Start here when not initializing.
  *idum=IA1*(*idum-k*IQ1)-k*IR1; // Compute idum=(IA1*idum) % IM1 without overflows by Schrage's 
                                 // method. 
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2; 
  idum2=IA2*(idum2-k*IQ2)-k*IR2; // Compute idum2=(IA2*idum) % IM2 likewise.
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;                     // Will be in the range 0..NTAB-1.
  iy=iv[j]-idum2;        // Here idum is shuffled, idum and idum2 are combined to generate output. 
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;

  double temp=AM*iy;
  if (temp >= 1) cout << "WARNING: random 1" << endl;
  //if ( temp > RNMX) return RNMX; // Because users don't expect endpoint values.
  //else return temp; 
  return temp;
} 

// -----------------------------------------------------------------------------------------
// inverse normal
//   see Peter J. Acklam  http://home.online.no/~pjacklam/notes/invnorm/
//   this is supposed to give full machine precision
// -----------------------------------------------------------------------------------------
static double inverse_normal(double p)
{
  static const double a1 = -3.969683028665376e+01;
  static const double a2 =  2.209460984245205e+02;
  static const double a3 = -2.759285104469687e+02;
  static const double a4 =  1.383577518672690e+02;
  static const double a5 = -3.066479806614716e+01;
  static const double a6 =  2.506628277459239e+00;
 
  static const double b1 = -5.447609879822406e+01;
  static const double b2 =  1.615858368580409e+02;
  static const double b3 = -1.556989798598866e+02;
  static const double b4 =  6.680131188771972e+01;
  static const double b5 = -1.328068155288572e+01;
  
  static const double c1 = -7.784894002430293e-03;
  static const double c2 = -3.223964580411365e-01;
  static const double c3 = -2.400758277161838e+00;
  static const double c4 = -2.549732539343734e+00;
  static const double c5 =  4.374664141464968e+00;
  static const double c6 =  2.938163982698783e+00;

  static const double d1 =  7.784695709041462e-03;
  static const double d2 =  3.224671290700398e-01;
  static const double d3 =  2.445134137142996e+00;
  static const double d4 =  3.754408661907416e+00;
  
  static const double p_low = 0.02425;
  static const double p_high = 1 - p_low;

  double q, r, x, e, u;

  if (p <= 0) {
    cout << "WARNING: inverse normal <=0 " << endl;
    return -HUGE_VAL; /* - infinity */
  }
  if (p >= 1) {
    cout << "WARNING: inverse normal >=1 " << endl;
    return HUGE_VAL;  /*   infinity */
  }
  
  /* Rational approximation for lower region */
  if (0 < p && p < p_low) {
    q = sqrt(-2*log(p));
    x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1);
  }

  /* Rational approximation for central region */
  else if (p_low <= p && p <= p_high) {
    q = p - 0.5;
    r = q*q;
    x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1);
  }

  /* Rational approximation for upper region */
  else {
    q = sqrt(-2*log(1-p));
    x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1);
  }

  /* One iteration of Halley's rational method */
  e = 0.5 * erfc(-x/sqrt(2.0)) - p;
  u = e * sqrt(2*M_PI) * exp(x*x/2);
  x = x - u/(1 + x*u/2);

  return x;
}

class ofApp : public ofBaseApp{

	public:
		void setup();
		void update();
		void draw();

		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);
   

        int inv_cum_dist(double x, double* dist, int len) 
        {
            int n = 1;
            while (x < dist[n-1] && n < len) n++;
            return n;
        }
        void rand_direction_2d(double x, ofPoint& pos)
        {
            if (x < 0.25)
                pos[0] += 1.0;
            else if (x < 0.5)
                pos[0] -= 1.0;
            else if (x < 0.75)
                pos[1] += 1.0;
            else
                pos[1] -= 1.0;
        }
        void rand_direction(double r1, double r2, ofPoint& pos)
        {
            pos[0] += inverse_normal(r1);
            pos[1] += inverse_normal(r2);
        }

        // State variables
        long n;
        double t;
        ofPoint sib_pos;
        ofPoint exp_pos;
        vector<ofPoint> sib_pos_hist;
        vector<ofPoint> exp_pos_hist;
        long sib_arrival;
        double sib_arrival_t;
        double sib_cd_dt;
        long sib_wait;
        long exp_arrival;
        double exp_arrival_t;
        double exp_cd_dt;
        long exp_wait;

        //ofMesh sib_graph;
        //ofMesh exp_graph;
        ofColor sib_color;
        ofColor exp_color;

        // Parameters
        double D_alpha;
        double alpha;
        double omega;
        double r;
        double dX;
        double dT;
        double base_dT;

        long max_n;
        
        double* survival_sib;
        double* survival_exp;

        double* sib_stats;
        double* exp_stats;
        double* sib_dist;
        long sjc;

        long seed;

        // Mouse drag parameters
        int mouse_x_start;

        
    	ofxPanel gui;
    	ofxIntSlider rate;
};
