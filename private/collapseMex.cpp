#include <mex.h>
#include <cmath>
#include <vector>
#ifdef USEOMP
#include <omp.h>
#endif

#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))
// I like degrees :)
#define sind(x) (sin((x) * M_PI / 180))
#define cosd(x) (cos((x) * M_PI / 180))

// return I[x,y] via bilinear interpolation
inline float interp( float *I, int h, int w, float x, float y ) {
  x = x<0 ? 0 : (x>w-1.001 ? w-1.001 : x);
  y = y<0 ? 0 : (y>h-1.001 ? h-1.001 : y);
  int x0=int(x), y0=int(y), x1=x0+1, y1=y0+1;
  float dx0=x-x0, dy0=y-y0, dx1=1-dx0, dy1=1-dy0;
  return I[x0*h+y0]*dx1*dy1 + I[x1*h+y0]*dx0*dy1 +
    I[x0*h+y1]*dx1*dy0 + I[x1*h+y1]*dx0*dy0;
}

// simple implementation of matlab's linspace
std::vector<float> linspace(float a, float b, int n) {
  std::vector<float> v(n);
  float step = (b-a)/(n-1);
  v[0] = a;
  for (int i = 1; i < n; i++)
    v[i] = v[i-1] + step;
  return v;
}

// V = mexFunction(W,thetas,gtWidth,nDists,nOrients,shrink,nThreads)
void mexFunction( int nl, mxArray *pl[], int nr, const mxArray *pr[] )
{
  if (nr != 7)
    mexErrMsgTxt("Wrong number of inputs");
  if (mxGetClassID(pr[0]) != mxSINGLE_CLASS ||
      mxGetClassID(pr[1]) != mxSINGLE_CLASS)
    mexErrMsgTxt("W,theta should be of type single");

  // get inputs
  float* W         = (float*) mxGetData(pr[0]);
  float* thetas    = (float*) mxGetData(pr[1]);
  int    gtWidth   = (int) mxGetScalar(pr[2]);
  int    nDists    = (int) mxGetScalar(pr[3]);
  int    nOrients  = (int) mxGetScalar(pr[4]);
  int    shrink    = (int) mxGetScalar(pr[5]);
  int    nThreads  = (int) mxGetScalar(pr[6]);

  // get dimensions
  const mwSize* wSize = mxGetDimensions(pr[0]);
  const int h = (int) wSize[0];
  const int w = (int) wSize[1];
  const int vDims[3] = {h, w, nOrients};

  // create output
  pl[0] = mxCreateNumericArray(3, vDims, mxSINGLE_CLASS, mxREAL);
  float* V = (float*) mxGetData(pl[0]);

  // make dists array like clustMasks.m
  int r = gtWidth/2-1;
  int centerIndex = (nDists-1)/2;
  std::vector<float> dists = linspace(-r, r, nDists);
  dists[centerIndex] = 0; // want exactly 0, not approx

  #ifdef USEOMP
  nThreads = min(nThreads,omp_get_max_threads());
  #pragma omp parallel for num_threads(nThreads)
  #endif
  for (int x = 0; x < w; x++) {
    for (int o = 0; o < nOrients; o++) {
      // make step vector in direction theta
      float theta = thetas[o];
      float dx = cosd(theta), dy = sind(theta);
      float m = max(std::abs(dx),std::abs(dy));
      dx /= m; dy /= m;
      for (int i = 0; i < nDists; i++) {
        float* W1 = W + (o*nDists + i)*h*w;
        float d = dists[i]/shrink;
        for (int y = 0; y < h; y++) {
          V[o*h*w + x*h + y] += d == 0 ?  W1[x*h + y] :
            interp(W1, h, w, x+d*dx, y-d*dy);
        }
      }
    }
  }
}
