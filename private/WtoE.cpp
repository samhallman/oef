#include <mex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef USEOMP
#include <omp.h>
#endif

typedef unsigned int uint32;
typedef unsigned short uint16;
typedef unsigned char uint8;
template<typename T> inline T min( T x, T y ) { return x<y ? x : y; }

// construct lookup array for mapping fids to channel indices
uint32* buildLookup( int *dims, int w ) {
  int c, r, z, n=w*w*dims[2]; uint32 *cids=new uint32[n]; n=0;
  for(z=0; z<dims[2]; z++) for(c=0; c<w; c++) for(r=0; r<w; r++)
    cids[n++] = z*dims[0]*dims[1] + c*dims[0] + r;
  return cids;
}

// Eorient = mexFunction(model,W,I)
void mexFunction( int nl, mxArray *pl[], int nr, const mxArray *pr[] )
{
  // get inputs
  mxArray *model = (mxArray*) pr[0];
  float *W = (float*) mxGetData(pr[1]);
  float *I = (float*) mxGetData(pr[2]);

  // extract relevant fields from model and options
  float *thrs = (float*) mxGetData(mxGetField(model,0,"thrs"));
  uint8 *segs = (uint8*) mxGetData(mxGetField(model,0,"segs"));
  uint16 *eBins = (uint16*) mxGetData(mxGetField(model,0,"eBins"));
  uint32 *eBnds = (uint32*) mxGetData(mxGetField(model,0,"eBnds"));
  mxArray *opts = mxGetField(model,0,"opts");
  const int imWidth = (int) mxGetScalar(mxGetField(opts,0,"imWidth"));
  const int gtWidth = (int) mxGetScalar(mxGetField(opts,0,"gtWidth"));
  const int nOrients = (int) mxGetScalar(mxGetField(opts,0,"nOrients"));
  const int stride = (int) mxGetScalar(mxGetField(opts,0,"stride"));
  int sharpen = (int) mxGetScalar(mxGetField(opts,0,"sharpen"));
  int nThreads = (int) mxGetScalar(mxGetField(opts,0,"nThreads"));
  const mwSize *segsSize = mxGetDimensions(mxGetField(model,0,"segs"));
  const int nDists = (int) mxGetScalar(mxGetField(opts,0,"nDists"));
  const int nClusts = (int) segsSize[2]; const int nBnds =
    int(mxGetNumberOfElements(mxGetField(model,0,"eBnds"))-1) / nClusts;
  const char *msgSharpen="Model supports sharpening of at most %i pixels!\n";
  if( sharpen>nBnds-1 ) { sharpen=nBnds-1; mexPrintf(msgSharpen,sharpen); }

  // get dimensions and constants
  const mwSize *imgSize = mxGetDimensions(pr[2]);
  const int h = (int) imgSize[0];
  const int w = (int) imgSize[1];
  const int Z = mxGetNumberOfDimensions(pr[2])<=2 ? 1 : imgSize[2];
  const mwSize *fidsSize = mxGetDimensions(mxGetField(model,0,"fids"));
  const int nTreeNodes = (int) fidsSize[0];
  const int nTrees = (int) fidsSize[1];
  const int h1 = (int) ceil(double(h-imWidth)/stride);
  const int w1 = (int) ceil(double(w-imWidth)/stride);
  const int h2 = h1*stride+gtWidth;
  const int w2 = w1*stride+gtWidth;
  const int imgDims[3] = {h,w,Z};
  const int outDims[3] = {h2,w2,nOrients};
  const int segDims[5] = {gtWidth,gtWidth,h1,w1,nClusts};

  // construct lookup tables
  uint32 *iids = buildLookup( (int*)imgDims, gtWidth );
  uint32 *eids = buildLookup( (int*)outDims, gtWidth );

  // create outputs
  pl[0] = mxCreateNumericArray(3,outDims,mxSINGLE_CLASS,mxREAL);
  float *E = (float*) mxGetData(pl[0]);
  if(nl>1) pl[1] = mxCreateNumericArray(5,segDims,mxUINT8_CLASS,mxREAL);
  uint8 *segsOut; if(nl>1) segsOut = (uint8*) mxGetData(pl[1]);

  // compute edge maps (avoiding collisions from parallel executions)
  if( !sharpen ) for( int c0=0; c0<gtWidth/stride; c0++ ) {
    #ifdef USEOMP
    #pragma omp parallel for num_threads(nThreads)
    #endif
    for( int c=c0; c<w1; c+=gtWidth/stride ) {
      for( int r=0; r<h1; r++ ) {
        for( int k=0; k<nClusts; k++ ) {
          float wk = W[ r + c*h1 + k*h1*w1 ];
          if( wk==0 ) continue; int orient = k/nDists;
          float *E1 = E + (r*stride) + (c*stride)*h2 + orient*w2*h2;
          int b0=eBnds[k*nBnds], b1=eBnds[k*nBnds+1];
          if( b0==b1 ) mexErrMsgTxt("bug: b0==b1");
          for( int b=b0; b<b1; b++ ) E1[eids[eBins[b]]] += wk;
          //if(nl>1) memcpy(segsOut+(r+c*h1+k*h1*w1)*gtWidth*gtWidth,
            //segs+k*gtWidth*gtWidth,gtWidth*gtWidth*sizeof(uint8));
        }
      }
    }
  }

  // computed sharpened edge maps, snapping to local color values
  if( sharpen ) {
    // compute neighbors array
    const int g=gtWidth; uint16 N[4096*4];
    for( int c=0; c<g; c++ ) for( int r=0; r<g; r++ ) {
      int i=c*g+r; uint16 *N1=N+i*4;
      N1[0] = c>0 ? i-g : i; N1[1] = c<g-1 ? i+g : i;
      N1[2] = r>0 ? i-1 : i; N1[3] = r<g-1 ? i+1 : i;
    }
    #ifdef USEOMP
    #pragma omp parallel for num_threads(nThreads)
    #endif
    for( int c=0; c<w1; c++ ) for( int r=0; r<h1; r++ ) {
      for( int k=0; k<nClusts; k++ ) {
        // get current segment and copy into S
        float wk = W[ r + c*h1 + k*h1*w1 ]; if( wk == 0 ) continue;
        uint8 S0[4096], *S=(nl<=1) ? S0 : segsOut+(r+c*h1+k*h1*w1)*g*g;
        memcpy(S,segs+k*g*g, g*g*sizeof(uint8));
        // compute color model for each segment using every other pixel
        int ci, ri, s, z; float ns[100], mus[1000];
        const float *I1 = I+(c*stride+(imWidth-g)/2)*h+r*stride+(imWidth-g)/2;
        for( s=0; s<2; s++ ) { ns[s]=0; for( z=0; z<Z; z++ ) mus[s*Z+z]=0; }
        for( ci=0; ci<g; ci+=2 ) for( ri=0; ri<g; ri+=2 ) {
          s = S[ci*g+ri]; ns[s]++;
          for( z=0; z<Z; z++ ) mus[s*Z+z]+=I1[z*h*w+ci*h+ri];
        }
        for(s=0; s<2; s++) for( z=0; z<Z; z++ ) mus[s*Z+z]/=ns[s];
        // update segment S according to local color values
        int b0=eBnds[k*nBnds], b1=eBnds[k*nBnds+sharpen];
        for( int b=b0; b<b1; b++ ) {
          float vs[10], d, e0=0, e1=0;
          for( z=0; z<Z; z++ ) vs[z]=I1[iids[eBins[b]]+z*h*w];
          for( z=0; z<Z; z++ ) {
            d=mus[0*Z+z]-vs[z]; e0 += d*d;
            d=mus[1*Z+z]-vs[z]; e1 += d*d;
          }
          S[eBins[b]] = e1 < e0;
        }
        // convert mask to edge maps (examining expanded set of pixels)
        const int orient = k/nDists, offset = orient*w2*h2;
        float *E1 = E + c*stride*h2 + r*stride; b1=eBnds[k*nBnds+sharpen+1];
        for( int b=b0; b<b1; b++ ) {
          int i=eBins[b]; uint8 s=S[i]; uint16 *N1=N+i*4;
          if( s!=S[N1[0]] || s!=S[N1[1]] || s!=S[N1[2]] || s!=S[N1[3]] )
            E1[ offset + eids[i] ] += wk;
        }
      }
    }
  }

  // free memory
  delete [] iids; delete [] eids;
}
