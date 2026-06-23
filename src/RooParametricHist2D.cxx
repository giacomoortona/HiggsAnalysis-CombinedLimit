/***************************************************************************** 
 *****************************************************************************/ 


#include "Riostream.h" 

#include "RooAbsData.h"
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "../interface/RooParametricHist2D.h" 

#include <math.h> 
#include "TMath.h" 
#include "RooFormulaVar.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooFit.h"

#include "TFile.h"
#include <typeinfo>

//using namespace RooFit ;

ClassImp(RooParametricHist2D) 

RooParametricHist2D::RooParametricHist2D( const char *name, 
                                          const char *title, 
                                          RooAbsReal& _x,
                                          RooAbsReal& _y,
                                          RooArgList& _pars,    // stores all bins - Ex. a 3x3 space [[1,2,3],[4,5,6],[7,8,9]] would be stored [1,2,3,4,5,6,7,8,9]
                                          const TH2 &_shape  ) :// need this to sort bins into x and y
    RooAbsPdf(name,title),
    x("observable_x","observable_x",this,_x),
    y("observable_y","observable_y",this,_y),
    pars("pars","pars",this),
    _coeffList("_coeffList","_coeffList",this)
    //SM_shape("SM_shape","SM_shape",this,_SM_shape),
{ 
    pars.add(_pars);

    if ( pars.getSize() != _shape.GetNbinsY()*_shape.GetNbinsX() ){
        std::cout << " Warning, number of parameters not equal to number of bins in shape histogram! " << std::endl;
    }
    
    initializeBins(_shape);
    
//  initializeNorm();
    cval = -1; 
    _smoothRegion = 0.;
    _hasMorphs = false;
} 

//_____________________________________________________________________________
RooParametricHist2D::RooParametricHist2D( const RooParametricHist2D& other, 
                                          const char* name) :
    RooAbsPdf(other, name),
    x("observable_x",this,other.x),
    y("observable_y",this,other.y),
    pars("_pars",this,RooListProxy()),
    _coeffList("_coeffList",this,RooListProxy())
{
    //std::cout << "Cloning RooParametricHist2D with name " <<other.GetName() << std::endl;
    N_bins = other.N_bins;
    N_bins_x = other.N_bins_x;
    N_bins_y = other.N_bins_y;
    _smoothRegion = other._smoothRegion;
    _hasMorphs = other._hasMorphs;
    //sum    = other.sum;

    pars.add(other.pars);
    _coeffList.add(other._coeffList);

    for(int i=0; i<=N_bins_x; ++i) {
        bins_x.push_back(other.bins_x[i]);
        if (i<N_bins_x) widths_x.push_back(other.widths_x[i]);
    }

    for(int i=0; i<=N_bins_y; ++i) {
        bins_y.push_back(other.bins_y[i]);
        if (i<N_bins_y) widths_y.push_back(other.widths_y[i]);
    }

    if (other._hasMorphs) {
        for (int i=0; i<N_bins; ++i) {
            std::vector<double> su;
            std::vector<double> di;
            for (int j=0; j<other._coeffList.getSize(); ++j) {
                su.push_back(other._sums[i][j]);
                di.push_back(other._diffs[i][j]);
            }
            _sums.push_back(su);
            _diffs.push_back(di);
        }
    }

    cval = other.cval; 

}
void RooParametricHist2D::initializeBins(const TH2 &shape) const {
    // std::cout << "Initializing RooParametricHist2D bins" << std::endl;
    ///std::vector<double> bins;
    //std::vector<double> widths;
    N_bins_x = shape.GetNbinsX();
    N_bins_y = shape.GetNbinsY();
    N_bins = N_bins_x*N_bins_y;


    for(int i=1; i<=N_bins_x+1; ++i) {
        bins_x.push_back(shape.GetXaxis()->GetBinLowEdge(i));
        if (i<=N_bins_x) widths_x.push_back(shape.GetXaxis()->GetBinWidth(i));
    }

    for(int j=1; j<=N_bins_y+1; ++j) {
        bins_y.push_back(shape.GetYaxis()->GetBinLowEdge(j));
        if (j<=N_bins_y) widths_y.push_back(shape.GetYaxis()->GetBinWidth(j));
    }
    // std::cout << "Initialize bins_x = ";
    // for (auto i = bins_x.begin(); i != bins_x.end(); ++i){
    //     std::cout << *i << ' ';
    // }
    // std::cout << std::endl;
    

    // std::cout << "Initialize bins_y = ";
    // for (auto i = bins_y.begin(); i != bins_y.end(); ++i){
    //     std::cout << *i << ' ';
    // }
    // std::cout << std::endl;

}

double RooParametricHist2D::getFullSum() const {
    double sum=0;

    for (int i = 0; i < pars.getSize(); ++i) {
        double thisVal = static_cast<RooAbsReal&>(pars[i]).getVal();
        if (_hasMorphs) thisVal*=evaluateMorphFunction(i);
        sum += thisVal;
    }

    return sum;
}

Int_t RooParametricHist2D::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet & analVars, const char*) const  {
    if (matchArgs(allVars,analVars,x,y)){
        return 1;
    }
    return 0;
}

Double_t RooParametricHist2D::analyticalIntegral(Int_t code, const char* rangeName) const    
{
    assert(code==1) ;

    // Case without range is trivial: p.d.f is by construction normalized 
    if (!rangeName) {
        //return 1;//getFullSum() ;
        return getFullSum();
    } else {
        std::cout << "Analytical integral for range " << rangeName << " in RooParametricHist2D not yet implemented" << std::endl;
    }
    
    // // Case with ranges, calculate integral explicitly 
    // double xmin = x.min(rangeName) ;
    // double xmax = x.max(rangeName) ;
    // double ymin = y.min(rangeName) ;
    // double ymax = y.max(rangeName) ;

    // double sum=0 ;
    // int i ;
    // for (i=1 ; i<=N_bins_x ; i++) {
    //     double binVal = (static_cast<RooAbsReal*>(pars.at(i-1))->getVal())/widths[i-1]; 
    //     if (bins[i-1]>=xmin && bins[i]<=xmax) {
    //         // Bin fully in the integration domain
    //         sum += (bins[i]-bins[i-1])*binVal ;
    //     } else if (bins[i-1]<xmin && bins[i]>xmax) {
    //         // Domain is fully contained in this bin
    //         sum += (xmax-xmin)*binVal ;
    //         // Exit here, this is the last bin to be processed by construction
    //         return sum/getFullSum() ;
    //     } else if (bins[i-1]<xmin && bins[i]<=xmax && bins[i]>xmin) {
    //         // Lower domain boundary is in bin
    //         sum +=  (bins[i]-xmin)*binVal ;
    //     } else if (bins[i-1]>=xmin && bins[i]>xmax && bins[i-1]<xmax) {
    //         sum +=  (xmax-bins[i-1])*binVal ;
    //         // Upper domain boundary is in bin
    //         // Exit here, this is the last bin to be processed by construction
    //         return sum ;
    //     }
    // }
    return 0;
}

void RooParametricHist2D::addMorphs(RooDataHist &hpdfU, RooDataHist &hpdfD, RooRealVar &cVar, double smoothRegion){

    if (!_hasMorphs) {
        for (int i=0; i<N_bins; ++i) {
            std::vector<double> difs;
            std::vector<double> sums;
            _diffs.push_back(difs);
            _sums.push_back(sums);
        }
    }
    for (int i=0; i<N_bins; ++i) {
        double f0 = static_cast<RooAbsReal*>(pars.at(i))->getVal();

        hpdfU.get(i); hpdfD.get(i);
        double dh = (hpdfU.weight()-f0);
        double dl = (hpdfD.weight()-f0);
        _diffs[i].push_back(dh-dl);
        _sums[i].push_back(dh+dl);
    }
    _coeffList.add(cVar);
    _hasMorphs = true;
    _smoothRegion = smoothRegion;
}

double RooParametricHist2D::evaluateMorphFunction(int j) const
{
    double scale=1.0;
    if (!_hasMorphs) return scale;

    int ndim = _coeffList.getSize();
    double f0 = static_cast<RooAbsReal*>(pars.at(j))->getVal();
    for (int i = 0; i < ndim; ++i) {
        double x = (dynamic_cast<RooRealVar*>(_coeffList.at(i)))->getVal();
        double a = 0.5*x, b = smoothStepFunc(x);
        scale *= 1+(1./f0) * a*(_diffs[j][i] + b*_sums[j][i]);
    }
    return scale;
}

double RooParametricHist2D::evaluatePartial() const
{
  auto itx = std::upper_bound(std::begin(bins_x), std::end(bins_x), x);
  if ( itx == std::begin(bins_x) ) {
    // underflow
    return 0;
  }
  else if ( itx == std::end(bins_x) ) {
    // overflow
    return 0;
  }
  auto ity = std::upper_bound(std::begin(bins_y), std::end(bins_y), y);
  if ( ity == std::begin(bins_y) ) {
    // underflow
    return 0;
  }
  else if ( ity == std::end(bins_y) ) {
    // overflow
    return 0;
  }
  size_t bin_ix = std::distance(std::begin(bins_x), itx) - 1;
  size_t bin_iy = std::distance(std::begin(bins_y), ity) - 1;

  int globalbin = N_bins_x * bin_iy + bin_ix;

  RooAbsReal *retVar = (RooAbsReal*)pars.at(globalbin);

  double ret = retVar->getVal() / (widths_x[bin_ix]*widths_y[bin_iy]);

  return ret;
}

double RooParametricHist2D::evaluateFull() const
{
  auto itx = std::upper_bound(std::begin(bins_x), std::end(bins_x), x);
  if ( itx == std::begin(bins_x) ) {
    // underflow
    return 0;
  }
  else if ( itx == std::end(bins_x) ) {
    // overflow
    return 0;
  }
  auto ity = std::upper_bound(std::begin(bins_y), std::end(bins_y), y);
  if ( ity == std::begin(bins_y) ) {
    // underflow
    return 0;
  }
  else if ( ity == std::end(bins_y) ) {
    // overflow
    return 0;
  }
  size_t bin_ix = std::distance(std::begin(bins_x), itx) - 1;
  size_t bin_iy = std::distance(std::begin(bins_y), ity) - 1;

  int globalbin = N_bins_x * bin_iy + bin_ix;

  double mVar = evaluateMorphFunction(globalbin);
  RooAbsReal *retVar = (RooAbsReal*)pars.at(globalbin);

  double ret = retVar->getVal()*mVar;
  ret /= widths_x[bin_ix]*widths_y[bin_iy];
  return ret;
}

Double_t RooParametricHist2D::evaluate() const
{
  double ret = _hasMorphs ? evaluateFull() : evaluatePartial() ;
  cval=ret;
  return ret > 0 ? ret : 0;
}
