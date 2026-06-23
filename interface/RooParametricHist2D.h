#ifndef ROOPARAMETRICHIST2D
#define ROOPARAMETRICHIST2D

#include <math.h>
#include <vector>

#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooAbsPdf.h"
#include "RooAddition.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
//#include "RooAbsData.h"
#include "TH1.h"
#include "TH2F.h"
#include "TH2D.h"
  
class RooParametricHist2D : public RooAbsPdf {
public:
  
  RooParametricHist2D() {} ;
  RooParametricHist2D (const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _y, RooArgList& _pars, const TH2 &_shape);
  
  RooParametricHist2D(const RooParametricHist2D& other, const char* name=0);
  TObject* clone(const char* newname) const override { return new RooParametricHist2D(*this,newname); }
  inline ~RooParametricHist2D () override{};
  Int_t getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars, const char* rangeName=0) const override ;
  Double_t analyticalIntegral(Int_t, const char* rangeName=0) const override ;   

  void addMorphs(RooDataHist&, RooDataHist&, RooRealVar&, double );
  
  //RooAddition & getYieldVar(){return sum;};
  
protected:
  
  RooRealProxy x;
  RooRealProxy y;
  //RooAddition sum;
  RooListProxy pars;
  RooListProxy _coeffList;
  mutable int N_bins; 
  mutable int N_bins_x;
  mutable int N_bins_y;
  mutable std::vector<double> bins_x; 
  mutable std::vector<double> bins_y;
  mutable std::vector<double> widths_x;
  mutable std::vector<double> widths_y;

  mutable double _smoothRegion;
  mutable bool   _hasMorphs;
  mutable std::vector<std::vector <double> > _diffs;
  mutable std::vector<std::vector <double> > _sums;
  double evaluateMorphFunction(int) const;

  void initializeBins(const TH2&) const;
  //void initializeNorm();
  
  double evaluatePartial() const ;
  double evaluateFull() const ;
  Double_t evaluate() const override ;
  double getFullSum() const ;
 
  mutable double cval;
  void update_cval(double r){cval=r;};

  inline double smoothStepFunc(double x) const {
    if (fabs(x) >= _smoothRegion) return x > 0. ? +1. : -1.;
    double xnorm = x/_smoothRegion, xnorm2 = xnorm*xnorm;
    return 0.125 * xnorm * (xnorm2 * (3.*xnorm2 - 10.) + 15.);
  }
private:
   ClassDefOverride(RooParametricHist2D, 2)
};

#endif
