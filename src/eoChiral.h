class EoS;
class EoSaux;
class EoSChiral : public EoS {
private:
 EoSaux *eossmall, *eosbig;

public:
 EoSChiral(void);
 ~EoSChiral(void);

 virtual void eos(double e, double nb, double nq, double ns, double &_T,
                  double &_mub, double &_muq, double &_mus, double &_p);
 virtual double p(double e, double nb, double nq, double ns);

 double dpe(double e, double nb, double nq, double ns);
 virtual double cs2(double e) { return dpe(e,0,0,0); }
};
