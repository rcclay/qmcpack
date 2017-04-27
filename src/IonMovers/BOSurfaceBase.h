#ifndef QMCPLUSPLUS_BOUPDATE_BASE_H
#define QMCPLUSPLUS_BOUPDATE_BASE_H


namespace qmcplusplus
{

class BOSurfaceBase 
{
  typedef std::vector<std::vector<double> > ParticlePos_t;
  public:
    BOSurfaceBase(){std::cout<<"Calling the BOSurfaceBase constructor\n";};  //constructor
    ~BOSurfaceBase(){};  //destructor

//    BOSurfaceBase(const BOSurfaceBase& u){}; //copy constructor
//    BOSurfaceBase& operator=(const BOSurfaceBase& u){}; //copy assignment
    
//    virtual bool put(xml){return 1;}
//    virtual bool compute(){return 1;};
//    virtual bool update(){return 1;};

    virtual bool E(const ParticlePos_t& R1, double& dE, double& err){return 0;};
    virtual bool dE(const ParticlePos_t& R0, const ParticlePos_t& R1, double& dE, double& err){return 1;};
    virtual bool dE(const ParticlePos_t& R0, const ParticlePos_t& R1, int iat, double& dE, double& err){return 1;};
    virtual bool F(const ParticlePos_t& R1, ParticlePos_t& F){return 1;};
  private:
    //nothing for now
};

}

#endif
