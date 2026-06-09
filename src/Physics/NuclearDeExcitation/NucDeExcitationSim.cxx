//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ March 05, 2008 - CA
   This event generation module was added in version 2.3.1. The initial
   implementation handles 16O only.
 @ Sep 15, 2009 - CA
   IsNucleus() is no longer available in GHepParticle. Use pdg::IsIon().
*/
//____________________________________________________________________________

#include <cstdlib>
#include <sstream>

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Physics/NuclearDeExcitation/NucDeExcitationSim.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"

using std::ostringstream;

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;
using namespace genie::controls;

//___________________________________________________________________________
NucDeExcitationSim::NucDeExcitationSim() :
EventRecordVisitorI("genie::NucDeExcitationSim")
{

}
//___________________________________________________________________________
NucDeExcitationSim::NucDeExcitationSim(string config) :
EventRecordVisitorI("genie::NucDeExcitationSim", config)
{

}
//___________________________________________________________________________
NucDeExcitationSim::~NucDeExcitationSim()
{

}
//___________________________________________________________________________
void NucDeExcitationSim::ProcessEventRecord(GHepRecord * evrec) const
{
  LOG("NucDeEx", pNOTICE) 
     << "Simulating nuclear de-excitation gamma rays";

  GHepParticle * nucltgt = evrec->TargetNucleus();
  if (!nucltgt) {
    LOG("NucDeEx", pINFO) 
      << "No nuclear target found - Won't simulate nuclear de-excitation";
    return;
  }

  if(nucltgt->Z()==8) this->OxygenTargetSim(evrec);

  LOG("NucDeEx", pINFO) 
     << "Done with this event";
}
//___________________________________________________________________________
void NucDeExcitationSim::OxygenTargetSim(GHepRecord * evrec) const
{
  LOG("NucDeEx", pNOTICE) 
     << "Simulating nuclear de-excitation gamma rays for Oxygen target";

  //LOG("NucDeEx", pNOTICE) << *evrec;

  GHepParticle * hitnuc = evrec->HitNucleon();
  if(!hitnuc) return;

  bool p_hole = (hitnuc->Pdg() == kPdgProton);
  double dt   = -1;

  RandomGen * rnd = RandomGen::Instance();

  //
  // ****** P-Hole 
  //
  if (p_hole) {
    // 
    // * Define all the data required for simulating deexcitations of p-hole states
    //

    // > probabilities for creating a p-hole in the P1/2, P3/2, S1/2 shells
    double Pp12 = 0.158;             // P1/2
    double Pp32 = 0.327;             // P3/2 
    double Ps12 = 0.200;             // S1/2 

    // > excited state energy levels & probabilities for P3/2-shell p-holes
    const int np32 = 3;
    double p32Elv[np32] = { 0.00632, 0.00993, 0.01070 };
    double p32Plv[np32] = { 0.899,   0.049,   0.052   }; 
    // - probabilities for deexcitation modes of P3/2-shell p-hole state '1' 
    double p32Plv1_1gamma  = 0.78;  // prob to decay via 1 gamma
    double p32Plv1_cascade = 0.22;  // prob to decay via gamma cascade

    // > excited state energy levels & probabilities for S1/2-shell p-holes
    // > many states measured by Kobayashi 2006 have no "detectable" (> 3 MeV) gamma-rays
    // > and therefore the probabilities will not sum to 1
    const int ns12 = 17;
    double s12Elv[ns12] = {   // J^pi
               0.00309,       // 1/2+
               0.00368,       // 3/2+
               0.00385,       // 5/2+
               0.00444,       // 2+
               0.00492,       // 0-
               0.00511,       // 2-
               0.00569,       // 1-
               0.00583,       // 3-
               0.00609,       // 1-
               0.00620,       // 1+
               0.00645,       // 3+
               0.00659,       // 0+
               0.00673,       // 3-
               0.00690,       // 0-
               0.00701,       // 2+
               0.00703,       // 2+
               0.00734 };     // 2-
    double s12Plv[ns12] = { 
               0.0300,        
               0.0420,        
               0.0460,        
               0.0580,        
               0.0520,        
               0.0000,        
               0.0450,        
               0.0054,      
               0.0000,        
               0.0000,        
               0.0280,        
               0.0000,       
               0.0043,        
               0.0000,        
               0.0335,        // 7.01 / 7.03 MeV states not individually resolvable
               0.0335,        // combined probability is 6.7%; split in two
               0.0570 };
    // - gamma energies and probabilities for S1/2-shell p-hole excited 

    //   states with >1 deexcitation modes
    //   sel_state = [2,6,7,9,10,12,14,16]
    const int ns12lv2 = 3;     // 3.85 MeV triple
    double s12Elv2[ns12lv2]    = { 0.00309, 0.00368, 0.00385 };
    double s12Plv2[ns12lv2]    = { 0.012,   0.363,   0.625   };

    const int ns12lv6 = 2;     // 5.69 MeV double
    double s12Elv6[ns12lv6]    = { 0.00338, 0.00569 };
    double s12Plv6[ns12lv6]    = { 0.639,   0.361   };

    const int ns12lv7 = 2;     // 5.83 MeV double
    double s12Elv7[ns12lv7]    = { 0.00511, 0.00583 };
    double s12Plv7[ns12lv7]    = { 0.629,   0.213   };

    const int ns12lv9 = 2;     // 6.20 MeV double
    double s12Elv9[ns12lv9]    = { 0.00389, 0.00620 };
    double s12Plv9[ns12lv9]    = { 0.769,   0.231   };

    const int ns12lv10 = 2;    // 6.45 MeV double
    double s12Elv10[ns12lv10]  = { 0.00511, 0.00644 };
    double s12Plv10[ns12lv10]  = { 0.081,   0.701   };

    const int ns12lv12 = 2;    // 6.73 MeV double
    double s12Elv12[ns12lv12]  = { 0.00609, 0.00673 };
    double s12Plv12[ns12lv12]  = { 0.036,   0.964   };

    const int ns12lv14 = 2;    // 7.01 MeV double
    double s12Elv14[ns12lv14]  = { 0.00609, 0.00701 };
    double s12Plv14[ns12lv14]  = { 0.014,   0.986   };

    const int ns12lv16 = 3;    // 7.34 MeV triple
    double s12Elv16[ns12lv16]  = { 0.00609, 0.00673, 0.00734 };
    double s12Plv16[ns12lv16]  = { 0.490,   0.343,   0.167   };
     
    // Select one of the P1/2, P3/2 or S1/2
    double rshell = rnd->RndDec().Rndm();
    //
    // >> P1/2 shell
    //
    if(rshell < Pp12) {
        LOG("NucDeEx", pNOTICE) 
          << "Hit nucleon left a P1/2 shell p-hole. Remnant is at g.s.";
        return;
    } 
    //
    // >> P3/2 shell
    //
    else
    if(rshell < Pp12 + Pp32) {
        LOG("NucDeEx", pNOTICE) 
            << "Hit nucleon left a P3/2 shell p-hole";
        // Select one of the excited states 
        double rdecmode  = rnd->RndDec().Rndm();        
        double prob_sum  = 0;
        int    sel_state = -1;
        for(int istate=0; istate<np32; istate++) {
            prob_sum += p32Plv[istate];
            if(rdecmode < prob_sum) {
              sel_state = istate;
              break;
            }
        }
        LOG("NucDeEx", pNOTICE) 
            << "Selected P3/2 excited state = " << sel_state;

        // Decay that excited state
        // >> 6.32 MeV state
        if(sel_state==0) { 
            this->AddPhoton(evrec, p32Elv[0], dt);
        } 
        // >> 9.93 MeV state
        else 
        if(sel_state==1) {    
            double r = rnd->RndDec().Rndm();        
            // >>> emit a single gamma 
            if(r < p32Plv1_1gamma) {
               this->AddPhoton(evrec, p32Elv[1], dt);
            }
            // >>> emit a cascade of gammas 
            // >>> cascade goes from 9.93 --> 6.32 --> g.s.
            else 
            if(r < p32Plv1_1gamma + p32Plv1_cascade) {
               this->AddPhoton(evrec, p32Elv[0],           dt);  // 6.32 MeV
               this->AddPhoton(evrec, p32Elv[1]-p32Elv[0], dt);  // 3.61 MeV (the two total 9.93 MeV)
            }
        }
        // >> 10.7 MeV state 
        else 
        if(sel_state==2) {
           // Above the particle production threshold - need to emit
           // a 0.5 MeV kinetic energy proton.
           // Will neglect that given that it is a very low energy
           // kinetic energy nucleon and the intranuke break-up nucleon
           // cross sections are already tuned.
           return;
        }
    } //p3/2
    //
    // >> S1/2 shell
    //
    else if (rshell < Pp12 + Pp32 + Ps12) {
        LOG("NucDeEx", pNOTICE) 
            << "Hit nucleon left an S1/2 shell p-hole";
        // Select one of the excited states caused by a S1/2 shell hole
        double rdecmode  = rnd->RndDec().Rndm();  
        double prob_sum  = 0;
        int    sel_state = -1;
        for(int istate=0; istate<ns12; istate++) {
            prob_sum += s12Plv[istate];
            if(rdecmode < prob_sum) {
              sel_state = istate;
              break;
            }
        }
        LOG("NucDeEx", pNOTICE) 
            << "Selected S1/2 excited state = " << sel_state;

        if(sel_state == -1) return;  // since the gamma-ray emitting states do not exhaust the total available s-shell states

        // Decay that excited state
        bool multiple_decay_modes = 
              (sel_state==2 || sel_state==6 || sel_state==7 || sel_state==9 || sel_state==10 || sel_state==12 || sel_state==14 || sel_state==16);
        if(!multiple_decay_modes) {
          this->AddPhoton(evrec, s12Elv[sel_state], dt); 
        } else {
          int ndec = -1;
          double * pdec = 0, * edec = 0;
          switch(sel_state) {
           case(2) : 
              ndec = ns12lv2;  pdec = s12Plv2;  edec = s12Elv2;
              break;
           case(6) : 
              ndec = ns12lv6;  pdec = s12Plv6;  edec = s12Elv6;
              break;
           case(7) : 
              ndec = ns12lv7;  pdec = s12Plv7;  edec = s12Elv7;
              break;
           case(9) : 
              ndec = ns12lv9;  pdec = s12Plv9;  edec = s12Elv9;
              break;
           case(10) : 
              ndec = ns12lv10; pdec = s12Plv10; edec = s12Elv10;
              break;
           case(12) : 
              ndec = ns12lv12;  pdec = s12Plv12;  edec = s12Elv12;
              break;
           case(14) : 
              ndec = ns12lv14;  pdec = s12Plv14;  edec = s12Elv14;
              break;
           case(16) : 
              ndec = ns12lv16;  pdec = s12Plv16;  edec = s12Elv16;
              break;
           default:
             return;
          }
          double r = rnd->RndDec().Rndm();  
          double decmode_prob_sum = 0;
          int sel_decmode = -1;
          for(int idecmode=0; idecmode < ndec; idecmode++) {
             decmode_prob_sum += pdec[idecmode];
             if(r < decmode_prob_sum) {
                 sel_decmode = idecmode;
                 break;
             }
          }
          if(sel_decmode == -1) return;
          this->AddPhoton(evrec, edec[sel_decmode], dt);  
        }//mult.dec.ch 

    } // s1/2
    else {       // continuum/other states above the shell model 
      LOG("NucDeEx", pNOTICE) 
          << "Hit nucleon populates continuum/other states; no photon emitted.";
      return;
    }
  } // p-hole

  //
  // ****** n-hole
  //
  else {
    // 
    // * Define all the data required for simulating deexcitations of n-hole states
    //

    // > probabilities for creating a n-hole in the P1/2, P3/2, S1/2 shells
    double Pp12 = 0.158;  // P1/2 
    double Pp32 = 0.310;  // P3/2 
    double Ps12 = 0.200;  // S1/2 
    //>
    double p32Elv = 0.00618;
    // > excited state energy levels & probabilities for S1/2-shell n-holes
    // > apply isospin symmetry to treat the n-hole states, as no direct measurements exist
    const int ns12 = 17;
    double s12Elv[ns12] = {   // J^pi
               0.00309,       // 1/2+
               0.00368,       // 3/2+
               0.00385,       // 5/2+
               0.00444,       // 2+
               0.00492,       // 0-
               0.00511,       // 2-
               0.00569,       // 1-
               0.00583,       // 3-
               0.00609,       // 1-
               0.00620,       // 1+
               0.00645,       // 3+
               0.00659,       // 0+
               0.00673,       // 3-
               0.00690,       // 0-
               0.00701,       // 2+
               0.00703,       // 2+
               0.00734 };     // 2-
    double s12Plv[ns12] = { 
               0.0300,        
               0.0420,        
               0.0460,        
               0.0580,        
               0.0520,        
               0.0000,        
               0.0450,        
               0.0054,      
               0.0000,        
               0.0000,        
               0.0280,        
               0.0000,       
               0.0043,        
               0.0000,        
               0.0335,  
               0.0335,
               0.0570 };
    // - gamma energies and probabilities for S1/2-shell p-hole excited 

    //   states with >1 deexcitation modes
    //   sel_state = [2,6,7,9,10,12,14,16]
    const int ns12lv2 = 3;     // 3.85 MeV triple
    double s12Elv2[ns12lv2]    = { 0.00309, 0.00368, 0.00385 };
    double s12Plv2[ns12lv2]    = { 0.012,   0.363,   0.625   };

    const int ns12lv6 = 2;     // 5.69 MeV double
    double s12Elv6[ns12lv6]    = { 0.00338, 0.00569 };
    double s12Plv6[ns12lv6]    = { 0.639,   0.361   };

    const int ns12lv7 = 2;     // 5.83 MeV double
    double s12Elv7[ns12lv7]    = { 0.00511, 0.00583 };
    double s12Plv7[ns12lv7]    = { 0.629,   0.213   };

    const int ns12lv9 = 2;     // 6.20 MeV double
    double s12Elv9[ns12lv9]    = { 0.00389, 0.00620 };
    double s12Plv9[ns12lv9]    = { 0.769,   0.231   };

    const int ns12lv10 = 2;    // 6.45 MeV double
    double s12Elv10[ns12lv10]  = { 0.00511, 0.00644 };
    double s12Plv10[ns12lv10]  = { 0.081,   0.701   };

    const int ns12lv12 = 2;    // 6.73 MeV double
    double s12Elv12[ns12lv12]  = { 0.00609, 0.00673 };
    double s12Plv12[ns12lv12]  = { 0.036,   0.964   };

    const int ns12lv14 = 2;    // 7.01 MeV double
    double s12Elv14[ns12lv14]  = { 0.00609, 0.00701 };
    double s12Plv14[ns12lv14]  = { 0.014,   0.986   };

    const int ns12lv16 = 3;    // 7.34 MeV triple
    double s12Elv16[ns12lv16]  = { 0.00609, 0.00673, 0.00734 };
    double s12Plv16[ns12lv16]  = { 0.490,   0.343,   0.167   };

    // Select one of the P1/2, P3/2 or S1/2
    double rshell = rnd->RndDec().Rndm();
    //
    // >> P1/2 shell
    //
    if(rshell < Pp12) {
        LOG("NucDeEx", pNOTICE) 
          << "Hit nucleon left a P1/2 shell n-hole. Remnant is at g.s.";
        return;
    } 
    //
    // >> P3/2 shell
    //
    else
    if(rshell < Pp12 + Pp32) {
        LOG("NucDeEx", pNOTICE) 
            << "Hit nucleon left a P3/2 shell n-hole";
        this->AddPhoton(evrec, p32Elv, dt); 
    } 
    //
    // >> S1/2 shell
    //
    else
    if(rshell < Pp12 + Pp32 + Ps12) {
        LOG("NucDeEx", pNOTICE) 
            << "Hit nucleon left a S1/2 shell n-hole";
        // Select one of the excited states caused by a S1/2 shell hole
        double rdecmode  = rnd->RndDec().Rndm();  
        double prob_sum  = 0;
        int    sel_state = -1;
        for(int istate=0; istate<ns12; istate++) {
            prob_sum += s12Plv[istate];
            if(rdecmode < prob_sum) {
              sel_state = istate;
              break;
            }
        }
        LOG("NucDeEx", pNOTICE) 
            << "Selected S1/2 excited state = " << sel_state;

        if(sel_state == -1) return;

        // Decay that excited state
        bool multiple_decay_modes = 
              (sel_state==2 || sel_state==6 || sel_state==7 || sel_state==9 || sel_state==10 || sel_state==12 || sel_state==14 || sel_state==16);
        if(!multiple_decay_modes) {
          this->AddPhoton(evrec, s12Elv[sel_state], dt); 
        } else {
          int ndec = -1;
          double * pdec = 0, * edec = 0;
          switch(sel_state) {
           case(2) : 
              ndec = ns12lv2;  pdec = s12Plv2;  edec = s12Elv2;
              break;
           case(6) : 
              ndec = ns12lv6;  pdec = s12Plv6;  edec = s12Elv6;
              break;
           case(7) : 
              ndec = ns12lv7;  pdec = s12Plv7;  edec = s12Elv7;
              break;
           case(9) : 
              ndec = ns12lv9;  pdec = s12Plv9;  edec = s12Elv9;
              break;
           case(10) : 
              ndec = ns12lv10; pdec = s12Plv10; edec = s12Elv10;
              break;
           case(12) : 
              ndec = ns12lv12;  pdec = s12Plv12;  edec = s12Elv12;
              break;
           case(14) : 
              ndec = ns12lv14;  pdec = s12Plv14;  edec = s12Elv14;
              break;
           case(16) : 
              ndec = ns12lv16;  pdec = s12Plv16;  edec = s12Elv16;
              break;
           default:
             return;
          }
          double r = rnd->RndDec().Rndm();  
          double decmode_prob_sum = 0;
          int sel_decmode = -1;
          for(int idecmode=0; idecmode < ndec; idecmode++) {
             decmode_prob_sum += pdec[idecmode];
             if(r < decmode_prob_sum) {
                 sel_decmode = idecmode;
                 break;
             }
          }
          if(sel_decmode == -1) return;
          this->AddPhoton(evrec, edec[sel_decmode], dt);  
        }//mult.dec.ch 
    }
    else {
      LOG("NucDeEx", pNOTICE) 
          << "Hit nucleon populates continuum/other states; no photon emitted.";
      return;
    }
  } //n-hole
}
//___________________________________________________________________________
void NucDeExcitationSim::AddPhoton(
                         GHepRecord * evrec, double E0, double dt) const
{
// Add a photon at the event record & recoil the remnant nucleus so as to
// conserve energy/momenta
//
  double E = (dt>0) ? this->PhotonEnergySmearing(E0, dt) : E0;

  LOG("NucDeEx", pNOTICE) 
    << "Adding a " << E/units::MeV << " MeV photon from nucl. deexcitation";

  GHepParticle * target  = evrec->Particle(1);
  GHepParticle * remnant = 0;
  for(int i = target->FirstDaughter(); i <= target->LastDaughter(); i++) {
    remnant  = evrec->Particle(i);
    if(pdg::IsIon(remnant->Pdg())) break;
  }

  TLorentzVector x4(0,0,0,0);
  TLorentzVector p4 = this->Photon4P(E);
  GHepParticle gamma(kPdgGamma, kIStStableFinalState,1,-1,-1,-1, p4, x4);  // note that this assigns the parent of the photon as the initial-state nucleon/nucleus.  (do we want that??)
  evrec->AddParticle(gamma);  


  remnant->SetPx     ( remnant->Px() - p4.Px() );
  remnant->SetPy     ( remnant->Py() - p4.Py() );
  remnant->SetPz     ( remnant->Pz() - p4.Pz() );
  remnant->SetEnergy ( remnant->E()  - p4.E()  );
}
//___________________________________________________________________________
double NucDeExcitationSim::PhotonEnergySmearing(double E0, double dt) const
{
// Returns the smeared energy of the emitted gamma
// E0 : energy of the excited state (GeV)
// dt: excited state lifetime (sec)
// 
  double dE = kPlankConstant / (dt*units::s);

  RandomGen * rnd = RandomGen::Instance();
  double E = rnd->RndDec().Gaus(E0 /*mean*/, dE /*sigma*/);        

  LOG("NucDeEx", pNOTICE) 
     << "<E> = " << E0 << ", dE = " << dE << " -> E = " << E;

  return E;
}
//___________________________________________________________________________
TLorentzVector NucDeExcitationSim::Photon4P(double E) const
{
// Generate a photon 4p 

  RandomGen * rnd = RandomGen::Instance();

  double costheta = -1. + 2. * rnd->RndDec().Rndm();
  double sintheta = TMath::Sqrt(TMath::Max(0., 1.-TMath::Power(costheta,2)));
  double phi      = 2*kPi * rnd->RndDec().Rndm();
  double cosphi   = TMath::Cos(phi);
  double sinphi   = TMath::Sin(phi);

  double px = E * sintheta * cosphi;
  double py = E * sintheta * sinphi;
  double pz = E * costheta;

  TLorentzVector p4(px,py,pz,E);
  return p4;
}
//___________________________________________________________________________
