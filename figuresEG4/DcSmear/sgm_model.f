      real function sgm_model(Ebeam,theta,Epr)
      implicit none
      real*8 Es,th,Ep,sig_el,sig_qe,sig_tot,
     &elra,qera,inel,eltail,qetail,smin
      real*8 sig_el_he,sel_he
      real*8 sig_el_p,sig_tot_p,sel_p,smin_p,sigman,sigmac
      real*8 He_qe,He_inel,He_el,He,proton,N15
      real Ebeam,theta,Epr
      real nu,q2,x,s,w2,ma,mp,x1,w21,rel,rel15n
      parameter (ma=12*0.93827)
      parameter (mp=0.93827)
      Es=dble(Ebeam)
      th=dble(theta)
      Ep=dble(Epr)
c Kinematic tests
      sgm_model=0.e0
      nu=Ebeam-Epr
c      if(nu.lt.0.e0) return
      q2=2.e0*Ebeam*Epr*(1.e0-cos(theta))
      x=q2/(2.e0*ma*nu)
      x1=q2/(2.e0*mp*nu)
c      if(x1.gt.1.e0) return
      s=ma**2+2.e0*Ebeam*ma
      w2=ma**2+q2*(1.e0/x-1.e0)
      w21=sqrt(mp**2+q2*(1.e0/x1-1.e0))
c      print*,'w2= ',w21,s,w2,nu,x1,q2
      
c      if(w2.lt.0.0.or.w2.gt.s) return
c      print*, q2,w21
c      call relat(q2,w21,rel)
c      print*, 's'
c      call relat_he(q2,w21,rel)
cccc Calculations for 12C ********************************************
c      call smrad(Es,th,Ep,sig_el,sig_qe,sig_tot,
c     &elra,qera,inel,eltail,qetail,smin)
cccc Calculations for 12C ********************************************

cccc Calculations for 4He ********************************************
c      call smrad_helium(Es,th,Ep,sig_el,sig_qe,sig_tot,
c     &elra,qera,inel,eltail,qetail,smin)      

c      He_qe=qera+qetail
c      He_inel=inel
c      call el_helium(Es,th,Ep,sig_el_he,sel_he)
c      He_el=sel_he
c      He=He_qe+He_inel+He_el
cccc Calculations for 4He ******************************************** 

cccc Calculations for proton ********************************************    
      call smrad_proton(Es,th,Ep,sig_el_p,sig_tot_p,sel_p,smin_p)
      proton=smin_p
cccc Calculations for proton ********************************************



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cccc Calculations for 15N ********************************************
      
c      call smrad_nitrogen(Es,th,Ep,sig_el,sig_qe,sig_tot,
c     &elra,qera,inel,eltail,qetail,smin)
cccc Calculations for 15N ********************************************

c      N15=(smin-(qera+inel+qetail))+15./12.*(qera+inel+qetail)
c      sgm_model=real(18./4.*0.3*He)
c      sgm_model=real(smin)
      sgm_model=real(proton)
c      sgm_model=real(He_el)
c       sgm_model=real(0.7*smin+0.7*3*proton+18./4.*0.3*He)
c       sgm_model=real(N15+3*proton)
c       sgm_model=real(0.7*N15+0.7*3*proton+18./4.*0.3*He)
c      sgm_model=real(0.7*smin+0.7*3*proton+18./4.*0.3*He)
c      sgm_model=real(0.7*rel*(15./12.)*smin+0.7*3*proton+18./4.*0.3*He)
      return
      end
