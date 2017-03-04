      subroutine smrad_proton(Es,theta,Ep,sig_el,sig_tot,sel,smin)
      implicit none
      real*8 Es,theta,Ep,Mp,eltail,sig_elastic_proton,varr
      real*8 peak_proton,elasrad_proton,elcor,innn
      real*8 sig_cont_proton,streg_tail
      real*8 straggling_elasic_proton,rad_tail_elastic_proton
      real*8 cor_cont_proton,t,tiw,tfw,delta_E,elra,E3
      real*8 sel,smin,sig_tot,sig_el
      real*8 dWp,nu,Q2,jack,sn2
      integer Z,A
	character*6 TG
        common/material/TG
	sn2=(sin(theta/2.0))**2
	Q2=4.0*Es*Ep*sn2
	nu=Es-Ep
c-----input values
      Mp=0.93827
      Z=1
      A=1
      E3=Es/(1.0+Es*(1.0-cos(theta))/Mp)
c-----elastic peak width==resolution
      dWp=0.001
      delta_E=(dWp+dWp**2/(2.0*Mp))/(1.0+2.0*Es/Mp*sn2)/5.
c      delta_E=0.03
c-----Target thickness
c      t=0.0188    ! 0.0175 NH3 + 0.0013 He
c      tiw=0.0012  ! 0.0012 Al
c      tfw=0.0014  ! 0.0012 Al + 0.00018 Kapton
c GDH
cccc      t=0.01342+0.011905 ! 0.6 g/cm^2 NH3 (Paraffin density assumed) + 9 mm He4
cccc      tiw=0.0008+0.00088   ! 71 um Al + 25 um Kapton
cccc      tfw=0.0008+0.00088   ! 71 um Al + 25 um Kapton

c      t=0.01342+0.00383 ! 0.6 g/cm^2 NH3 (Paraffin density assumed) + 0.36 g/cm^2 (24.9 mm) He4
c	if(TG.eq.'2mmC'.OR.TG.eq.'pureC') then
c      t=0.01218+0.00383 ! 0.52 g/cm^2 12C  + 0.36 g/cm^2 (24.9 mm) He4
c	elseif(TG.eq.'empty') then
c      t=0.00383 ! 0.36 g/cm^2 (24.9 mm) He4
c        endif
      t=0.0181              ! 0.6 cm NH3
      tiw=0.0079            ! 7.45 mm He
      tfw=0.0079            ! 7.45 mm He
c Zero
c      t=0.0d+0
c      tiw=0.0d+0
c      tfw=0.0d+0

      elra=0.0d+0
      innn=0.0d+0
      eltail=0.0d+0
c-----cross section
      if(Ep.gt.(E3-2.0*delta_E)) then
cccc      if(Ep.lt.(E3+2.0*delta_E)) then
c      varr=abs(E3-Ep)
      varr=abs(nu-Q2/(2.0*Mp))
      jack=Es/E3
      
ccc      print*,'esaltic: ',Es,E3,Ep,delta_E,theta*57.3
      
c      print*,'s1'
      
      elcor=elasrad_proton(Es,theta,t,delta_E,Z,A,tiw,tfw)
      
c      print*,'s2'
      
      elra=sig_elastic_proton(Es,theta,Mp,Z)*elcor*
     &peak_proton(varr,dWp)*jack
      
      
c      print*,'s3'
      
      innn=0.0d+0
      eltail=0.0d+0
ccc      else
ccc      elra=0.0d+0
ccc      innn=0.0d+0
ccc      eltail=0.0d+0
ccc      endif
      else
      elra=0.0d+0
      
c      print*,'s4'
      
ckp      innn=cor_cont_proton(Es,Ep,Mp,theta,t,delta_E,Z,tiw,tfw)
      
c      print*,'s5'
      
      streg_tail=straggling_elasic_proton(Z,Mp,t,Es,Ep,theta,delta_E)
ckp      eltail=rad_tail_elastic_proton(Es,Ep,Mp,theta)+streg_tail
ckp      eltail=streg_tail      !kp
ckp      eltail=rad_tail_elastic_proton(Es,Ep,Mp,theta) !kp
ccc      print*,'H tail: ',Ep,theta,streg_tail,eltail   !kp
      
c      print*,'s6'
      
      endif
      
c Not Radiative cross section
      if(Ep.gt.(E3-2.0*delta_E)) then
c      varr=abs(E3-Ep)
      varr=abs(nu-Q2/(2.0*Mp))
      jack=Es/E3
      sig_el=sig_elastic_proton(Es,theta,Mp,Z)*
     &peak_proton(varr,dWp)*jack
      endif
      sig_tot=sig_cont_proton(Es,Ep,theta)
      
      
      
      sel=elra+eltail
c     &+straggling_elasic_proton(Z,Mp,t,Es,Ep,theta,delta_E)
      smin=elra+innn+eltail
c     &+straggling_elasic_proton(Z,Mp,t,Es,Ep,theta,delta_E)
      return
      END

c##################################################################
c cole smith program for radiative correction in the elastic region			     
      real*8 function elasrad_proton(es,theta,t,delta,
     &Zvvv,Atnum,tiw,tfw)
      implicit none
      double precision spence,arg
      real*8 me,mp,pi,alpha,z,t,es,theta,delta
      real*8 cst1,eta,eel,qs
      real*8 znuc,deltac(27),del_mo,delta_t,idel
      real*8 arg11,arg15,arg19,arg23
      real*8 epr,e1,e3,e4,beta4
      real*8 radcor
      real*8 snth,tiw,tfw
      real*8 stragg_peak_proton
      integer Zvvv,Atnum
      data me/0.000511/
      data alpha/7.2993e-3/
      mp=0.93827*FLOAT(Atnum)
      pi=acos(-1.0d+0)
      znuc      = FLOAT(Zvvv)
      z		= 1.0d+0
      snth	= sin(theta)
      cst1	= 1.0-cos(theta)
      eel	= es/(1.0+es*cst1/mp)
      epr	= es+mp-eel
      e1	= es
      e3	= eel
      e4	= epr
      beta4	= sqrt(e4**2-mp**2)/e4
      eta	= es/eel
      qs	= 2.0*es*eel*cst1
      deltac(1)=28./9.-13./6.*log(qs/me**2)
      deltac(2)=(log(qs/me**2)-1.+2.*znuc*log(eta))
     &*(2.*log(e1/delta)-3.*log(eta))
      arg=(e3-e1)/e3
      deltac(3)=-spence(arg)
      deltac(4)=-znuc**2*log(e4/mp)
      deltac(5)=znuc**2*log(mp/eta/delta)*
     &(log((1.+beta4)/(1.-beta4))/beta4-2.)
      deltac(6)=znuc**2/beta4*(log((1.0+beta4)/(1.-beta4))*
     &log((e4+mp)/2.0/mp)/2.0)
      arg=sqrt((e4-mp)*(1.+beta4)/(e4+mp)/(1.-beta4))
      deltac(7)=-znuc**2*spence(-arg)/beta4
      arg=(e3-mp)/e1
      deltac(8)=znuc*spence(arg)
      arg=(mp-e3)*mp/(2.*e3*e4-mp*e1)
      deltac(9)=-znuc*spence(arg)
      arg=2.*e3*(mp-e3)/(2.0*e3*e4-mp*e1)
      deltac(10)=znuc*spence(arg)
      arg11=(2.0*e3*e4-mp*e1)/e1/(mp-2.0*e3)
      arg11=abs(arg11)
      deltac(11)=znuc*log(arg11)*log(mp/2.0/e3)
      arg=(e3-e4)/e3
      deltac(12)=-znuc*spence(arg)
      arg=(e4-e3)*mp/(2.*e1*e4-mp*e3)
      deltac(13)=znuc*spence(arg)
      arg=2.*e1*(e4-e3)/(2.0*e1*e4-mp*e3)
      deltac(14)=-znuc*spence(arg)
      arg15=(2.0*e1*e4-mp*e3)/e3/(mp-2.0*e1)
      arg15=abs(arg15)
      deltac(15)=-znuc*log(arg15)*log(mp/2.0/e1)
      arg=(e1-mp)/e1
      deltac(16)=-znuc*spence(arg)
      arg=(mp-e1)/e1
      deltac(17)=znuc*spence(arg)
      arg=2.*(mp-e1)/mp
      deltac(18)=-znuc*spence(arg)
      arg19=abs(mp/(2.0*e1-mp))
      deltac(19)=-znuc*log(arg19)*log(mp/2.0/e1)
      arg=(e3-mp)/e3
      deltac(20)=znuc*spence(arg)
      arg=(mp-e3)/e3
      deltac(21)=-znuc*spence(arg)
      arg=2.0*(mp-e3)/mp
      deltac(22)=znuc*spence(arg)
      arg23=abs(mp/(2.0*e3-mp))
      deltac(23)=znuc*log(arg23)*log(mp/2.0/e3)
      arg=(e1-e3)/e1
      deltac(24)=-spence(arg)
      arg=(e4-mp)*(1.0-beta4)/(e4+mp)/(1.0+beta4)
      arg=sqrt(arg)
      deltac(25)=znuc**2*spence(arg)/beta4
      arg=(e4-mp)/(e4+mp)
      arg=sqrt(arg)
      deltac(26)=-znuc**2*spence(arg)/beta4
      deltac(27)=znuc**2*spence(-arg)/beta4
      del_mo=0.0d+0
      do idel=1,27
      del_mo=del_mo+deltac(idel)
      enddo
      del_mo=-alpha*del_mo/pi
c------Straggling correction delta_t--------------------------
c      delta_t=stragg_peak_proton(znuc,t,es,eel,eta,delta,tiw,tfw)
      delta_t=0.0d+0
c-------------------------------------------------------------
      radcor=del_mo+delta_t
      elasrad_proton=exp(radcor)
      return
      end
c######################################################
      real*8 function stragg_peak_proton(z,t,es,eel,eta,delta,tiw,tfw)
      implicit none
      real*8 z,t,es,eel,eta,delta,bt,bfunc,b
      real*8 tiw,tfw,btiw,btfw
      b=bfunc(z)
      bt=b*t
      btiw=b*tiw
      btfw=b*tfw
      stragg_peak_proton=-((btiw+bt/2.0)*log(es/eta**2/delta)+
     &(btfw+bt/2.0)*log(eel/delta))
      return
      end
c######################################################
