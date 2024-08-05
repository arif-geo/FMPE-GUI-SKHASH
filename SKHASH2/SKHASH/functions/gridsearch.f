
C FILE: GRIDSEARCH.F
      subroutine FOCALAMP_MC_WT(p_azi_mc,p_the_mc,sp_amp,p_pol,
     &    p_qual,nmc,dang,maxout,nextra,ntotal,qextra,
     &    qtotal,min_amp,npsta,
     &    nf,strike,dip,rake,faults,slips)
C           Performs grid search to find acceptable focal mechanisms,
C           for multiple trials of ray azimuths and takeoff angles.
C           Acceptable mechanisms are those with less than "ntotal"
C           misfit polarities, or the minimum plus "nextra" if this
C           is greater.
C
C  Inputs:
C           p_azi_mc(npsta,nmc)  =  azimuth to station from event (deg. E of N)
C           p_the_mc(npsta,nmc)  =  takeoff angle (from vert, up=0, <90 upgoing, >90 downgoing)
C           sp_amp(npsta)  =  amplitude ratios
C           p_pol(npsta)  =  first motion, 1=up, -1=down
C           p_qual(npsta) =  quality, used as weight, range 0 to 1
C           npsta  =  number of first motions
C           nmc    =  number of trials
C           dang   =  desired angle spacing for grid search
C           maxout =  maximum number of fault planes to return:
C                     if more are found, a random selection will be returned
C           nextra =  number of additional misfits allowed above minimum
C           ntotal =  total number of allowed misfits
C           qextra =  additional amplitude misfit allowed above minimum
C           qtotal =  total allowed amplitude misfit
C  Outputs:
C           nf     =  number of fault planes found
C           strike(min(maxout,nf)) = strike
C           dip(min(maxout,nf))    = dip
C           rake(min(maxout,nf))   = rake
C           faults(3,min(maxout,nf)) = fault normal vector
C           slips(3,min(maxout,nf))  = slip vector
C
C    Example call for creating the signature file that will be called by SKHASH:
C        python3 -m numpy.f2py -c gridsearch.f -m gridsearch
C
C    This subroutine was originally written by Hardebeck & Shearer for HASH.
C    Modifications made by Skoumal.
C    Last modified Feb 2023

C ntab = of tables
C npick0 = maximum number of picks per event
C nmc0 = maximum number of trials of location/take-off angles
C nmax0 = maximum number of acceptable mechanisms output
C ncoor = number of test mechanisms
C dang0 = minimum grid spacing, in degrees
      integer, parameter :: ntab=180, npick0=15000
      integer, parameter :: nmc0=500, nmax0=500
      integer, parameter :: ncoor=31032
      real, parameter :: dang0=5.0

C input and output arrays
      REAL, dimension(npick0,nmc0) :: p_azi_mc,p_the_mc
      INTEGER npsta, nf
      INTEGER nmc,dang,maxout
      integer p_pol(npsta)
      real p_a1(npick0),p_a2(npick0),p_a3(npick0),p_qual(npsta)
      real faultnorm(3),slip(3),faults(3,nmax0),slips(3,nmax0)
      real strike(nmax0),dip(nmax0),rake(nmax0),sp_amp(npsta)
      real min_amp
      save nrot,b1,b2,b3,thetable,phitable,amptable

C Directives for f2py
Cf2py intent(in) p_azi_mc
Cf2py intent(in) p_the_mc
Cf2py intent(in) sp_amp
Cf2py intent(in) p_pol
Cf2py intent(in) p_qual
Cf2py intent(in) nmc
Cf2py intent(in) dang
Cf2py intent(in) maxout
Cf2py intent(in) nextra
Cf2py intent(in) ntotal
Cf2py intent(in) qextra
Cf2py intent(in) qtotal
Cf2py intent(in) min_amp
Cf2py intent(in) npsta
Cf2py intent(out) nf
Cf2py intent(out) strike
Cf2py intent(out) dip
Cf2py intent(out) rake
Cf2py intent(out) faults
Cf2py intent(out) slips


C coordinate transformation arrays
      real b1(3,ncoor),bb1(3)
      real b2(3,ncoor),bb2(3)
      real b3(3,ncoor),bb3(3)
C P and S amplitude arrays
      real amptable(2,ntab,2*ntab)
      real phitable(2*ntab+1,2*ntab+1)
      real thetable(2*ntab+1)

C fit arrays
      real fit(ncoor),afit(ncoor)
      real weight_fit(ncoor),weight_afit(ncoor)
      integer irotgood(ncoor),irotgood2(ncoor)
      real fran
      real nextra,ntotal
      real qextra,qtotal
      real nmiss,ncount
      real nmissmin,nmissmax
      save irotgood,irotgood2
      save fit,afit,weight_fit,weight_afit

      pi=3.1415927
      degrad=180./pi
      if (maxout.gt.nmax0) then
        maxout=nmax0
      end if

C Set up array with direction cosines for all coordinate transformations
      irot=0
      do 5 ithe=0,int(90.1/dang)
         the=real(ithe)*dang
         rthe=the/degrad
         costhe=cos(rthe)
         sinthe=sin(rthe)
         fnumang=360./dang
         numphi=nint(fnumang*sin(rthe))
         if (numphi.ne.0) then
            dphi=360./float(numphi)
         else
            dphi=10000.
         end if
         do 4 iphi=0,int(359.9/dphi)
            phi=real(iphi)*dphi
            rphi=phi/degrad
            cosphi=cos(rphi)
            sinphi=sin(rphi)
            bb3(3)=costhe
            bb3(1)=sinthe*cosphi
            bb3(2)=sinthe*sinphi
            bb1(3)=-sinthe
            bb1(1)=costhe*cosphi
            bb1(2)=costhe*sinphi
            bb2(1)=bb3(2)*bb1(3)-bb3(3)*bb1(2)
            bb2(2)=bb3(3)*bb1(1)-bb3(1)*bb1(3)
            bb2(3)=bb3(1)*bb1(2)-bb3(2)*bb1(1)
            do 3 izeta=0,int(179.9/dang)
               zeta=real(izeta)*dang
               rzeta=zeta/degrad
               coszeta=cos(rzeta)
               sinzeta=sin(rzeta)
               irot=irot+1
               if (irot.gt.ncoor) then
                  print *,'***Error: # of rotations too big'
                  return
               end if
               b3(3,irot)=bb3(3)
               b3(1,irot)=bb3(1)
               b3(2,irot)=bb3(2)
               b1(1,irot)=bb1(1)*coszeta+bb2(1)*sinzeta
               b1(2,irot)=bb1(2)*coszeta+bb2(2)*sinzeta
               b1(3,irot)=bb1(3)*coszeta+bb2(3)*sinzeta
               b2(1,irot)=bb2(1)*coszeta-bb1(1)*sinzeta
               b2(2,irot)=bb2(2)*coszeta-bb1(2)*sinzeta
               b2(3,irot)=bb2(3)*coszeta-bb1(3)*sinzeta
3           continue
4        continue
5     continue
      nrot=irot

      astep=1./real(ntab)
      do 150 i=1,2*ntab+1
        bbb3=-1.+real(i-1)*astep
        thetable(i)=acos(bbb3)
        do 140 j=1,2*ntab+1
          bbb1=-1.+real(j-1)*astep
          phitable(i,j)=atan2(bbb3,bbb1)
          if (phitable(i,j).lt.0.) then
            phitable(i,j)=phitable(i,j)+2.*pi
          end if
140     continue
150   continue

      do 250 i=1,2*ntab
        phi=real(i-1)*pi*astep
        do 240 j=1,ntab
          theta=real(j-1)*pi*astep
          amptable(1,j,i)=abs(sin(2*theta)*cos(phi))
C Ensures amps are greater than the specified min value
          if (amptable(1,j,i).lt.min_amp) then
            amptable(1,j,i)=0.0
          end if

          s1=cos(2*theta)*cos(phi)
          s2=-cos(theta)*sin(phi)
          amptable(2,j,i)=sqrt(s1*s1+s2*s2)

C Ensures amps are greater than the specified min value
          if (amptable(2,j,i).lt.min_amp) then
            amptable(2,j,i)=0.0
          end if


240     continue
250   continue

      do irot=1,nrot
        irotgood(irot)=0
      end do

C loop over multiple trials
      do 430 im=1,nmc

C  Convert data to Cartesian coordinates
      do 40 i=1,npsta
        p_a1(i)=1*sin(p_the_mc(i,im)/degrad)*cos(p_azi_mc(i,im)/degrad)
        p_a2(i)=1*sin(p_the_mc(i,im)/degrad)*sin(p_azi_mc(i,im)/degrad)
        p_a3(i)=-1*cos(p_the_mc(i,im)/degrad)
40    continue

C  find misfit for each solution and minimum misfit
         qmissmin=99999.0
         nmissmin=99999.0
         weight_qmissmin=99999.0
         weight_nmissmin=99999.0
         do 420 irot=1,nrot
            qmiss=0.
            qcount=0.
            nmiss=0.
            ncount=0.

            do 400 ista=1,npsta
             p_b1= b1(1,irot)*p_a1(ista)
     &              +b1(2,irot)*p_a2(ista)
     &              +b1(3,irot)*p_a3(ista)
             p_b3= b3(1,irot)*p_a1(ista)
     &              +b3(2,irot)*p_a2(ista)
     &              +b3(3,irot)*p_a3(ista)
              if (sp_amp(ista).gt.0.) then

               p_proj1=p_a1(ista)-p_b3*b3(1,irot)
               p_proj2=p_a2(ista)-p_b3*b3(2,irot)
               p_proj3=p_a3(ista)-p_b3*b3(3,irot)
               plen=sqrt(p_proj1*p_proj1+p_proj2*p_proj2+
     &                    p_proj3*p_proj3)
               p_proj1=p_proj1/plen
               p_proj2=p_proj2/plen
               p_proj3=p_proj3/plen
               pp_b1=b1(1,irot)*p_proj1+b1(2,irot)*p_proj2
     &                +b1(3,irot)*p_proj3
               pp_b2=b2(1,irot)*p_proj1+b2(2,irot)*p_proj2
     &              +b2(3,irot)*p_proj3
               i=nint((p_b3+1.)/astep)+1
               theta=thetable(i)
               i=nint((pp_b2+1.)/astep)+1
               j=nint((pp_b1+1.)/astep)+1
               phi=phitable(i,j)
               i=nint(phi/(pi*astep))+1
               if (i.gt.2*ntab) i=1
               j=nint(theta/(pi*astep))+1
               if (j.gt.ntab) j=1
               p_amp=amptable(1,j,i)
               s_amp=amptable(2,j,i)
               if (p_amp.eq.0.0) then
                 sp_ratio=4.0
               else if (s_amp.eq.0.0) then
                 sp_ratio=-2.0
               else
                 sp_ratio=real(log10(4.9*s_amp/p_amp))
               end if

               qmiss=qmiss+abs(sp_amp(ista)-sp_ratio)
               qcount = qcount +1.0

             end if
             if (p_pol(ista).ne.0) then
               prod=p_b1*p_b3
               ipol=-1
               if (prod.gt.0.) ipol=1    ! predicted polarization
               if (ipol.ne.p_pol(ista)) then
                  nmiss=nmiss+p_qual(ista)
               end if
               ncount=ncount+p_qual(ista)
             end if
400         continue ! end sta loop

            fit(irot)=(nmiss) ! sum of misfit polarities
            if (ncount.eq.0.0) then
              fit(irot)=0.0
            else
              weight_fit(irot)=(nmiss/ncount)  ! weighted fraction misfit polarities
            end if

            afit(irot)=qmiss  !  sum misfit amp ratio
            if (qcount.eq.0.0) then
              weight_afit(irot)=0.0
            else
              weight_afit(irot)=qmiss/qcount !  fraction misfit amp ratio
            end if

            if ((weight_fit(irot)).lt.weight_nmissmin) then
              weight_nmissmin=weight_fit(irot)
            end if
            if ((weight_afit(irot)).lt.weight_qmissmin) then
              weight_qmissmin=weight_afit(irot)
            end if
            if (nmiss.lt.nmissmin) then
              nmissmin=nmiss
            end if

            if (qmiss.lt.qmissmin) then
              qmissmin=(qmiss)
            end if



420      continue ! end irot loop


C choose fit criteria
         weight_nmissmax=weight_ntotal
         if (weight_nmissmax.lt.(weight_nmissmin+weight_nextra)) then
            weight_nmissmax=weight_nmissmin+weight_nextra
         end if
         weight_qamissmax=weight_qtotal
         if (weight_qamissmax.lt.(weight_qmissmin+weight_qextra)) then
            weight_qamissmax=weight_qmissmin+weight_qextra
         end if


         nmissmax=ntotal
         if (nmissmax.lt.(nmissmin+nextra)) then
            nmissmax=nmissmin+nextra
         end if
         qmissmax=qtotal
         if (qmissmax.lt.(qmissmin+qextra)) then
            qmissmax=qmissmin+qextra
         end if

C loop over rotations - find those meeting fit criteria
425      nadd=0
         do irot=1,nrot
            nmiss=fit(irot)
            qmiss=afit(irot)


            if ((nmiss.le.nmissmax).and.(qmiss.le.qmissmax)) then
              irotgood(irot)=1
              nadd=nadd+1
            end if
        end do

        if (nadd.eq.0) then  ! if there are no solutions that meet
          qmissmin=99999.0  ! the criteria, loosen the amplitude criteria
          do irot=1,nrot
            nmiss=fit(irot)
            qmiss=afit(irot)
            if ((nmiss.le.nmissmax).and.
     &           (qmiss.lt.qmissmin)) then
              qmissmin=qmiss
            end if
          end do
          qmissmax=qtotal
          if (qmissmax.lt.qmissmin+qextra) then
             qmissmax=qmissmin+qextra
          end if
          goto 425
        end if

430     continue

        nfault=0
        do irot=1,nrot
          if (irotgood(irot).gt.0) then
            nfault=nfault+1
            irotgood2(nfault)=irot
          end if
        end do

C  Select output solutions
        nf=0
        if (nfault.le.maxout) then
          do i=1,nfault
            irot=irotgood2(i)
            nf=nf+1
            faultnorm(1)=b3(1,irot)
            faultnorm(2)=b3(2,irot)
            faultnorm(3)=b3(3,irot)
            slip(1)=b1(1,irot)
            slip(2)=b1(2,irot)
            slip(3)=b1(3,irot)
            do m=1,3
              faults(m,nf)=faultnorm(m)
              slips(m,nf)=slip(m)
            end do
            call FPCOOR(s1,d1,r1,faultnorm,slip,2)
            strike(nf)=s1
            dip(nf)=d1
            rake(nf)=r1
          end do
        else
          do 441 i=1,99999
            fran=rand(0)
            iscr=nint(fran*float(nfault)+0.5)
            if (iscr.lt.1) iscr=1
            if (iscr.gt.nfault) iscr=nfault
            if (irotgood2(iscr).le.0) goto 441
            irot=irotgood2(iscr)
            irotgood2(iscr)=-1
            nf=nf+1
            faultnorm(1)=b3(1,irot)
            faultnorm(2)=b3(2,irot)
            faultnorm(3)=b3(3,irot)
            slip(1)=b1(1,irot)
            slip(2)=b1(2,irot)
            slip(3)=b1(3,irot)
            do m=1,3
              faults(m,nf)=faultnorm(m)
              slips(m,nf)=slip(m)
            end do
            call FPCOOR(s1,d1,r1,faultnorm,slip,2)
            strike(nf)=s1
            dip(nf)=d1
            rake(nf)=r1
            if (nf.eq.maxout) go to 445
441       continue
445       continue
        end if

      return
      end


      subroutine FPCOOR(strike,dip,rake,fnorm,slip,idir)
C           FPCOOR gets fault normal vector,fnorm, and slip
C           vector, slip, from (strike,dip,rake) or vice versa.
C           idir = 1 compute fnorm,slip
C           idir = C compute strike,dip,rake
C           Reference:  Aki and Richards, p. 115
C           uses (x,y,z) coordinate system with x=north, y=east, z=down
      REAL strike,dip,rake
      REAL fnorm(3),slip(3),phi,del,lam,a,clam,slam
      INTEGER idir
      degrad=180./3.1415927
      pi=3.1415927
      phi=strike/degrad
      del=dip/degrad
      lam=rake/degrad
      if (idir.eq.1) then
         fnorm(1)=-sin(del)*sin(phi)
         fnorm(2)= sin(del)*cos(phi)
         fnorm(3)=-cos(del)
         slip(1)= cos(lam)*cos(phi)+cos(del)*sin(lam)*sin(phi)
         slip(2)= cos(lam)*sin(phi)-cos(del)*sin(lam)*cos(phi)
         slip(3)=-sin(lam)*sin(del)
      else
         if ((1.-abs(fnorm(3))).le.1e-7) then
           del=0.
           phi=atan2(-slip(1),slip(2))
           clam=cos(phi)*slip(1)+sin(phi)*slip(2)
           slam=sin(phi)*slip(1)-cos(phi)*slip(2)
           lam=atan2(slam,clam)
         else
           phi=atan2(-fnorm(1),fnorm(2))
           a=sqrt(fnorm(1)*fnorm(1)+fnorm(2)*fnorm(2))
           del=atan2(a,-fnorm(3))
           clam=cos(phi)*slip(1)+sin(phi)*slip(2)
           slam=-slip(3)/sin(del)
           lam=atan2(slam,clam)
           if (del.gt.(0.5*pi)) then
             del=pi-del
             phi=phi+pi
             lam=-lam
           end if
         end if
         strike=phi*degrad
         if (strike.lt.0.) strike=strike+360.
         dip=del*degrad
         rake=lam*degrad
         if (rake.le.-180.) rake=rake+360.
         if (rake.gt.180.) rake=rake-360.
      end if
      return
      end
