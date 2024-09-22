!==================================================================================================
! One step of base irradiation
!==================================================================================================
      subroutine onestep(time, dt, nnn)

      use globals
      implicit none
      integer, parameter :: mcol=5*maxr+2
      integer nnn,j,i,ieh,ier,iez,ii,indxa(mcol),ir,ish,isr,isz,jj,k
      real*8 time,dt,pi,cpcool,area,tmp(maxtab),polp,tdwn,pow(maxz),kcool,rocool,vicool, &
      &      htc_l,htc_t,tc,pecool,fff,xx,ttt,kclad_,qss,rrr,vvv,clamb,qf(maxr,maxz),    &
      &      error,hgap,err1,sto,kfuel_,ak,an,d,de,dsign,hclad,hclad0, &
      &      hfuel,hfuel0,pss,sg,sgfr,sgh,sgr,sgz,sss,vf,vf0,yng,gaphtc,flamb,ftexp,fdens, &
      &      fswel,felmod,fsigf,fcreep,ctexp,cswel,celmod,ccreep,a(mcol,mcol),b(mcol), &
      &      fpoir,cpoir,ql2
      
      pi=2.0d0*acos(0.d0)

!     save some variables at the beginning of the time step
      do j=1,nzz
!        save burnup and fuel swelling
         bup0(j)=bup(j)
         do i=1,nf-1
            efs0(i,j)=efs(i,j)
         end do
!        save amount of fission gases generated and released
         fggen0=fggen
         fgrel0=fgrel
!        save fuel outer surface hoop and axial strains
         efh0(j)=efh(nf,j)
         efz0(j)=efz(j)
!        save inner clad surface hoop and axial strains
         eh0(j)=eh(1,j)
         ez0(j)=ez(j)
!        save average clad temperature (for clad failure estimate) 
         tcave0(j)=tcave(j)
!        save average fuel temperature
         tfave0(j)=tfave(j)
         do i=1,nc-1
!           save clad plastic strains
            epe0(i,j)=epe(i,j)
            ehp0(i,j)=ehp(i,j)
            ezp0(i,j)=ezp(i,j)
            erp0(i,j)=erp(i,j)
         end do
         do i=1,nc-1
!           save clad creep strains
            ece0(i,j)=ece(i,j)
            ehc0(i,j)=ehc(i,j)
            ezc0(i,j)=ezc(i,j)
            erc0(i,j)=erc(i,j)
         end do
         do i=1,nf-1
!           save fuel plastic strains
            epef0(i,j)=epef(i,j)
            ehpf0(i,j)=ehpf(i,j)
            ezpf0(i,j)=ezpf(i,j)
            erpf0(i,j)=erpf(i,j)
         end do
         do i=1,nf-1
!           save fuel creep strains
            ecef0(i,j)=ecef(i,j)
            ehcf0(i,j)=ehcf(i,j)
            ezcf0(i,j)=ezcf(i,j)
            ercf0(i,j)=ercf(i,j)
         end do
      end do

!     calculate coolant specific heat
      if(ctype.eq.'he')then
         cpcool=5193.d0
      else if(ctype.eq.'pbbi')then
         cpcool=146.5d0
      else if(ctype.eq.'pb')then
         cpcool=147.3d0
      else if(ctype.eq.'na')then
         cpcool=1270.0d0
      else
         cpcool=0.0
         write(*,*)'coolant not available for base irradiation:',ctype
      end if

!     fuel pellet cross section area
      area=pi*(rfo0**2-rfi0**2)
      do j=1,nzz
!        power density
         do i=1,nqv
            tmp(i)=qv(j,i)
         end do
         qqv1(j)=polp(time,tmp,tqv,nqv)*polp(time,qvmlt,tqvmlt,nqvmlt)
!        linear power (W/m) 
         ql(j)=qqv1(j)*area
!        power (W) 
         pow(j)=ql(j)*dz0(j)  
!        calculate fuel burnup (MWd/kgU)
         bup(j)=bup0(j)+(qqv1(j)/1.d6)/rof0*(dt/8.64d4)
!        calculate coolant temperature
         if(j.eq.1)then
            tcool(j)=tcoolin+0.5*pow(j)/flowr/cpcool
         else
            tdwn=tcool(j-1)
            tcool(j)=tcool(j-1)+0.5*(pow(j-1)+pow(j))/flowr/cpcool
         end if

!        calculate coolant properties
         if(ctype.eq.'he')then
!           thermal conductivity
            kcool=2.639d-3*tcool(j)**0.7085
!           density
            rocool=press(j)*4.d-3/rmu/tcool(j)
            vicool=(7.0019d-6+4.5244d-8*tcool(j)-6.28d-12*tcool(j)**2)/rocool
!           velocity
            vcool(j)=flowr/rocool/xarea
      
!           clad-to-coolant heat transfer coefficient
!           1. laminar forced convection
            htc_l=4.36*kcool/dhyd
      
!           2. turbulent forced convection
            htc_t=0.023d0*kcool/dhyd*(rocool*vcool(j)*dhyd/vicool)**0.8*(vicool*cpcool/kcool)**0.4
            htc(j)=dmax1(htc_l,htc_t)
      
         else if(ctype.eq.'pbbi' .or. ctype.eq.'pb')then
!           thermal conductivity
            if(ctype.eq.'pbbi')then
               kcool=3.9021d0+0.0123d0*tcool(j)
!              density
               rocool=10735.0 - 1.375*(tcool(j)-273.15d0)
            else
               tc=tcool(j)-273.15d0
               if(tc.le.550.0d0)then
                  kcool=14.30d0+0.0023d0*tc
               else
                  kcool=34.89d0-0.0730d0*tc+6.9d-5*tc**2
               end if
!              density
               rocool=11072.0d0-1.2d0*tc
            end if
!           velocity
            vcool(j)=flowr/rocool/xarea
            pecool=vcool(j)*dhyd*rocool*cpcool/kcool
            if(pecool.lt.4000.0d0) then
               htc(j)=kcool/dhyd*(5.0d0+0.025*pecool**0.8)
            else
               htc(j)=kcool/dhyd*(7.5d0+0.005d0*pecool)
            end if
         else if(ctype.eq.'na')then
!           thermal conductivity
            kcool=1.1641922e2*EXP(-7.0184779e-4*tcool(j))
!           density
            fff=1.0-tcool(j)/2503.7
            rocool=219.0+275.32*fff+511.58*fff**0.5
!           velocity
            vcool(j)=flowr/rocool/xarea
            pecool=vcool(j)*dhyd*rocool*cpcool/kcool
!           pitch-to-diameter ratio
            xx=1.1
!           clad-to-coolant heat transfer coefficient
            htc(j)=kcool/dhyd*0.047*(1.0-exp(-3.8*(xx-1.0)))*(pecool**0.77 + 250.0)
         end if
!
!        calculate outer cladding temperature
         if(ntco .gt. 0)then
            do i=1,ntco
               tmp(i)=tco(j,i)
            end do
            tem(nf+nc,j) = polp(time,tmp,ttco,ntco)
         else
!           heat flux from the clad outer surface
            qs(j)=pow(j)/(2.d0*pi*rco0*dz0(j))
            tem(nf+nc,j) = tcool(j)+qs(j)/htc(j)
         end if
!        calculate other clad temperatures
         do i=nf+nc-1,nf+1,-1
            ttt=0.5*(tem(i,j)+tem(i+1,j))
            kclad_=clamb(ttt,cmat)
            qss=pow(j)/(2.d0*pi*rad_0(i)*dz0(j))
            tem(i,j)=tem(i+1,j)+qss*drc0/kclad_
         end do
!        calculate heat flux between fuel and cladding
         rrr=(rfo0+rci0)/2.0
         qf(nf,j)=pow(j)/(2.d0*pi*rrr*dz0(j))
!        calculate other heat fluxes calculated between nodes of pellet
         do i=nf-1,1,-1
            if(i .eq. nf-1)then
               vvv = rfo0**2 - rad_0(i)**2
            else
               vvv = rad_0(i+1)**2 - rad_0(i)**2
            end if
            qf(i,j)=(qf(i+1,j)*rad_0(i+1) - qqv1(j)*vvv/2.0)/rad_0(i)
         end do
      end do

!     fission gas release
      call fgr()

!     gas plenum temperature
      tple=tcoolin

      error = 2.0
      niter = 0
!     outer iteration loop
      do while(error > 1.0)
         niter = niter + 1
         error = 0.0
         if(niter > 49)then
            !write(*,*)'Too many iterations.'
            go to 50
            !stop
         end if
         do j=1,nzz
!--------------------------------------------------------------------------------------------------
!           1. Given gap + pfc, calculate hgap
!           pellet-clad gap heat conductance
            hgap = gaphtc(j)
!            err1 = min(abs((hgapt(j)-hgap)/hgap/rtol),
!                       abs((hgapt(j)-hgap)/atol))
!            if(err1 > error)dueto = 'hgap'
!            error = max(error,err1)
            hgapt(j)=hgap
            
!--------------------------------------------------------------------------------------------------
!           2. Given hgap, calculate fuel temperatures
!           outer fuel temperature
            ttt = tem(nf+1,j)+qf(nf,j)/hgapt(j)
            err1 = min(abs((tem(nf,j)-ttt)/ttt/rtol), abs((tem(nf,j)-ttt)/atol))
            if(err1 > error)dueto = 'tem'
            error = max(error,err1)
            tem(nf,j)=ttt
!           table of fuel stoichiometry vs fuel burnup (MWd/kg)
            if(nfsto.gt.0)then
               sto=polp(bup(j),stob,bsto,nfsto)
            else
               sto=sto0
            end if   
            
!           other fuel temperatures
            do i=nf-1,1,-1
               ttt=0.5*(tem(i,j)+tem(i+1,j))
               kfuel_=flamb(rof0,bup(j),ttt, pucont,sto,fmat)
               ttt = tem(i+1,j) + qf(i,j)*drf0/kfuel_
               err1 = min(abs((tem(i,j)-ttt)/ttt/rtol), abs((tem(i,j)-ttt)/atol))
               if(err1 > error)dueto = 'tem'
               error = max(error,err1)
               tem(i,j) = ttt
            end do

!           fuel average temperature
            tfave(j)=0.d0
            do i=1,nf
               tfave(j)=tfave(j)+tem(i,j)*az0(i)
            end do
            tfave(j)=tfave(j)/pi/(rfo0**2-rfi0**2)
!           clad average temperature
            tcave(j)=0.d0
            do i=nf+1,nf+nc
               tcave(j)=tcave(j)+tem(i,j)*az0(i)
            end do
            tcave(j)=tcave(j)/pi/(rco0**2-rci0**2)

!--------------------------------------------------------------------------------------------------
!           3. Given fuel strains (efh, efr, efz), calculate rad, drf, dzf
            do i=1,nf
!              fuel nodes radii
               rad(i,j)=rad0(i)*(1.d0+efh(i,j))
               if(rad(i,j).lt.0.0)then
                  write(*,*)rad(i,j), ': rad(i,j) < 0 (',i,j,')'
                  stop
               end if
            end do
            do i=1,nf-1
!             fuel nodes thickness
              drf(i,j)=drf0*(1.d0+efr(i,j))
              if(drf(i,j).lt.0.0)then
                 write(*,*)drf(i,j), ': drf(i,j) < 0 (',i,j,')'
                 stop
              end if
            end do
            dzf(j)=dz0(j)*(1.d0+efz(j))
            if(dzf(j).lt.0.0)then
               write(*,*)dzf(j), ': dzf(j) < 0 (',j,')'
               stop
            end if
            rfi(j)=rad(1,j)
            rfo(j)=rad(nf,j)
            
!--------------------------------------------------------------------------------------------------
!           4. Given clad strains (eh, er, ez), calculate rad, drc, dzc
            do i=1,nc
!              clad nodes radii
               rad(nf+i,j)=rad0(nf+i)*(1.d0+eh(i,j))
               if(rad(nf+i,j).lt.0.0)then
                  write(*,*)rad(nf+i,j), ': rad(nf+i,j) < 0 (',i,j,')'
                  stop
               end if
            end do
            do i=1,nc-1
!              clad nodes thickness
               drc(i,j)=drc0*(1.d0+er(i,j))
               if(drc(i,j).lt.0.0)then
                  write(*,*)drc(i,j), ': drc(i,j) < 0 (',i,j,')'
                  stop
               end if
            end do
            dzc(j)=dz0(j)*(1.d0+ez(j))
            if(dzc(j).lt.0.0)then
               write(*,*)dzc(j), ': dzc(j) < 0 (',j,')'
               stop
            end if
            rci(j)=rad(nf+1,j)
            rco(j)=rad(nf+nc,j)

!--------------------------------------------------------------------------------------------------
!           5. Given rad, calculate gap width (gap), gap state (flag) and contact pressure (pfc)
            gap(j)=rci(j)-rfo(j)
            if(gap(j) .le. ruff+rufc)flag(j) = 'clos'

!           contact pressure
            if(flag(j) .eq. 'clos')pfc(j)=dabs(sigr(1,j))
            
            if(pfc(j).lt.gpres)then
               flag(j)='open'
!              contact pressure
               pfc(j)=0.d0
            end if

!--------------------------------------------------------------------------------------------------
!           6. Given tem, rad and sig, calculate fuel strain components (eft, efd, efs, ehpf, erpf, ezpf, ehcf, ercf, ezcf) and effective stress (sigf)
            do i=1,nf-1
               ttt=0.5d0*(tem(i,j)+tem(i+1,j))
!              adjust fuel pellet density
               vf0=(rad0(i+1)**2-rad0(i)**2)*dz0(j)
               vf=(rad(i+1,j)**2-rad(i,j)**2)*dzf(j)
               rof(i,j)=rof0*vf0/vf
!              fuel thermal expansion
               eft(i,j)=ftexp(ttt,pucont,fmat)
!              fuel densification
               efd(i,j)=fdens(rof(i,j),pucont,tsint, bup(j))
!              fuel swelling 
               if(ifswel.ne.0)then
                  efs(i,j)=efs0(i,j)+fswel(i,j)/3.0d0
               end if                        
               sgr=0.5d0*(sigfr(i,j)+sigfr(i+1,j))
               sgh=sigfh(i,j)
               sgz=sigfz(i,j)
               sss=(sgh+sgz+sgr)/3.d0
!              elasticity modulus 
               yng=felmod(ttt,rof(i,j),pucont,fmat)
!              effective stress
               sg=dsqrt((sgh-sgz)**2+(sgh-sgr)**2+(sgz-sgr)**2)/sqrt(2.0)
               sigf(i,j)=sg
!              fuel plasticity
               if(ifplas .ne. 0 .and. sg .gt. 0.0)then
                  sgfr=fsigf(ttt,pucont,rof(i,j),fmat)
                  de=dmax1((sg-sgfr)/yng,0.0d0)
                  epef(i,j)=epef0(i,j)+de
!                 Prandtl-Reuss flow rules
                  ehpf(i,j)=ehpf0(i,j)+1.5*de*(sgh-sss)/sg
                  ezpf(i,j)=ezpf0(i,j)+1.5*de*(sgz-sss)/sg
                  erpf(i,j)=erpf0(i,j)+1.5*de*(sgr-sss)/sg
               end if
!              fuel creep 
               if(ifcreep .ne.0 .and. sg .gt. 0.0)then
                  de=fcreep(i,j)*dt
                  ecef(i,j)=ecef0(i,j)+de
!                 Prandtl-Reuss flow rules
                  ehcf(i,j)=ehcf0(i,j)+1.5*de*(sgh-sss)/sg
                  ezcf(i,j)=ezcf0(i,j)+1.5*de*(sgz-sss)/sg
                  ercf(i,j)=ercf0(i,j)+1.5*de*(sgr-sss)/sg
               end if
            end do

!--------------------------------------------------------------------------------------------------
!           7. Given tem, rad and sig, calculate clad strain components (et, ecs, ehp, erp, ezp, ehc, erc, ezc) and effective stress (sig)
            
!           adjust clad density
            roc(j)=roc0*(rco0**2-rci0**2)*dz0(j)/((rco(j)**2-rci(j)**2)*dzc(j))
            do i=1,nc-1
               ttt=0.5d0*(tem(nf+i,j)+tem(nf+i+1,j))
!              clad thermal expansion
               et(i,j)=ctexp(ttt,cmat)
!              clad swelling
               if(icswel.ne.0)then
                  ecs(i,j)=cswel(icswel,ttt,bup(j),cmat)/3.0d0
               end if   
               sgr=0.5d0*(sigr(i,j)+sigr(i+1,j))
               sgh=sigh(i,j)
               sgz=sigz(i,j)
               sss=(sgh+sgz+sgr)/3.d0
!              elasticity modulus 
               yng=celmod(ttt,cmat)
!              effective stress
               sg=dsqrt((sgh-sgz)**2+(sgh-sgr)**2+(sgz-sgr)**2)/sqrt(2.0)
               sig(i,j)=sg
!              clad plasticity
               if(icplas .ne.0 .and. sg .gt. 0.0)then
!                 parameters of stress-strain diagram
                  call ckmn(ttt,cmat,yng,rtol,ak,an,sigy(i,j),sigb(j),uelon(j))
                  dsign=sg-ak*(epe(i,j)+sg/yng)**an
                  de=dmax1(dsign/yng,0.0d0)
                  epe(i,j)=epe0(i,j)+de
!                 Prandtl-Reuss flow rules
                  ehp(i,j)=ehp0(i,j)+1.5*de*(sgh-sss)/sg
                  ezp(i,j)=ezp0(i,j)+1.5*de*(sgz-sss)/sg
                  erp(i,j)=erp0(i,j)+1.5*de*(sgr-sss)/sg
               end if
!              clad creep
               if(iccreep .ne.0 .and. sg .gt. 0.0)then
!                 effective clad creep increment
                  de=ccreep(ttt,sg,cmat)*dt
                  ece(i,j)=ece0(i,j)+de
!                 Prandtl-Reuss flow rules
                  ehc(i,j)=ehc0(i,j)+1.5*de*(sgh-sss)/sg
                  ezc(i,j)=ezc0(i,j)+1.5*de*(sgz-sss)/sg
                  erc(i,j)=erc0(i,j)+1.5*de*(sgr-sss)/sg
               end if
            end do

!--------------------------------------------------------------------------------------------------
!           8. Given rad, drf, drc, gpres, fuel strain components (eft, efd, efs, ehpf, erpf, ezpf, ehcf, ercf, ezcf) and clad strain components (et, ecs, ehp, erp, ezp, ehc, erc, ezc), 
!              calculate total fuel strains (efh, efr, efz), clad strains (eh, er, ez), fuel stress (sigfh, sigfr, sigfz) and clad stress (sigh, sigr, sigz)
            
!           calculate fuel and clad stresses and strains
            do ii=1,mcol
            do jj=1,mcol
               a(ii,jj)=0.d0
            end do
            end do
            
            ir=0
            
!           FUEL...
            ieh=ir
            ier=ieh+nf
            iez=ier+nf-1
            ish=iez+1
            isr=ish+nf-1
            isz=isr+nf
!           FUEL: Hooke's eq. for eh:yng*eh-sigh+pss*sigr+pss*sigz=yng*eh_pct
            do i=1,nf-1
               ir=ir+1
               ttt=tem(i,j)
               yng=felmod(ttt,rof(i,j),pucont,fmat)
               pss=fpoir(fmat)
               a(ir,ieh+i)=0.5d0*yng
               a(ir,ieh+i+1)=0.5d0*yng
               a(ir,ish+i)=-1.0d0
               a(ir,isr+i)=0.5d0*pss
               a(ir,isr+i+1)=0.5d0*pss
               a(ir,isz+i)=pss
               b(ir)=yng*(eft(i,j)+efd(i,j)+efs(i,j)+ehcf(i,j)+ehpf(i,j))
            end do
!           FUEL: Hooke's eq. for er:yng*er+pss*sigh-sigr+pss*sigz=yng*er_pct
            do i=1,nf-1
               ir=ir+1
               ttt=tem(i,j)
               yng=felmod(ttt,rof(i,j),pucont,fmat)
               pss=fpoir(fmat)
               a(ir,ier+i)=yng
               a(ir,isr+i)=-0.5d0
               a(ir,isr+i+1)=-0.5d0
               a(ir,ish+i)=pss
               a(ir,isz+i)=pss
               b(ir)=yng*(eft(i,j)+efd(i,j)+efs(i,j)+ercf(i,j)+erpf(i,j))
            end do
!           FUEL: Hooke's eq. for ez:yng*ez+pss*sigh+pss*sigr-sigz=yng*ez_pct
            do i=1,nf-1
               ir=ir+1
               ttt=tem(i,j)
               yng=felmod(ttt,rof(i,j),pucont,fmat)
               pss=fpoir(fmat)
               a(ir,iez+1)=yng
               a(ir,isz+i)=-1.d0
               a(ir,ish+i)=pss
               a(ir,isr+i)=0.5d0*pss
               a(ir,isr+i+1)=0.5d0*pss
               b(ir)=yng*(eft(i,j)+efd(i,j)+efs(i,j)+ezcf(i,j)+ezpf(i,j))
            end do
!           FUEL: strain compatibility eq.:(rr-rl)*er+rl0*ehl-rr0*ehr=0
            do i=1,nf-1
               ir=ir+1
               a(ir,ier+i)=drf(i,j)
               a(ir,ieh+i)=rad(i,j)
               a(ir,ieh+i+1)=-rad(i+1,j)
               b(ir)=0.d0
            end do
!           FUEL: stress equilibrity eq.:(rr-rl)*sigh-rl0*sigrl+rr0*sigrr=0
            do i=1,nf-1
               ir=ir+1
               a(ir,ish+i)=drf(i,j)
               a(ir,isr+i)=rad(i,j)
               a(ir,isr+i+1)=-rad(i+1,j)
               b(ir)=0.d0
            end do
!           FUEL: axial stress equilibrity eq.:
            ir=ir+1
            if(flag(j).eq.'clos')then
!              sum(sigz*da)=press**rco**2
               do i=1,nf-1
                  a(ir,isz+i)=rad(i+1,j)**2-rad(i,j)**2
               end do
               do i=1,nc-1
                  k = (isz+nf-1)+nc+(nc-1)+1+(nc-1)+nc+i
                  a(ir,k)=rad(nf+i+1,j)**2-rad(nf+i,j)**2
               end do
               b(ir)=-press(j)*rco(j)**2
            else
!              open gap
!              sum(sigz*da)=-gpress**rfo**2
               do i=1,nf-1
                  a(ir,isz+i)=rad(i+1,j)**2-rad(i,j)**2
               end do
               b(ir)=-gpres*(rfo(j)**2-rfi(j)**2)
            end if
            
!           FUEL: boundary conditions
            ir=ir+1
            if(rad0(1).eq.0.d0)then
!              no central hole. symmetry: sigr=sigh
               a(ir,isr+1)=1.d0
               a(ir,ish+1)=-1.d0
               b(ir)=0.0d0
            else
!              central hole. gas pressure: sigr=-gpres
               a(ir,isr+1)=1.d0
               b(ir)=-gpres
            end if
            ir=ir+1
            if(flag(j).eq.'open')then
!              open gap: sigr=-gpres
               a(ir,isr+nf)=1.d0
               b(ir)=-gpres
            else
!              closed gap: D(ez_fuel)=D(ez_clad)
               a(ir,iez+1)=1.d0
               a(ir,(isz+nf-1)+nc+(nc-1)+1)=-1.d0
               b(ir)= efz0(j)-ez0(j)
            end if
            
!           CLAD...
            ieh=ir
            ier=ieh+nc
            iez=ier+nc-1
            ish=iez+1
            isr=ish+nc-1
            isz=isr+nc
!           CLAD: Hooke's eq. for eh:yng*eh-sigh+pss*sigr+pss*sigz=yng*eh_pct
            do i=1,nc-1
               ir=ir+1
               ttt=tem(nf+i,j)
               yng=celmod(ttt,cmat)
               pss=cpoir(cmat)
               a(ir,ieh+i)=0.5d0*yng
               a(ir,ieh+i+1)=0.5d0*yng
               a(ir,ish+i)=-1.0d0
               a(ir,isr+i)=0.5d0*pss
               a(ir,isr+i+1)=0.5d0*pss
               a(ir,isz+i)=pss
               b(ir)=yng*(ehp(i,j)+ehc(i,j)+et(i,j)+ecs(i,j))
            end do
!           CLAD: Hooke's eq. for er:yng*er+pss*sigh-sigr+pss*sigz=yng*er_pct
            do i=1,nc-1
               ir=ir+1
               ttt=tem(nf+i,j)
               yng=celmod(ttt,cmat)
               pss=cpoir(cmat)
               a(ir,ier+i)=yng
               a(ir,isr+i)=-0.5d0
               a(ir,isr+i+1)=-0.5d0
               a(ir,ish+i)=pss
               a(ir,isz+i)=pss
               b(ir)=yng*(erp(i,j)+erc(i,j)+et(i,j)+ecs(i,j))
            end do
!           CLAD: Hooke's eq. for ez:yng*ez+pss*sigh+pss*sigr-sigz=yng*ez_pct
            do i=1,nc-1
               ir=ir+1
               ttt=tem(nf+i,j)
               yng=celmod(ttt,cmat)
               pss=cpoir(cmat)
               a(ir,iez+1)=yng
               a(ir,isz+i)=-1.d0
               a(ir,ish+i)=pss
               a(ir,isr+i)=0.5d0*pss
               a(ir,isr+i+1)=0.5d0*pss
               b(ir)=yng*(ezp(i,j)+ezc(i,j)+et(i,j)+ecs(i,j))
            end do
!           CLAD: strain compatibility eq.:(rr-rl)*er+rl0*ehl-rr0*ehr=0
            do i=1,nc-1
               ir=ir+1
               a(ir,ier+i)=drc(i,j)
               a(ir,ieh+i)=rad(nf+i,j)
               a(ir,ieh+i+1)=-rad(nf+i+1,j)
               b(ir)=0.d0
            end do
!           CLAD: stress equilibrity eq.:(rr-rl)*sigh-rl0*sigrl+rr0*sigrr=0
            do i=1,nc-1
               ir=ir+1
               a(ir,ish+i)=drc(i,j)
               a(ir,isr+i)=rad(nf+i,j)
               a(ir,isr+i+1)=-rad(nf+i+1,j)
               b(ir)=0.d0
            end do
!           CLAD: axial stress equilibrity eq.:
            ir=ir+1
            if(flag(j).eq.'clos')then
!              sigr_fuel=sigr_clad
               a(ir,isr+1)=1.d0
               a(ir,nf+(nf-1)+1+(nf-1)+nf)=-1.d0
               b(ir)=0.d0
            else
!              open gap: sum(sigz*da)=gpres*rci**2-press**rco**2
               do i=1,nc-1
                  a(ir,isz+i)=rad(nf+i+1,j)**2-rad(nf+i,j)**2
               end do
               b(ir)=gpres*rci(j)**2-press(j)*rco(j)**2
            end if
!           CLAD: boundary conditions
            ir=ir+1
            a(ir,isr+nc)=1.d0
            b(ir)=-press(j)
            ir=ir+1
            if(flag(j).eq.'open')then
!              open gap: sigr=-gpres
               a(ir,isr+1)=1.d0
               b(ir)=-gpres
            else
!              closed gap: D(eh_fuel)=D(eh_clad)
               a(ir,ieh+1)=-1.d0
               a(ir,nf)=1.d0
               b(ir)=efh0(j)-eh0(j)
            end if
            
            call ludcmp(a,ir,mcol,indxa,d)
            call lubksb(a,ir,mcol,indxa,b)
            
!           FUEL: read the solution
            ir=0
            do i=1,nf
               ir=ir+1
               err1 = min(abs((efh(i,j)-b(ir))/b(ir)/rtol),abs((efh(i,j)-b(ir))/atol))
               if(err1 > error)dueto = 'efh'
               error = max(error,err1)
               efh(i,j)=b(ir)
            end do
            do i=1,nf-1
               ir=ir+1
               err1 = min(abs((efr(i,j)-b(ir))/b(ir)/rtol),abs((efr(i,j)-b(ir))/atol))
               if(err1 > error)dueto = 'efr'
               error = max(error,err1)
               efr(i,j)=b(ir)
            end do
            ir=ir+1
            err1 = min(abs((efz(j)-b(ir))/b(ir)),abs((efz(j)-b(ir)/atol)))
            if(err1 > error)dueto = 'efz'
            error = max(error,err1)
            efz(j)=b(ir)
            do i=1,nf-1
               ir=ir+1
!               err1 = min(abs((sigfh(i,j)-b(ir))/b(ir)/rtol),abs((sigfh(i,j)-b(ir))/atol))
!               if(err1 > error)dueto = 'sigfh'
!               error = max(error,err1)
               sigfh(i,j)=b(ir)
            end do
            do i=1,nf
               ir=ir+1
!               err1 = min(abs((sigfr(i,j)-b(ir))/b(ir)/rtol),abs((sigfr(i,j)-b(ir))/atol))
!               if(err1 > error)dueto = 'sigfr'
!               error = max(error,err1)
               sigfr(i,j)=b(ir)
            end do
            do i=1,nf-1
               ir=ir+1
!               err1 = min(abs((sigfz(i,j)-b(ir))/b(ir)/rtol),abs((sigfz(i,j)-b(ir))/atol))
!               if(err1 > error)dueto = 'sigfz'
!               error = max(error,err1)
               sigfz(i,j)=b(ir)
            end do
            
!           CLAD: read the solution
            do i=1,nc
               ir=ir+1
               err1 = min(abs((eh(i,j)-b(ir))/b(ir)/rtol),abs((eh(i,j)-b(ir))/atol))
               if(err1 > error)dueto = 'eh'
               error = max(error,err1)
               eh(i,j)=b(ir)
            end do
            do i=1,nc-1
               ir=ir+1
               err1 = min(abs((er(i,j)-b(ir))/b(ir)/rtol),abs((er(i,j)-b(ir))/atol))
               if(err1 > error)dueto = 'er'
               error = max(error,err1)
               er(i,j)=b(ir)
            end do
            ir=ir+1
            err1 = min(abs((ez(j)-b(ir))/b(ir)/rtol),abs((ez(j)-b(ir))/atol))
            if(err1 > error)dueto = 'ez'
            error = max(error,err1)
            ez(j)=b(ir)
            do i=1,nc-1
               ir=ir+1
               err1 = min(abs((sigh(i,j)-b(ir))/b(ir)/rtol),abs((sigh(i,j)-b(ir))/atol))
               if(err1 > error)dueto = 'sigh'
               error = max(error,err1)
               sigh(i,j)=b(ir)
            end do
            do i=1,nc
               ir=ir+1
               err1 = min(abs((sigr(i,j)-b(ir))/b(ir)/rtol),abs((sigr(i,j)-b(ir))/atol))
               if(err1 > error)dueto = 'sigr'
               error = max(error,err1)
               sigr(i,j)=b(ir)
            end do
            do i=1,nc-1
               ir=ir+1
               err1 = min(abs((sigz(i,j)-b(ir))/b(ir)/rtol),abs((sigz(i,j)-b(ir))/atol))
               if(err1 > error)dueto = 'sigz'
               error = max(error,err1)
               sigz(i,j)=b(ir)
            end do

         end do !end of loop over axial layers

!        calculation of clad failure criteria
!         call cfail()
!--------------------------------------------------------------------------------------------------
!        Given rad and tem, calculate inner gas pressure (gpres)
         call precal()
      end do ! end of outer iteration loop
50    continue
      
!     initial and new fuel stack height
      hfuel0=0.d0
      hfuel=0.d0
      do j=1,nzz
         hfuel0=hfuel0+dz0(j)
         hfuel=hfuel+dzf(j)
      end do
!     fuel stack elongation
      dlfuel=hfuel/hfuel0-1.d0
      
!     initial and new clad height
      hclad0=0.d0
      hclad=0.d0
      do j=1,nzz
         hclad0=hclad0+dz0(j)
         hclad=hclad+dzc(j)
      end do
!     clad elongation
      dlclad=hclad/hclad0-1.d0

      write(*,'(a19,1pe12.5,a3,1pe12.5,a3,i3,a14,a10)')'>>>FRED BASE IRRAD ',time/3600./24.,' d ',dt/3600./24.,' d ', &
     &       niter, ' iters due to ', dueto
      ql2=0.0
      fff=0.0
      do j=1,nzz
         ql2=ql2+ql(j)/1.d3*dz0(j)
         fff=fff+dz0(j)
      end do
      ql2=ql2/fff
      j=jaxial
      write(702,'(20(1x,1pe12.5),1x,a10)') &
     &      time, &
     &      time/86400., &
     &      tem(1,j)-273.15, &
     &      tem(nf,j)-273.15, &
     &      tem(nf+1,j)-273.15, &
     &      tem(nf+nc,j)-273.15, &
     &      tfave(j)-273.15, &
     &      tcave(j)-273.15, &
     &      hgapt(j), &
     &      rfi(j)*1.d3, &
     &      rfo(j)*1.d3, &
     &      rci(j)*1.d3, &
     &      rco(j)*1.d3, &
     &      gap(j)*1.d6, &
     &      gapth(j)*1.d6, &
     &      gpres/1.d6, &
     &      pfc(j)/1.d6, &
     &      htc(j), &
     &      fgrpp, &
     &      efs(nf-1,j)*100., &
     &      flag(j)

      call outfrd(time,nnn)

      return
      end

!==================================================================================================
! Returns clad material effective creep rate (1/s)
!     input arguments :
!        tk -- temperature (k)
!        nflux -- neutron flux (n/cm^2*s)
!        mnen -- mean neutron energy (MeV)
!        sg -- effective stress (pa)
!        cmat -- clad material
!==================================================================================================
      real*8 function ccreep(tk,sg,cmat)
      implicit none
      real*8 tk,sg
      character*6 cmat
     if(cmat.eq.'aim1')then
!        B.K. Choudhary and E. Isaac Samuel, Journal of Nuclear Materials, 412 (2011)
!        effective creep = Thermal creep + Irradiaion induced creep
         ccreep = 2.3d14*dexp(-84600.0d0/(1.986d0*tk))*dsinh((39.72d0*(sg/1.d6))/(1.986d0*tk))
!         ccreep = ccreep + (3.2d-24*mnen*nflux*(sg/1.d6))
      else
         ccreep=4.5d4*(sg/1.d6)**12.9*dexp(-6.21d5/8.314d0/tk)
      end if
      if(ccreep.gt.1.0d-4)then
         ccreep=1.0d-4
      end if
      ccreep=ccreep/3.6d3
      return
      end

!==================================================================================================
! Calculation of the clad Young's modulus
!                unit: Pa
!     input arguments :
!        tk .... temperature (K)
!==================================================================================================
      real*8 function celmod(tk,cmat)
      implicit none
      real*8 tk,tc
      character*6 cmat

      celmod = 0.0
      if(cmat.eq.'ss823')then
         celmod=2.705d10-1.143d7*tk
      else if(cmat.eq.'aim1')then
         tc=tk-273.15
         celmod=2.027d11-8.167d7*tc
      else if(cmat.eq.'t91')then
         tc=tk-273.15
         if(tc .le. 500.0)then
            celmod=2.073d11-64.58d6*tc
         else if(tc .le. 600.0)then
            celmod=2.95d11-240.0d6*tc
         end if   
      else if(cmat.eq.'ss316')then
         celmod = (205.91-2.6913E-2*tk-4.1876E-5*tk**2)*1.e9      
      else
         write(*,*)'wrong clad material:',cmat
         stop
      end if
      return
      end

!==================================================================================================
! Calculation of clad failure criteria
!==================================================================================================
      subroutine cfail()

      use globals
      implicit none
      integer j,i
      real*8 da,sigma,sum,tb1(12),sb1(12),tc,el0
      data tb1/297.0,327.0,371.0,450.0,509.0,580.0,642.0,714.0,795.0,874.0,926.0,960.0/
      data sb1/742.7,687.2,629.3,548.3,497.4,458.0,423.3,393.2,360.8,335.3,326.0,316.7/

      if(ifail.ne.0)return

      do j=1,nzz

! CRITERION 1: rupture temperature versus cladding stress

         if(cmat.eq.'ss316'.or.cmat.eq.'t91')then
!           find effective clad stress
            sigma=0.d0  
               sum=0   
            do i=1,nc-1
               tc=tem(i,j)-273.15
               da=rad(nf+i+1,j)**2-rad(nf+i,j)**2      
               sigma=sigma+sig(i,j)*da
               sum=sum+da
            end do
            sigma=sigma/sum
            tc=tcave(j)-273.15
            sigb(j)=(-5.0e-9*tc**4+2.0e-5*tc**3-0.0269*tc**2+13.576*tc-1397.8)*1.0d6
            fcrit(1,j)=sigma/sigb(j)
         else
            fcrit(1,j)=0.d0
         end if

! CRITERION 2: creep rate/time-to-failure criterion
         if(cmat.eq.'ss316')then
            sum=0.d0
            do i=1,nc-1
               da=rad(nf+i+1,j)**2-rad(nf+i,j)**2
               sum=sum+epe(i,j)*da
            end do
            sum=sum/(rco(j)**2-rci(j)**2)
            tc=tcave(j)-273.15
            el0=(2.0e-12*tc**5-1.0e-8*tc**4+2.0E-5*tc**3-0.0173*tc**2+7.9583*tc-1437.6)/100.0
            fcrit(2,j)=sum/el0
         else
            fcrit(2,j)=0.d0
         end if

! CRITERION 3: strain energy density

! CRITERION 4: cladding meting

         if(cmat.eq.'ss316')then
            fcrit(4,j)=tcave(j)/1642.0
         else
            fcrit(4,j)=0.d0
         end if

         do i=1,4
            fcrit(i,j)=dmin1(fcrit(i,j),1.d0)
            if(fcrit(i,j).ge.1.0d0)then
               ifail=i
               return
            end if
         end do
      end do
      return
      end

!==================================================================================================
! Calculates parameters for the cladding stress-strain curve as a function of temperature
!      input arguments
!      t - cladding temperature (k)
!      cmat - clad material index
!      yng - Young's modulus (Pa)
!      rtol - relative tolerance
!      output arguments
!      ak - strength coefficient (pa)
!      an - strain hardening exponent (unitless)
!==================================================================================================
      subroutine ckmn(t,cmat,yng,rtol,ak,an,ssy0,ssb0,el0)
      implicit none
      integer icon,iter
      real*8 t,yng,rtol,ak,an,ssb0,ssy0,el0,g1,csigb,csigy,cuelon
      character*6 cmat

      ssb0=csigb(t,cmat)
      ssy0=csigy(t,cmat)
      el0=cuelon(t,cmat)
      ak=(ssb0-ssy0)/el0
      an=1.d0
      g1=0.0
      icon=0
      iter=0
      do while(icon.eq.0)
         iter=iter+1
         if(iter.gt.1000)stop 'too many ckmn iterations'
         an=dlog(ssy0/ak)/dlog(ssy0/yng)
         g1=ssb0/(el0+ssb0/yng)**an
         if(dabs(ak-g1).lt.rtol*g1)icon=1
         ak=g1
      end do
      return
      end

!==================================================================================================
! Calculation of the clad thermal conductivity
!                unit: w/(m*k)
!     input arguments :
!        tk .... temperature (K)
!==================================================================================================
      real*8 function clamb(tk,cmat)
      implicit none
      real*8 tk, tc
      character*6 cmat

      if(cmat.eq.'zrn')then
         clamb = 15.0636d0*dexp(4.61843d-4*tk)
      else if(cmat.eq.'zry')then
         clamb=7.511d0+tk*(2.088d-2+tk*(-1.450d-5+tk*7.668d-09))
      else if(cmat.eq.'ss823')then
         clamb = 16.d0+0.0046d0*tk
      else if(cmat.eq.'aim1')then
         tc=tk-273.15d0
         clamb=13.95d0+1.163d-2*tc
! PSI+Na_CABRI_PA     
      else if(cmat.eq.'ss316')then
         if(tk.LT.1642.)then
            clamb = 9.248 + 1.571e-2*tk
           else
              clamb = 20.0
           end if
! PSI-Na_CABRI_PA      
      else if(cmat.eq.'t91')then
         tc=tk-273.15d0
         clamb = 23.71 + 0.01718*tc - 1.45E-05*tc**2
      else
         write(*,*)'wrong clad material:',cmat
         stop
      end if
      return
      end

!==================================================================================================
! Calculation of the clad Poisson ratio
!                unitless
!     input arguments :
!        tk .... temperature (k)
!==================================================================================================
      real*8 function cpoir(cmat)
      implicit none
      character*6 cmat

      if(cmat.eq.'aim1'.or.cmat.eq.'ss823')then
         cpoir=0.289d0
      else if(cmat.eq.'ss316'.or.cmat.eq.'t91')then
         cpoir=0.28d0
      else
         write(*,*)'wrong clad material:',cmat
         stop
      end if
      return
      end

!==================================================================================================
! Calculates burst stress for cladding
!==================================================================================================
      real*8 function csigb(tk,cmat)
      implicit none
      real*8 polp,tk,tc,tem823(14),sb823(14)
      character*6 cmat

!    Unirradiated data for EP-823 (Belomitzev):
      data tem823/293.d0,573.d0,623.d0,673.d0,723.d0,773.d0,873.d0,923.d0,973.d0,1053.d0,1173.d0,1273.d0,1373.d0,1473.d0/
      data sb823 /975.d0,764.5d0,745.d0,686.d0,642.0d0,558.5d0,367.5d0,274.5,172.d0,72.d0,47.d0,30.d0,23.d0,13.d0/

      if(cmat.eq.'ss823')then
         csigb=polp(tk,sb823 ,tem823,14)*1.0d6
      else if(cmat.eq.'aim1'.or.cmat.eq.'t91')then
         csigb=1.5957d9 - 4.7253d6*tk + 9.8851d3*tk**2 - 8.8864d0*tk**3 + 2.5538d-3*tk**4
      else if(cmat.eq.'ss316')then
         tc=tk-273.15
         csigb=(-5.0e-9*tc**4+2.0e-5*tc**3-0.0269*tc**2+13.576*tc-1397.8)*1.0d6
      else
         csigb = 0.0
      end if
      return
      end

!==================================================================================================
! Calculate yield stress for cladding
!==================================================================================================
      real*8 function csigy(tk,cmat)
      implicit none
      real*8 polp,tk,tc,tem823(14),sy823(14)
      character*6 cmat

!    Unirradiated data for EP-823 (Belomitzev):
      data tem823/293.d0,573.d0,623.d0,673.d0,723.d0,773.d0,873.d0,923.d0,973.d0,1053.d0,1173.d0,1273.d0,1373.d0,1473.d0/
      data sy823 /960.5d0,749.5d0,720.5d0,676.d0,632.0d0,553.5d0,362.5d0,264.5d0,162.d0,68.d0,42.d0,27.d0,21.d0,12.d0/

      if(cmat.eq.'ss823')then
         csigy=polp(tk,sy823 ,tem823,14)*1.0d6
      else if(cmat.eq.'t91')then
         csigy=1.3109d9 - 3.6916d6*tk + 7.8909d3*tk**2 - 7.3551d0*tk**3 + 2.1966d-3*tk**4
      else if(cmat.eq.'aim1')then
        tc = tk - 273.15
        if(tc.lt.600)then
            csigy=5.555d8 - 0.25d0*tc
        else if(tc.ge.600 .and. tc.le.1000)then
            csigy=4.055d8 - 0.755d0*(tc-600.0d0)
        else
            csigy=3.455d8 - 0.25d0*tc
        end if
      else if(cmat.eq.'ss316')then
         tc=tk-273.15
         csigy=(848.65-0.7665*tc)*1.0d6
      else
         csigy = 0.0
      end if
      return
      end

!==================================================================================================
! Calculation of the clad steady-state volume swelling
!     input arguments :
!        t .... temperature (k)
!        bup... burnup (MWd/kgU)
!        dpa .. displacement per atom
!==================================================================================================
      real*8 function cswel(icswel,t,bup,cmat)
      implicit none
      integer icswel
      real*8 t,dpa,ttt,GO,DDEL,DDELP,dpa_f,Roche_f,bup,fnf
      character*6 cmat
      parameter(Roche_f=0.4)

      if(icswel .eq. 0)then
         cswel=0.0
         return
      end if
!     clad dpa linearly depends on fuel burnup, according to Michael Schikorr July 2012
      dpa=1.4087 * bup
!     temperature in celsium
      ttt=t-273.15d0
      if(cmat.eq.'ss823')then
         cswel=2.8d-6*dexp(-5.5d-4*(ttt-440.d0)**2)*dpa**1.8
      else if(cmat.eq.'t91')then
!        Roche correlation from Michael Schikorr July 2012
         GO=max(0.0017*exp(-((ttt-491.0)/84.4)**2),1.0d-6)
         DDEL=max(111.0*exp(-((ttt-531.0)/108.0)**2),1.0d-6)
         DDELP=max(DDEL*1.36765,1.0d-6)
!        French way to determine dpa
         dpa_f=dpa*1.35
         if(dpa_f.le.DDELP)then
            cswel = GO*(DDELP-DDEL)*(dpa_f/DDELP)**3.72*Roche_f
         else
            cswel = GO*(dpa_f-DDEL)*Roche_f
         end if
      else if(cmat.eq.'aim1')then
!        L. Luzzi, et al., Modeling and Analysis of Nuclear Fuel Pin Behavior for Innovative Lead Cooled FBR,Report RdS/PAR2013/022
         dpa=0.4156*bup ! Based on CAPRIX data bup in (MWd/kgiHM)
!        Fast neutron fluence (n/cm^2) divided by 10^22
         fnf=2*dpa/10.0
         cswel=1.3d-5*dexp(-1*((ttt-490.0)/100.0)**2)*fnf**3.9 
      else
         cswel=0.d0
      end if
      return
      end

!==================================================================================================
! Calculates the clad thermal expansion strain
!==================================================================================================
      real*8 function ctexp(tk,cmat)
      implicit none
      real*8 tk, tc
      character*6 cmat

      if(cmat.eq.'ss316'.or.cmat.eq.'ss823'.or.cmat.eq.'t91')then
         ctexp=-0.2177d-2 + 6.735d-6*tk + 5.12d-9*tk**2 - 2.248d-12*tk**3+3.933d-16*tk**4
      else if(cmat.eq.'aim1')then
         tc = tk - 273.15d0
         ctexp=-3.101d-4 + 1.545d-5*tc + 2.75d-9*tc**2
      else
         write(*,*)'wrong cmat option ',cmat
         stop
      end if
      return
      end

!==================================================================================================
! Calculates ultimate elongation for cladding
!==================================================================================================
      real*8 function cuelon(tk,cmat)
      implicit none
      real*8 polp,tk,tem823(14),el823(14)
      character*6 cmat

!    Unirradiated data for EP-823 (Belomitzev):
      data tem823/293.d0,573.d0,623.d0,673.d0,723.d0,773.d0,873.d0,923.d0,973.d0,1053.d0,1173.d0,1273.d0,1373.d0,1473.d0/
      data el823 /14.d0,11.d0,10.5d0,11.d0,13.5d0,14.d0,21.d0,24.d0,21.5d0,30.d0,30.d0,30.d0,30.d0,30.d0/

      if(cmat.eq.'ss823')then
         cuelon= polp(tk,el823 ,tem823,14)/100.d0
      else if(cmat.eq.'aim1'.or.cmat.eq.'t91')then
         if(tk.le.720.d0)then
            cuelon=-0.58258 + 8.4018d-3*tk - 3.2807d-5*tk**2 + 4.9989d-8*tk**3 - 2.6347d-11*tk**4
         else
            cuelon=-0.85401 + 4.5753d-3*tk - 8.2202d-6*tk**2 + 6.1983d-9*tk**3 - 1.6897d-12*tk**4
         end if           
      else if(cmat.eq.'ss316')then
         cuelon=(2.0e-12*tk**5-1.0e-8*tk**4+2.0E-5*tk**3-0.0173*tk**2+7.9583*tk-1437.6)/100.0
      else
         cuelon = 0.0
      end if
      return
      end

!==================================================================================================
! Calculation of fuel specific heat
!                  unit : j/(kg*k)
!       input arguments :
!          tk .... temperature (k)
!==================================================================================================
      real*8 function fcp (tk,pucont,fmat)
      implicit none
      real*8 tk,pucont,tuo2(18),cpuo2(18),tuc(17),cpuc(17),polp,tet,ucp_,pucp_,tPopov(29),cpuo2Popov(29),cppuo2Popov(29),tliq,tf
      character*3 fmat
      data tuo2  /293.0d0,500.0d0,700.0d0,900.0d0,1100.d0,1300.d0, &
     &            1500.d0,1700.d0,1900.d0,2100.d0,2300.d0,2500.d0, &
     &            2700.d0,2900.d0,3095.d0,3100.d0,3105.d0,5000.d0/
      data cpuo2 /280.0d0,287.0d0,302.0d0,310.0d0,314.0d0,319.0d0, &
     &            320.0d0,328.0d0,340.0d0,364.0d0,390.0d0,426.0d0, &
     &            470.0d0,520.0d0,594.0d0,55467.5d0,503.d0,503.0d0/
      data tuc   /293.0d0,500.0d0,700.0d0,900.0d0,1100.d0,1300.d0, &
     &            1500.d0,1700.d0,1900.d0,2100.d0,2300.d0,2500.d0, &
     &            2700.d0,2775.d0,2780.d0,2785.d0,5000.d0/
      data cpuc  /207.0d0,216.0d0,225.0d0,234.0d0,243.0d0,253.0d0, &
     &            264.0d0,275.0d0,288.0d0,301.0d0,316.0d0,332.0d0, &
     &            350.0d0,357.0d0,54443.d0,357.d0,357.0d0/

! S.G. Popov, J.J. Carbajo, V.K. Ivanov, G.L. Yoder. Thermophysical properties of MOX and UO2 fuels including the effects of irradiation. ORNL/TM-2000/351 http://web.ornl.gov/~webworks/cpr/v823/rpt/109264.pdf
      data tPopov     /300.,400.,500.,600.,700.,800.,900., &
     &                1000.,1100.,1200.,1300.,1400.,1500.,1600., &
     &                1700.,1800.,1900.,2000.,2100.,2200.,2300., &
     &                2400.,2500.,2600.,2700.,2800.,2900.,3000., &
     &                3100./
      data cpuo2Popov /235.51,265.79,282.14,292.21,299.11,304.24,308.32, &
     &                 311.74,314.76,317.59,320.44,323.60,327.41,332.31, &
     &                 338.77,347.30,358.40,372.54,390.12,411.46,436.78, &
     &                 466.21,499.80,537.50,579.17,624.63,673.63,725.88, &
     &                 781.08/
      data cppuo2Popov/203.71,225.79,235.77,241.33,244.92,247.47,249.44, &
     &                 251.04,252.43,253.70,254.96,256.32,257.93,259.93, &
     &                 262.47,265.67,269.57,274.17,279.35,284.96,290.80, &
     &                 296.70,302.47,308.00,313.21,318.08,322.59,326.76, &
     &                 330.64/

      if(fmat.eq.'uo2'.or.fmat.eq.'mox')then
!         fcp = polp (tk,cpuo2,tuo2,18)
! new correlation for OECD XADS beam trip benchmark
!         ucp_=81.825d0+0.78695d0*tk-1.1552d-3*tk**2+9.9037d-7*tk**3-5.198d-10*tk**4+1.5241d-13*tk**5-1.7906d-17*tk**6
!         pucp_=-4.9236d6/tk**2+240.89d0+0.32556d0*tk-3.5398d-4*tk**2+1.512d-7*tk**3-1.9707d-11*tk**4
!         fcp=ucp_*(1.d0-pucont)+pucp_*pucont
!
         ucp_=polp(tk,cpuo2Popov,tPopov,29)
         pucp_=polp(tk,cppuo2Popov,tPopov,29)
         fcp=ucp_*(1.d0-pucont)+pucp_*pucont
!        melting point
         tliq=3120.0-388.1*pucont-30.4*pucont**2
         if(tk.ge.tliq-1.0)then
            fcp=fcp+1.d5*abs(tk-tliq+1.0)
         end if
      else if(fmat.eq.'uc')then
         fcp = polp (tk,cpuc,tuc,17)
      else if(fmat.eq.'un')then
         tk=dmin1(dmax1(tk,293.d0),4000.d0)
         tet = 367.5d0
!        specific heat in j/mol-k
         fcp=51.14d0*(tet/tk)**2*dexp(tet/tk)/(dexp(tet/tk)-1.d0)**2+9.491d-3*tk+2.642d11*dexp(-1.8081d4/tk)/tk**2
!        convert in j/kg-k
         fcp=fcp/0.252d0
      else if(fmat.eq.'bn')then
!        borrowed from TRACE v5
!        1) w.l. kirchner,"reflood heat transfer in a light water reactor,"u.s. nuclear regulatory commission report nureg-0106 (1976).
!        2) n. fujita,et al.,"a prediction of the semiscale blowdown heat transfer test s-02-8 (nrc standard problem five),"electric power research institute report epri np-212 (october 1976).  (cpw)
         tf = tk*9./5. - 459.67
         fcp = 760.59 + tf*(1.7955 + tf*(tf*1.5896e-7 - 8.6704e-4))
      else
         write(*,*)'wrong fuel material:',fmat
         stop
      end if
      return
      end

!==================================================================================================
! Returns fuel material effective creep rate (1/s)
!==================================================================================================
      real*8 function fcreep(i,j)
      use globals
      implicit none
      integer i,j
      real*8 G,fdot,sg,ttt,tden,rrr,ddd,fcreep1,fcreep2,fcreep3,fcreep4

!     universal gas constant in cal/K/mol 
      rrr=1.98722
      fcreep=0.0
      if(fmat.eq.'mox'.and.sigf(i,j).gt.0.0)then
!        LIFE-GCFR model ANL-81-4 http://www.iaea.org/inis/collection/NCLCollectionStore/_Public/14/738/14738447.pdf     

!        grain size (mkm)
         G=20.0

!        fission rate (fis/cm3/s): 1 fiss/s = 200 MeV/s = 200 * 1.0e6 * 1.60214d-19 W
         fdot=(qqv1(j)*1.d-6)/(1.60214d-13*200.d0)

!        stress in psi 
         sg=sigf(i,j)*0.000145
         if(sg.lt.1000.0)return
         if(sg.gt.16000.0)sg=16000.0

         ttt=0.5d0*(tem(i,j)+tem(i+1,j))
         if(ttt.lt.1673.0)return
         if(ttt.gt.1948.0)ttt=1948.0

         tden=11460.d0*pucont+10960.d0*(1.d0-pucont)
         ddd=(rof(i,j)/tden)*100.0
         if(ddd.lt.88.0)ddd=88.0
         if(ddd.gt.97.0)ddd=97.0

         fcreep1=2.23d7*sg/(G**2.0d0)*dexp(-92.5d3/rrr/ttt)
         fcreep1=fcreep1*(1.0d0+2.11d0*(97.0d0-ddd))

         fcreep2=1.0d-3*(sg**4.4d0)*dexp(-136.8d3/rrr/ttt)
         fcreep2=fcreep2*(1.0d0+0.22d0*(97.0d0-ddd))

         fcreep3=1.96d-19*sg*fdot*dexp(-13.7d3/rrr/ttt)

         fcreep4=3.72d-23*sg*fdot

         fcreep=fcreep1+fcreep2+fcreep3+fcreep4
         if(fcreep.gt.5.d-6)fcreep=5.d-6
         fcreep=fcreep/3.6d3 !convert from 1/h to 1/s

!        SAS4A model ? 
!         fcreep=8.0d-2*(sigf(i,j)/1.d6)**3.0d0*dexp(-2.85d4/rrr/ttt)

!        K. Lassmann, A. Moreno, atke, bd.30 (1977) 207-215
!         if(ddd.lt.92.0)ddd=92.0

!         fcreep1 = (9.73d6 + 3.24d-6*fdot)/(ddd-87.7) * (sg/G**2)*exp(-90000./rrr/ttt)
!         fcreep1 = fcreep1/3600.0 !convert from 1/h to 1/s

!         fcreep2 = (1.38d-4 + 4.6d-17*fdot)/(ddd-90.5) * (sg**4.5)*exp(-132000./rrr/ttt)
!         fcreep2 = fcreep2/3600.0 !convert from 1/h to 1/s

!         fcreep3 = 7.0e-23*sg*fdot
!         fcreep3 = fcreep3/3600.0 !convert from 1/h to 1/s

!         fcreep = fcreep1 + fcreep2 + fcreep3

      end if
      return
      end
!==================================================================================================
! Linear fuel densification (MATPRO model)
!  den -- fuel density (kg/m3)
!  cont -- plutonium content (-)
!  tsint -- sintering temperature (K), i.e. temperature at which the fuel was sintered during manufacturing
!  B -- burnup (MWd/kgHM)
!==================================================================================================
      real*8 function fdens(den,cont,tsint,B)
      implicit none
      integer iter
      real*8 den,cont,tsint,B
      real*8 tden,por,D_assy,D,dB0,B0,e1,e2

      if(tsint .EQ. 0.0) then
         fdens = 0.0
         return
      end if

!     theoretical density (kg/m3)
      tden = 11460.d0*cont + 10960.d0*(1.d0-cont)
!     porosity (%)
      por = (1.0 - den/tden) * 100.

!     assymptotic densification (%), i.e. densification at bup --> infinity
      D_assy = -22.2*por/(tsint-1453.0)

!     densification as a function of burnup bup has the following form: D = D_assy + 0.93 * exp(-3.*(B+B0)) + 2.07 * exp(-35.*(B+B0))
!     now we have to find B0 such that D = 0 at B = 0
      B0 = 0.0
      dB0 = 1.e5
      iter = 10
      do while( (dB0.GT.1.0e-5) .AND. (iter.GT.0) )
         e1 = exp(-3.*B0)
         e2 = exp(-35.*B0)
         dB0 = (D_assy + 0.93 * e1 + 2.07 * e2) / (2.79*e1 + 72.45*e2)
         B0 = B0 + dB0
         iter = iter - 1
      end do
      if(iter .LE. 0)then
        write(*,*)'Too many iterations'
        write(*,*)'in fuel densification routine fdens'
        stop
      end if
      D = D_assy + 0.93 * exp(-3.*(B+B0)) + 2.07 * exp(-35.*(B+B0))
      fdens = D * 0.01
      return
      end

!==================================================================================================
! Calculation of the fuel Young's modulus, unit: Pa
!     input arguments :
!        tk .... temperature (k)
!==================================================================================================
      real*8 function felmod(tk,den,cont,fmat)
      implicit none
      real*8 tk,den,cont,tden,por
      character*6 fmat

      if(fmat.eq.'uo2'.or.fmat.eq.'mox')then
         tden=11460.d0*cont+10960.d0*(1.d0-cont)
         por=1.d0-den/tden
!        MATPRO correlation
         felmod=(2.334d11*(1.0-2.752*por)*(1.d0-1.0915d-4*tk))*(1.d0+0.15d0*cont)
      else if(fmat.eq.'bn')then
         felmod=4.7e10
      else
         write(*,*)'wrong fuel material:',fmat
         stop
      end if
      return
      end

!==================================================================================================
! Calculation of fission gas release
!==================================================================================================
      subroutine fgr()

      use globals
      implicit none
      integer j
      real*8 vfz,vcz,fiss,polp,buper,fgr_,tc,bufree,pi,a,fgrR,fgrU,fgg,bu

      if(ifgr .eq. 0)return
      pi=2.0d0*acos(0.d0)
      fgrel=fgrel0
      fggen=fggen0

      do j=1,nzz
         vfz=pi*(rfo0**2-rfi0**2)*dz0(j)
         vcz=pi*(rco0**2-rci0**2)*dz0(j)

!        fuel burnup in percents of heavy metals
         buper=bup(j)*bcnvr
!        number of fissions (1MWd=8.64d10J,1.60214d-13J=1MeV,200MeV=1fiss)
         fiss=(bup(j)-bup0(j))*rof0*vfz*8.64d10/1.60214d-13/200.0
!        26 atoms of FGP are supposed produced per 100 fissions
!        1mol=6.0247d23 atoms
!        fission gas production in moles
         fgg=0.26d0*fiss/avo
         fggen=fggen+fgg

         fgr_ = 0.0
         if(nfgrt.gt.0)then
!           fission gas release as a table of fuel temperature
            fgr_=polp(tfave(j),fgrt,tfgr,nfgrt)

         else if(fmat.eq.'un')then
!           fission gas release
            fgr_=(-0.331d0-0.496d0*buper+0.409d0*buper**2)/100.0

         else if(fmat.eq.'uc')then
            tc = tem(1,j)-273.15d0
!           the data for UC are taken from PREUSSER, T., Modelling of Carbide Fuel Rods, Nucl. Technol. 57 (1982) 343.
            if(tc.le.1455.d0)then
!              bufree - percent of burnup at which fission gas release begins
               bufree=2.d0
            else if (tc.gt.1455.d0.and.tc.le.2325.d0)then
               bufree=-0.0023d0*tc+5.3504d0
            else
               bufree=0.d0
            end if

!           fission gas release in percents
            if(buper.ge.bufree) then
               if(tc.gt.1000.d0.and.tc.le.2070.d0)then
                  fgr_=0.000467d0*tc-0.467d0
               else if(tc.gt.2070.d0) then
                  fgr_=0.741918d0*dlog(0.7675d0*tc)-4.968477d0
               else
                  fgr_ = 0.d0
               end if
            else
               fgr_ = 0.d0
            end if
            fgr_=fgr_*(1.d0-dexp(-1.5d0*(buper-bufree)))

         else if(fmat.eq.'uo2'.or.fmat.eq.'mox')then
!           Waltar and Reynolds model
            bu=bup(j)
            if(bu .gt. 0.0)then
!              calculate release fraction for the restructured fuel
               a = 4.7/bu*(1.0-exp(-bu/5.9))
               if(a .gt. 1.0)then
                  fgrR = 0.0
               else
                  fgrR = 1.0 - a
               end if

!              calculate release fraction for the unrestructured fuel
               if(bu .le. 3.5d0)then
                   a = 100.0
               else    
                   a = 25.6/(bu-3.5d0)*(1.0-exp(-bu/3.5d0)-1.0d0)*exp(-0.0125d0*(ql(j)/1.d3))
                   if(bu .ge. 49.2d0)then
                      a = a*exp(-0.3d0*(bu-49.2d0))
                   end if
               end if    
               if(a .gt. 1.0)then
                  fgrU = 0.0
               else
                  fgrU = 1.0 - a
               end if
!              assume that all fuel is restructured 
               fgr_=fgrR
!              assume that all fuel is unrestructured 
               fgr_=fgrU
            else
               fgr_=0.0
            end if
         end if
         if(fgr_.gt.1.d0)fgr_=1.d0
!        cumulative fission gas release in moles
         fgrel=fgrel+dmax1(fgr_,0.d0)*fgg
      end do

!     fission gas release from fuel in %%
      if(fggen.gt.0.d0)fgrpp=fgrel/fggen*100.d0
      return
      end

!==================================================================================================
! Calculation of the fuel local thermal conductivity, unit : w/(m*k)
!     input arguments :
!        dens - density (kg/m3)
!        bup - burn-up (mwd/kgu)
!        tk - temperature (k)
!        pucon - Pu fraction (-)
!        sto - fuel stoichiometry
!        fmat - fuel material index
!==================================================================================================
      real*8 function flamb(dens,bup,tk,pucon,sto,fmat)
      implicit none
      real*8 dens,bup,tk,pucon,ulm095,ulmb,ak1,ak2,ak3,ak4,b,temp,tden,por,bbb,ac,tc,sto,tf
      character*3 fmat

      if(fmat.eq.'uo2')then
!        thermal conductivity for fresh fuel with 95% density
         ulm095=1.d2/(2.58d-2*tk+3.77d0)+1.1d-4*tk+1.01d-11*tk**3*exp(7.2d-4*tk)

!        thermal conductivity for fresh fuel with the given density
         ulmb=2.158d0*ulm095*dens/(3.291d4-dens)

!        coefficient to account for fuel composition
         ak1=(5.3d-2+2.2d-4*tk)/(5.3d-2+1.71d-3*bup+(2.2d-4-5.33d-8*bup)*tk)

!        coefficient to account for fuel porosity
         if(tk.le.1773.d0)then
            temp=tk-273.15d0
            b=4.4d0-3.2d-3*temp+4.d-7*temp**2
         else
            b=0.5
         end if
         ak2=1.d0-1.d-3*b*bup
         flamb=ulmb*ak1*ak2

!        frapcon-3 model
!        Harding and Martin (1989)
         ulmb=1.d0/(0.0375+2.165d-4*tk)+4.715d9*exp(-16361.d0/tk)/tk**2
         if(bup.gt.0.d0)then
            bbb=bup*0.123
!           effect of burnup
            ak1=1.09/bbb**3.265+0.0643d0*sqrt(tk/bbb)
            ak1=ak1*atan(1.d0/ak1)
!           effect of precipitated fission products
            ak2=1.0d0+0.019d0*bbb/(3.d0-0.019d0*bbb)/(1.d0+exp(12.-tk/100.))
         else
            ak1=1.d0
            ak2=1.d0
         end if
!        effect of fuel porosity
         por=1.d0-dens/10980.d0
         ak3=(1.0-por)/(1.0+0.5d0*por)
!        radiation effect
         ak4=1.d0-0.2d0/(1.0d0+exp((tk-900.0)/80.0d0))
         flamb=ulmb*ak1*ak2*ak3*ak4

      else if(fmat.eq.'mox')then
!        stoichometry
         ac=1.320d0*((2.0d0- sto)+0.0093d0)**0.5d0 - 0.091d0 + 0.0038*bup/9.33d0
         tden=11460.d0*pucon+10960.d0*(1.d0-pucon)
!        fuel porosity
         por=1.d0-dens/tden
!        FBR MOX fuel  thermal conductivity in [W/mK] according to Y. Philipponneau, J. Nuclear Matter., 188 (1992) 194-197
         flamb=(1.0d0/(ac + 2.493d-4*tk) + 88.4d-12*tk**3)*(1.d0-por)/(1.0d0+2.0d0*por)

      else if(fmat.eq.'uc')then
!         flamb=1.d1*tk/(1.0d0+0.5d0*tk) + 1.d-3*tk
!         model based on the GFR report GCFR-DEl-011 25.12.2009 mikityuk https://www.gofastr.org/gcfr/GCFR-DEL-011.pdf/view
          tc=tk-273.15
!         UC thermal conductivity
          if(tk.le.973.0)then
             ak1=21.7 - 3.04e-3*tc + 3.61e-6*tc**2
          else
             ak1=20.2 + 1.48e-3*tc
          end if   
!         PuC thermal conductivity
          ak2=7.45 - 4.04e-3*tc + 1.2e-5*tc**2
          if(ak2.ge.ak1)ak2=ak1
!         mixed (U-Pu)C thermal conductivity
          flamb=(1.0d0-pucon)*ak1+pucon*ak2
!         fuel porosity effect
          tden=13580.0
          por=1.d0-dens/tden
          flamb=flamb*(1.d0-por)/(1.0d0+2.0d0*por)

      else if(fmat.eq.'un')then
!        theoretical density
         tden=14320.d0*(1.d0-pucon)+14230.d0*pucon
!        fuel porosity
         por=1.d0-dens/tden
         flamb=(1.37d0-1.6d0*pucon+1.142d0*pucon**2)*tk**0.41*(1.d0-por)/(1.d0+por)
      else if(fmat.eq.'bn')then
!        borrowed from TRACE v5
!        1) w.l. kirchner,"reflood heat transfer in a light water reactor,"u.s. nuclear regulatory commission report nureg-0106 (1976).
!        2) n. fujita,et al.,"a prediction of the semiscale blowdown heat transfer test s-02-8 (nrc standard problem five),"electric power research institute report epri np-212 (october 1976).  (cpw)
         tf = tk*9./5. - 459.67
         flamb=25.27-1.365e-03*tf
      else
         write(*,*)'wrong fuel material:',fmat
         stop
      end if
      return
      end

!==================================================================================================
! Calculation of the fuel Poisson ratio, unitless
!==================================================================================================
      real*8 function fpoir(fmat)
      implicit none
      character*6 fmat

!     MATPRO model
      if(fmat.eq.'uo2')then
         fpoir=0.316d0
      else if(fmat.eq.'mox')then
         fpoir=0.276d0
      else if(fmat.eq.'bn')then
         fpoir=0.3
      else
         write(*,*)'wrong fuel material:',fmat
         stop
      end if
      return
      end

!==================================================================================================
! Fracture strength of fuel
!==================================================================================================
      real*8 function fsigf(tk,cont,dens,fmat)
      implicit none
      real*8 tk,cont,dens,rmu,tden,por,ttt
      character*3 fmat
      parameter(rmu=8.314d0)

      fsigf = 0.0
      if(fmat.eq.'uo2' .or. fmat .eq. 'mox')then
!        MATPRO model
         tden=11460.d0*cont+10960.d0*(1.d0-cont)
         por=1.d0-dens/tden
         ttt=max(tk,1000.0)
         fsigf=1.7d8*(1-2.62d0*por)**0.5*dexp(-1590.0/rmu/ttt)
      end if
      return
      end

!==================================================================================================
! Volumetric fuel swelling (-) for the burnup step (bup-bup0)
!==================================================================================================
      real*8 function fswel(i,j)
      use globals
      implicit none
      integer i,j
      real*8 buper,buper0,polp,fiss,fiss0,ttt 

      fswel=0.0
      if(ifswel .eq. 0)return
!     if not specified, use MWd/kg-to-%% conversion factor for u-235 fission by thermal neutrons
      if(bcnvr.eq.0.d0)bcnvr=0.123d0
!     fuel burnup in percents of heavy metals
      buper=bup(j)*bcnvr
      buper0=bup0(j)*bcnvr
      ttt=tem(i,j)

      if(fmat.eq.'uo2'.or.fmat.eq.'mox')then
!        number of fissions per unit volume (1MWd=8.64d10J, 1MeV=1.60214d-13J, 1fiss=200MeV)
         fiss=bup(j)*rof0*8.64d10/(1.60214d-13*200.d0)
         fiss0=bup0(j)*rof0*8.64d10/(1.60214d-13*200.d0)
!        MATPRO model for swelling due to solid FPs
         fswel=fswel+2.5e-29*(fiss-fiss0)
!        MATPRO model for swelling due to gaseous FPs
         if(ttt .lt. 2800.0)then
            fswel=fswel+8.8d-56*(fiss-fiss0)*(2800.0-ttt)**11.73*exp(-0.0162*(2800.0-ttt)) * exp(-8.0e-27*fiss)
         end if
      end if
!     table of additional fuel swelling (%) vs fuel burnup (MWd/kg)
      if(nfswb.gt.0)then
         fswel=polp(bup(j),fswb,bfsw,nfswb)*(bup(j)-bup0(j))/100.d0
      end if
      return
      end

!==================================================================================================
! Calculation of fuel thermal expansion, unit : m/m
!       input arguments :
!          tk .... temperature (k)
!==================================================================================================
      real*8 function ftexp(tk,cont,fmat)
      implicit none
      real*8 tk,cont,tc,ftexpUO2,ftexpPuO2
      character*3 fmat

      if(fmat.eq.'uo2'.or.fmat.eq.'mox')then
!        MATPRO model
         ftexpUO2  = 1.0d-5*tk - 3.0d-3 + 4.d-2*dexp(-5000./tk)
         ftexpPuO2 = 0.9d-5*tk - 2.7d-3 + 7.d-2*dexp(-5072./tk)
         ftexp = ftexpUO2*(1.0-cont) + ftexpPuO2*cont 
      else if(fmat.eq.'uc')then
         ftexp=8.51d-6*tk+1.70d-9*tk**2-2.639d-3
!        model based on the GFR report GCFR-DEl-011 25.12.2009 mikityuk
         tc = tk - 273.15 
         ftexp=1.004d-5*tc + 1.17d-9*tc**2
      else if(fmat.eq.'un')then
         ftexp=(-0.224+7.41d-4*tk+0.82d-7*tk**2)/100.d0
      else if(fmat.eq.'bn')then
         ftexp=12.0e-6*(tk-293.0)
      else
         write(*,*)'wrong fuel material:',fmat
         stop
      end if
      return
      end

!==================================================================================================
! Open and closed fuel/clad gap conductance
!  j - axial slice
!==================================================================================================
      real*8 function gaphtc(j)

      use globals
      implicit none
      integer j,i,jj
      real*8 gtemp,fe,sbc,emissf,emissc,cmhard,r,conf,conc,fkm,flamb,clamb,ahe,pgas,gask_(4), &
     &       xmol_(4),sumx,rmol,h0_,tf,tc,hgap(3),pkon,pi,mw(4),gap_
!             he        ar        kr      xe
      data mw/4.0026d0, 39.948d0, 83.8d0, 131.30d0/

      pi=2.0d0*acos(0.d0)
!     Stefan-Boltzmann constant (w/m**2-k**4)
      data sbc/5.6697d-8/

      tf=tem(nf,j)
      tc=tem(nf+1,j)

!     average gap temperature
      gtemp=(tf+tc)/2.d0
!     helium thermal conductivity
      gask_(1)=2.639d-3*gtemp**0.7085
!     argon thermal conductivity
      gask_(2)=2.986d-4*gtemp**0.7224
!     kripton thermal conductivity
      gask_(3)=8.247d-5*gtemp**0.8363
!     xenon thermal conductivity
      gask_(4)=4.351d-5*gtemp**0.8618
!     number of moles of technological gas
      rmol=gpres0*vfree0/293.d0/rmu
      if(gmat.eq.'he')then
         xmol_(1)=rmol
         xmol_(2)=0.d0
      else if(gmat.eq.'ar')then
         xmol_(1)=0.d0
         xmol_(2)=rmol
      else
         write(*,*)'wrong gas gap material ',gmat
         write(*,*)'only he and ar allowed'
         stop
      end if
!     number of moles of released fission gases (0.8846xe,0.0769kr,0.0385he)
      xmol_(1)=xmol_(1)+0.0385d0*fgrel
      xmol_(3)=0.0769d0*fgrel
      xmol_(4)=0.8846d0*fgrel
      sumx=0.d0
      do i=1,4
         sumx=sumx+xmol_(i)
         xmol(i)=xmol_(i)
      end do

      gask(j)=1.0d0
      do i=1,4
         gask(j)=gask(j)*gask_(i)**(xmol_(i)/sumx)
      end do

!      gask(j)=0.0d0
!      sumx=0.d0
!      do i=1,4
!         sumx=sumx+xmol_(i)
!         do jj=1,4
!            if(i.ne.jj)then
!               term=mw(i)/mw(jj)
!               sumx=sumx+xmol_(jj)*(1.0d0+dsqrt(gask_(i)*dsqrt(term)/gask_(jj)))**2/dsqrt(8.0d0*(1.0d0+term))
!            end if
!         end do
!         if(sumx.gt.0.d0)gask(j)=gask(j)+xmol_(i)*gask_(i)/sumx
!      end do
      
!     open gap conductance
      ahe=0.425d0-2.3d-4*dmin1(gtemp,1000.d0)
      pgas=dmax1(gpres,0.1d6)
!     accomodation distance
      ajump(j)=0.024688d0*gask(j)*dsqrt(gtemp)*2.d0/(pgas*ahe)

!     using FRAPCON model reduce the gap width due to fuel cracking while calculating the gap conductance 
      reloc(j)=0.25 + dmax1( 5.0, dmin1(10.,ql(j)/4000.) )*( dmin1( 1.0, bup(j)/5.) + 1.0 )/100.

!     effective roughness
      r=sqrt(rufc**2+ruff**2)
      gap_=max(gap(j), r)
 
      gapth(j)=max(gap_*(1.0 - reloc(j)), 0.0) + ajump(j)
      hgap(1)=gask(j)/gapth(j)

!     radiation conductance
!     set fuel and clad surface emissivities to reasonable value.
!     according to KIT recommendation (Struwe+Schikorr) 15.06.2011 
      emissf = 0.8d0
      emissc = 0.9d0

      fe=1.0/(1.0/emissf+(rfo(j)/rci(j))*(1.0/emissc-1.0))
      hgap(2)=sbc*fe*(tf**2+tc**2)*(tf+tc)

      if(flag(j).eq.'open')then
         hgap(3)=0.d0
      else
!        closed gap conductance
         tc=tcave(j)
!        Meyer hardness for SS (from SAS4A hgap routine)
         if(tc.le.893.9203)then
            cmhard=5.961d9*tc**(-0.206d0)
         else
            cmhard=2.750d28*tc**(-6.530d0)
         end if
         if(cmhard.lt.1.e5)cmhard=1.e5

         pkon=pfc(j)-gpres

!        effective roughness
         r=sqrt(rufc**2+ruff**2)

!        effective conductivity
         conf=flamb(rof0,bup(j),tf,pucont,sto0,fmat)
         conc=clamb(tc,cmat)
         fkm=2.d0*conf*conc/(conf+conc)

!        closed gap conductance
         hgap(3)=33.3d0*fkm*pkon/(cmhard*sqrt(r))
      end if
      gaphtc=hgap(1)+hgap(2)+hgap(3)
      hgapi(1,j)=hgap(1)
      hgapi(2,j)=hgap(2)
      hgapi(3,j)=hgap(3)
      return
      end

!==================================================================================================
! Solves the set of N linear equations : A * x = B. Here A is input,
! not as the matrix A but rather as its LU decomposition, determined
! by the routine LUDCMP. INDX is input as the permutation vector
! returned by LUDCMP. B is input as right-hand side vector, and returns
! with the solution vector x. A, N, NP and INDX are not modified by
! this routine and can be left in place for successive calls with
! different right-hand side vector B. This routines takes into account
! the possibility that B will begin with many zero elements, so it is
! efficient for use in matrix inversion.
!==================================================================================================
      SUBROUTINE LUBKSB(A, N, NP, INDX, B)

      IMPLICIT NONE
      INTEGER N, NP, INDX(N), II, I, LL, J
      REAL*8 A(NP,NP), B(N), SUM

      II = 0

      DO  I = 1, N
          LL = INDX(I)
          SUM = B(LL)
          B(LL) = B(I)

          IF(II .NE. 0) THEN
             DO  J = II, I-1
                 SUM = SUM - A(I,J)*B(J)
             END DO

          ELSE IF(SUM .NE. 0.0D0) THEN

!            A NON-ZERO ELEMENT WAS ENCOUNTERED, SO FROM NOW ON WE WILL
!            HAVE TO DO THE SUMS IN THE LOOP ABOVE
             II = I

          ENDIF
          B(I) = SUM
      END DO

!     BACK SUBSTITUTION
      DO  I = N, 1, -1
          SUM = B(I)
          IF(I .LT. N) THEN
             DO  J = I+1, N
                 SUM = SUM - A(I,J)*B(J)
             END DO
          END IF

!         STORE A COMPONENT OF THE SOLUTION VECTOR X
          B(I) = SUM/A(I,I)
      END DO
      RETURN
      END

!==================================================================================================
! Given  an N*N matrix A, with physical dimension NP, this routine
! replaces it by the LU decomposition of a rowwise permutation of
! itself. A and N are input. A is output; INDX is an output vector
! which record the row permutation effected by the partial pivoting;
! D is output as +-1 depending on whether the number of row interchanges
! was even or odd, respectively. This routine is used in combination
! with LUBKSB to solve linear equations or invert matrix.
! nmax - largest expected N; tiny - small number; vv stores the implicit
! scalling of each row
!==================================================================================================
      SUBROUTINE LUDCMP(A, N, NP, INDX, D)

      IMPLICIT NONE
      INTEGER NMAX
      REAL*8 TINY
      PARAMETER (NMAX=3000,TINY=1.0D-20)
      INTEGER N, NP, INDX(N), I, J, K, IMAX
      REAL*8 A(NP,NP), VV(NMAX), D, AAMAX, SUM, DUM

      IF(N.GT.NMAX) THEN
         WRITE(*,*)'LU ERROR: NUMBER OF EQUATIONS IS OUT OF RANGE',N
         STOP
      END IF

!     NO ROW INTERCHANGES YET
      D = 1.0D0

!     CHECK MATRIX SINGULARITY
      DO  I = 1, N
          K = 0
          DO  J = 1, N
              IF(A(J,I) .NE. 0.0D0) K = 1
          END DO
          IF(K .EQ. 0) THEN
           WRITE(*,*)'LU ERROR: SINGULAR MATRIX'
           WRITE(*,*)'ALL ZEROS IN COLUMN ', I
           STOP
          END IF
      END DO

!     LOOP OVER ROWS TO GET THE IMPLICIT SCALING INFORMATION
      DO  I = 1, N

          AAMAX = 0.0D0
          DO  J = 1, N
              IF(DABS(A(I,J)) .GT. AAMAX) AAMAX = DABS(A(I,J))
          END DO

!         NO NONZERO LARGEST ELEMENT
          IF(AAMAX .EQ. 0.0D0) THEN
             WRITE(*,*)'LU ERROR: SINGULAR MATRIX'
            WRITE(*,*)'ALL ZEROS IN ROW ', I
           STOP
          END IF

!         SAVE THE SCALING
          VV(I) = 1.0D0/AAMAX
      END DO

!     THIS IS THE LOOP OVER COLUMNS OF CROUT'S METHOD
      DO  J = 1, N

          IF(J .GT. 1) THEN
             DO  I = 1, J-1
                 SUM = A(I,J)
                 IF(I .GT. 1) THEN
                    DO  K = 1, I-1
                        SUM = SUM - A(I,K)*A(K,J)
                    END DO
                    A(I,J) = SUM
                 END IF
             END DO
          ENDIF

!         INITIALIZE FOR THE SEARCH FOR LARGEST PIVOT ELEMENT
          AAMAX = 0.0D0

          DO  I = J, N
              SUM = A(I,J)
              IF(J .GT. 1) THEN
                 DO  K = 1, J-1
                     SUM = SUM - A(I,K)*A(K,J)
                 END DO
                 A(I,J) = SUM
              END IF

!             FIGURE OF MERIT FOR THE PIVOT
              DUM = VV(I)*DABS(SUM)

!             IS IT BETTER THEN THE BEST SO FAR?
              IF(DUM .GE. AAMAX) THEN
                 IMAX = I
                 AAMAX = DUM
              END IF
          END DO

!         DO WE NEED TO INTERCHANGE ROWS?
          IF (J .NE. IMAX) THEN

              DO  K = 1, N
!                 YES DO SO...
                  DUM = A(IMAX,K)
                  A(IMAX,K) = A(J,K)
                  A(J,K) = DUM
              END DO

!             ...AND CHANGE THE PARITY OF D
              D = -D

!              VV(IMAX) = VV(J)
          END IF

          INDX(J) = IMAX
          IF(J .NE. N) THEN
             IF(A(J,J) .EQ. 0.0D0) A(J,J) = TINY

!            NOW, FINALLY, DIVIDE BY THE PIVOT ELEMENT
             DUM = 1.0D0/A(J,J)

             DO  I = J+1, N
                 A(I,J) = A(I,J)*DUM
             END DO
          END IF

!     GO BACK FOR THE NEXT COLUMN IN THE REDUCTION
      END DO

      IF(A(N,N) .EQ. 0.0D0) A(N,N) = TINY
      RETURN
      END

!==================================================================================================
! Write outfrd file
!==================================================================================================
      subroutine outfrd(time,nnn)

      use globals
      implicit none
      integer nnn
      integer i,j
      real*8 time
      real*8 flamb
      character*12 s
      character*20 title

      write(s,'(i12.12)') nnn
      open(700,file='outfrd' // s)
      write(700,*)
      title = "time (s)            "
      write(700,2002)title,time

      write(700,*)
      title = "time (d)            "
      write(700,2002)title,time/86400.

      write(700,'(500a)')('-',j=1,100)
      title = "rof0 (kg/m3)        "
      write(700,2001)title,rof0

      title = "fggen (mol)         "
      write(700,2001)title,fggen

      write(700,*)
      title = "fgrel (mol)         "
      write(700,2001)title,fgrel

      write(700,*)
      title = "gpres (Pa)          "
      write(700,2001)title,gpres


      write(700,'(500a)')('-',j=1,100)
      title = "z (m)               "
      write(700,2001)title,(dz0(j)*(real(j)-0.5),j=1,nzz)

      write(700,*)
      title = "tfin (C)            "
      write(700,2001)title,(tem(1,j)-273.15,j=1,nzz)

      write(700,*)
      title = "tfout (C)           "
      write(700,2001)title,(tem(nf,j)-273.15,j=1,nzz)

      write(700,*)
      title = "tcin (C)            "
      write(700,2001)title,(tem(nf+1,j)-273.15,j=1,nzz)

      write(700,*)
      title = "tcout (C)           "
      write(700,2001)title,(tem(nf+nc,j)-273.15,j=1,nzz)

      write(700,*)
      title = "tfave (C)           "
      write(700,2001)title,(tfave(j)-273.15,j=1,nzz)

      write(700,*)
      title = "tcave (C)           "
      write(700,2001)title,(tcave(j)-273.15,j=1,nzz)

      write(700,*)
      title = "tcool (C)           "
      write(700,2001)title,(tcool(j)-273.15,j=1,nzz)

      write(700,*)
      title = "rfi (m)             "
      write(700,2001)title,(rfi(j),j=1,nzz)

      write(700,*)
      title = "rfo (m)             "
      write(700,2001)title,(rfo(j),j=1,nzz)

      write(700,*)
      title = "rci (m)             "
      write(700,2001)title,(rci(j),j=1,nzz)

      write(700,*)
      title = "rco (m)             "
      write(700,2001)title,(rco(j),j=1,nzz)

      write(700,*)
      title = "dzf (m)             "
      write(700,2001)title,(dzf(j),j=1,nzz)

      write(700,*)
      title = "dzc (m)             "
      write(700,2001)title,(dzc(j),j=1,nzz)

      write(700,*)
      title = "reloc (-)           "
      write(700,2001)title,(reloc(j),j=1,nzz)

      write(700,*)
      title = "qs (W/m2)           "
      write(700,2001)title,(qs(j),j=1,nzz)

      write(700,*)
      title = "qv (W/m3)           "
      write(700,2001)title,(qqv1(j),j=1,nzz)

      write(700,*)
      title = "ql (W/cm)           "
      write(700,2001)title,(ql(j)/100.,j=1,nzz)

      write(700,*)
      title = "ql2 (W/cm)          "
      write(700,2001)title,(ql(j)/100.*dz0(j)/dzf(j),j=1,nzz)

      write(700,*)
      title = "bup (MWd/kg)        "
      write(700,2001)title,(bup(j),j=1,nzz)

      write(700,*)
      title = "kfuel (W/mK)        "
      write(700,2001)title,(flamb(rof0,bup(j),tfave(j),pucont,sto0,fmat),j=1,nzz)

      write(700,*)
      title = "gap (m)             "
      write(700,2001)title,(gap(j),j=1,nzz)

      write(700,*)
      title = "gapth (m)           "
      write(700,2001)title,(gapth(j),j=1,nzz)

      write(700,*)
      title = "hgap (W/m2K)        "
      write(700,2001)title,(hgapt(j),j=1,nzz)

      write(700,*)
      title = "hgap1 (W/m2K)       "
      write(700,2001)title,(hgapi(1,j),j=1,nzz)

      write(700,*)
      title = "hgap2 (W/m2K)       "
      write(700,2001)title,(hgapi(2,j),j=1,nzz)

      write(700,*)
      title = "hgap3 (W/m2K)       "
      write(700,2001)title,(hgapi(3,j),j=1,nzz)

      write(700,*)
      title = "ajump (m)           "
      write(700,2001)title,(ajump(j),j=1,nzz)

      write(700,*)
      title = "gask (W/m-K)        "
      write(700,2001)title,(gask(j),j=1,nzz)

      write(700,*)
      title = "pfc (Pa)            "
      write(700,2001)title,(pfc(j),j=1,nzz)

      write(700,*)
      title = "htc (W/m-K)         "
      write(700,2001)title,(htc(j),j=1,nzz)

      write(700,*)
      title = "vcool (m/s)         "
      write(700,2001)title,(vcool(j),j=1,nzz)

      write(700,'(500a)')('-',j=1,100)
      title = "r (m)               "
      write(700,2001)title,(rad0(i),i=1,nf+nc)
      
      write(700,*)
      title = "temperature (C)     "
      do j=1,nzz
         write(700,2000)title,j,(tem(i,j)-273.15,i=1,nf+nc)
      end do

      write(700,*)
      title = "sig r (Pa)          "
      do j=1,nzz
         write(700,2000)title,j,(sigfr(i,j),i=1,nf),(sigr(i,j),i=1,nc)
      end do

      write(700,*)
      title = "sig h (Pa)          "
      do j=1,nzz
         write(700,2000)title,j,(sigfh(i,j),i=1,nf-1),(sigh(i,j),i=1,nc-1)
      end do

      write(700,*)
      title = "sig z (Pa)          "
      do j=1,nzz
         write(700,2000)title,j,(sigfz(i,j),i=1,nf-1),(sigh(i,j),i=1,nc-1)
      end do

      write(700,*)
      title = "eps r total (%)     "
      do j=1,nzz
         write(700,2000)title,j,(efr(i,j)*100.,i=1,nf-1),(er(i,j)*100.,i=1,nc-1)
      end do

      write(700,*)
      title = "eps h total (%)     "
      do j=1,nzz
         write(700,2000)title,j,(efh(i,j)*100.,i=1,nf),(eh(i,j)*100.,i=1,nc)
      end do

      write(700,*)
      title = "eps z total (%)    "
      do j=1,nzz
         write(700,2000)title,j,efz(j)*100.,ez(j)*100.
      end do

      write(700,*)
      title = "eps r plastic (%)   "
      do j=1,nzz
         write(700,2000)title,j,(erpf(i,j)*100.0,i=1,nf-1),(erp(i,j)*100.0,i=1,nc-1)
      end do

      write(700,*)
      title = "eps h plastic (%)   "
      do j=1,nzz
         write(700,2000)title,j,(ehpf(i,j)*100.0,i=1,nf-1),(ehp(i,j)*100.0,i=1,nc-1)
      end do

      write(700,*)
      title = "eps z plastic (%)   "
      do j=1,nzz
         write(700,2000)title,j,(ezpf(i,j)*100.0,i=1,nf-1),(ezp(i,j)*100.0,i=1,nc-1)
      end do

      write(700,*)
      title = "eps r creep (%)     "
      do j=1,nzz
         write(700,2000)title,j,(ercf(i,j)*100.0,i=1,nf-1),(erc(i,j)*100.0,i=1,nc-1)
      end do

      write(700,*)
      title = "eps h creep (%)     "
      do j=1,nzz
         write(700,2000)title,j,(ehcf(i,j)*100.0,i=1,nf-1),(ehc(i,j)*100.0,i=1,nc-1)
      end do

      write(700,*)
      title = "eps z creep (%)     "
      do j=1,nzz
         write(700,2000)title,j,(ezcf(i,j)*100.0,i=1,nf-1),(ezc(i,j)*100.0,i=1,nc-1)
      end do

      write(700,*)
      title = "eps thermal lin (%) "
      do j=1,nzz
         write(700,2000)title,j,(eft(i,j)*100.0,i=1,nf-1),(et(i,j)*100.0,i=1,nc-1)
      end do

      write(700,*)
      title = "eps swell lin (%)   "
      do j=1,nzz
         write(700,2000)title,j,(efs(i,j)*100.0,i=1,nf-1),(ecs(i,j)*100.0,i=1,nc-1)
      end do

      write(700,*)
      title = "eps densific    (%) "
      do j=1,nzz
         write(700,2000)title,j,(efd(i,j)*100.0,i=1,nf-1)
      end do

      write(700,*)
      title = "density (kg/m3)     "
      do j=1,nzz
         write(700,2000)title,j,(rof(i,j),i=1,nf-1),(roc(j),i=1,nc-1)
      end do
      close(700)
2000  format(a20,1x,"iz:",i2,500(1x,1pe12.5))
2001  format(a20,1x,"   ",2x,500(1x,1pe12.5))
2002  format(a20,1x,"   ",2x,500(1x,1pe12.5))

      return
      end

!==================================================================================================
! polp is a linear interpolation function whose value is equal to
! interpolated value.
! tc-is the value of the independent variable
! crf-dependent variable array
! crft-independent variable array
! n-number of points in variable arrays
!==================================================================================================
      function polp (tc,crf,crft,n)

      implicit real*8 (a-h,o-z)
      dimension crf(n),crft(n)
      i=1
      if(tc.ge.crft(1))then
         i=n
         if(tc.lt.crft(n))then
            do j=1,n
               i=j
               if(tc-crft(j) < 0.0)then
                  polp=crf(i-1)+(crf(i)-crf(i-1))*(tc-crft(i-1))/(crft(i)-crft(i-1))
                  return
               else if(tc-crft(j) == 0.0)then
                  polp=crf(i)
                  return
               end if
            end do
         end if
      end if
      polp=crf(i)
      return
      end

!==================================================================================================
! Internal gas pressure calculation
!==================================================================================================
      subroutine precal()

      use globals
      implicit none
      integer j
      real*8 v_t,tgap,thol,vpl,tfree0,pi,h0

      pi=2.0d0*acos(0.d0)

      if(ifail.ne.0)then
         gpres=press(1)
         return
      end if

!     initial fuel stack height
      h0=0.d0
      do j=1,nzz
         h0=h0+dz0(j)
      end do

!     as-fabricated temperature
      tfree0=293.d0

!     current free volume and volume-to-temperature ratio
      vfree=0.d0
      v_t=0.d0
      do j=1,nzz
!       gap volume + hole volume + open porosity volume
        vfree=vfree+pi*dzf(j)*(rci(j)**2-rfo(j)**2+rfi(j)**2)
!       gap and hole temperature
        tgap=0.d0
        thol=0.d0
        tgap=(tem(nf,j)+tem(nf+1,j))/2.d0
        thol=tem(1,j)
!       gap vol/tem + hole vol/tem + fuel open porosity/ave.tem
        v_t=v_t+pi*dzf(j)*((rci(j)**2-rfo(j)**2)/tgap+rfi(j)**2/thol)
      end do

!     plenum volume and temperature
      vpl=vgp*(1.d0-dlfuel)

!     add current plenum volume
      vfree=vfree+vpl
     
!     add current plenum volume/plenum tem.
      v_t=v_t+vpl/tple

!     free volume temperature
      tfree=vfree/v_t
!     free volume gas pressure
      gpres=(gpres0*vfree0/tfree0+fgrel*rmu)*tfree/vfree

      return
      end

!==================================================================================================
! Reads input deck and assigns initial values
!==================================================================================================
      subroutine rdinp(inp)

      use globals
      implicit none
      integer maxrow
      parameter (maxrow=500)
      integer numlin,i1,i2,i,j,nsk,idvkq,n2012,n300001
      real*8, intent(inout) :: inp(20)
      real*8 v,r1,press_,vtot,pi,h0
      character*1 c
      character*20 w2

      pi=2.0d0*acos(0.d0)

!     table dimensions
      nqv=0
      ntco=0
      nsk=0
      idvkq=0
      n2012=0
      n300001=0
      nfswb=0      
      ncswb=0      
      nfsto=0      
      nzz=0
! default values for options      
      jaxial=1
      ifgr=0
      ifcreep=0
      ifswel=0
      iccreep=0
      icswel=0
      icplas=0
      ifplas=0
      tsint=0.0

      open(700, file='fred.inp')
!      first symbol of the card
       read(700,'(a1)',err=1000) c
!      input deck line number
       numlin = 1

!      $ - terminal card
       do while(c .ne. '$')
!         * - comment card
          if(c .ne. '*') then
             backspace 700
!            first integer of the card
             read(700,*,err=1000)i1

!            time beginning card
             if(i1.eq.000000) then
                backspace 700
!               rtol - relative tolerance parameter (scalar)
!               atol - absolute tolerance parameter (scalar)
                read(700,*,err=1000)i1,rtol,atol

!            options card
             else if(i1.eq.000002) then
                backspace 700
                read(700,*,err=1000)i1,w2
                
                if(w2 .eq.'OUTPUT_AXIAL_LAYER')then
                   backspace 700
                   read(700,*,err=1000)i1,w2,jaxial
                   
                else if(w2 .eq. 'FGR')then
                   ifgr=1
                   
                else if(w2 .eq. 'FUEL_CREEP')then
                   ifcreep=1
                   
                else if(w2 .eq. 'FUEL_SWEL')then
                   ifswel=1
                   
                else if(w2 .eq. 'CLAD_CREEP')then
                   iccreep=1
                   
                else if(w2 .eq. 'CLAD_SWEL')then
                   icswel=1
                   
                else if(w2 .eq. 'CLAD_PLAS')then
                   icplas=1
                   
                else if(w2 .eq. 'FUEL_PLAS')then
                   ifplas=1
                   
                else if(w2 .eq. 'FUEL_DENS')then
                   backspace 700
                   read(700,*,err=1000)i1,w2,tsint
                end if   

!            initial temperature card
             else if(i1.eq.000006) then

                backspace 700
                read(700,*,err=1000)i1,tem0

!            fuel rod axial division card
             else if(i1.eq.000003)then
                backspace 700
!               dz0 - axial slice height
                
                read(700,*,err=1000)i1,r1,i2
                if(i2.gt.maxz)stop 'card 000003:increase maxz'
                nzz=i2
                do j=1,i2
                   dz0(j)=r1
                end do

!            fuel card
             else if(i1.eq.100001)then
                backspace 700

!               fmat  - fuel material (UO2, MOX, UC, UN, none)
!               rof - fuel density (kg/m3)
!               pucont - plutonium content
!               rfi0  - inner fuel radius (m)
!               rfo0  - outer fuel radius (m)
!               ruff - arithmetic mean roughness height of fuel (m)
!               sto0 - initial fuel stoichiometry
!               nf    - number of radial fuel nodes

               read(700,*,err=1000) i1,fmat,rof0,pucont,rfi0,rfo0,ruff,sto0,nf
               if(fmat.eq.'none')nf=0
               do j=1,nzz
               do i=1,nf-1
                  rof(i,j)=rof0
               end do
               end do

 !           gap card
             else if(i1.eq.100002)then
                backspace 700
!               gmat - gap material (HE,AR,PB)
!               gap  - gap width (m)
!               gpres0 - initial inner gas pressure (Pa)
!               vgp - gas plenum volume
                read(700,*,err=1000)i1,gmat,gap(1),gpres0,vgp
                do j=2,nzz
                   gap(j)=gap(1)
                   gap0(j)=gap(1)
                end do
                tple=293.15

!            clad card
             else if(i1.eq.100003)then
                backspace 700
!               cmat - clad material (ZRN, ZRY, SS)
!               rco0  - outer clad radius (m)
!               roc - clad density (kg/m3)
!               rufc - arithmetic mean roughness height of cladding (m)
                read(700,*,err=1000)i1,cmat,rco0,roc0,rufc,nc
                if(cmat.eq.'none')nc=0
                do j=1,nzz
                   roc(j)=roc0
                end do

!            MWd/kg-to% conversion factor
             else if(i1.eq.100006)then
                backspace 700
                read(700,*,err=1000)i1,bcnvr

!            clad swelling vs fuel burnup card
             else if(i1.eq.100007) then
                ncswb=ncswb+1

                if(ncswb.gt.maxtab) stop 'increase maxtab'
                backspace 700
!               ncswb - number of clad swelling vs fuel burnup entries
!               bcsw - burnup array corresponding to fcwb (%)
!               cswb - array of clad swelling vs fuel burnup (%)
                read(700,*,err=1000)i1,bcsw(ncswb),cswb(ncswb)

!            clad outer temperature card
             else if(i1.eq.200000)then
                ntco=ntco+1
                if(ntco.gt.maxtab) stop 'increase maxtab'
                backspace 700
!               ntco  - number of clad outer temperature table entries
!               tco - array of fuel rod clad outer temperature vs time ttco (K)
!               ttco - time array corresponding to tco (s)
                read(700,*,err=1000)i1,ttco(ntco),(tco(i,ntco),i=1,nzz)
                
!            table 'power distribution mutiplier vs time'
             else if(i1.eq.300001)then
                n300001=1
                backspace 700
!               qvmlt - power distribution multiplier
                nqvmlt=nqvmlt+1
                if(nqvmlt.gt.maxtab)then
                   write(*,*)'input line:',numlin
                   stop 'increase maxtab'
                end if
                read(700,*,err=1000)i1,tqvmlt(nqvmlt),qvmlt(nqvmlt)

!            power card
             else if(i1.eq.300000)then
                nqv=nqv+1
                if(nqv.gt.maxtab) stop 'increase maxtab'
                backspace 700
!               nqv  - number of power table entries
!               qv - array of fuel rod power density vs time tqv (kW/m)
!               tqv - time array corresponding to qv (s)
                read(700,*,err=1000)i1,tqv(nqv),(qv(i,nqv),i=1,nzz)

!            fission gas release vs temperature table
             else if(i1.eq.400003) then
                nfgrt=nfgrt+1
                if(nfgrt.gt.maxtab) stop 'increase maxtab'
                backspace 700
!               nfgrt - number of FGR vs. temperature table entries
!               tfgr - temperature array corresponding to fgrt (K)
!               fgrt - array of FGR vs. temperature (1/1)
                read(700,*,err=1000)i1,tfgr(nfgrt),fgrt(nfgrt)
                if(tfgr(nfgrt).lt.273.15d0)then
                   write(*,*)'wrong input temp. in table: line',numlin
                   stop
                end if
                if(fgrt(nfgrt).gt.1.d0)fgrt(nfgrt)=1.d0
                if(fgrt(nfgrt).lt.0.d0)fgrt(nfgrt)=0.d0

!            fuel swelling vs fuel burnup table
             else if(i1.eq.400004) then
                nfswb=nfswb+1
                if(nfswb.gt.maxtab) stop 'increase maxtab'
                backspace 700
!               nfswb - number of fuel swelling vs fuel burnup entries
!               bfsw - burnup array corresponding to fswb (%)
!               fswb - array of fuel swelling vs fuel burnup (%)
                read(700,*,err=1000)i1,bfsw(nfswb),fswb(nfswb)

!            fuel stoichiometry vs fuel burnup table
             else if(i1.eq.400005) then
                nfsto=nfsto+1
                if(nfsto.gt.maxtab) stop 'increase maxtab'
                backspace 700
!               nfsto - number of fuel stoichiometry vs fuel burnup entries
!               bsto - burnup array corresponding to stob
!               stob - array of fuel stoichiometry vs fuel burnup
                read(700,*,err=1000)i1,bsto(nfsto),stob(nfsto)

!            coolant card
             else if(i1.eq.500000)then
                backspace 700
!               ctype - coolant type
!               press - coolant pressure (Pa)
!               dhyd - hydraulic diameter of the channel (m)
!               xarea - cross-sectional area of the channel (m2)
!               flowr - flowrate (kg/s)
!               tcoolin - inlet coolant temperature (k) 
                read(700,*,err=1000)i1,ctype,dhyd,xarea,press_,flowr,tcoolin
                do j=1,nzz
                   press(j)=press_
                end do

!            base-irradiation time card
             else if(i1.eq.600000)then
                backspace 700

!               tend - end time
!               dtout - output time step
                read(700,*,err=1000)i1,tend,dtout

             else
                write(*,*)'input error: wrong card number',i1
                stop
             end if
          end if
!         first symbol of the next string
          read(700,'(a1)',err=1000)c
          numlin=numlin+1
       end do
      close(700)

!     some preparations...

!     set initial geometry
      do j=1,nzz
         rci0=rfo0+gap(j)
         rfi(j)=rfi0
         rfo(j)=rfo0
         rco(j)=rco0
         rci(j)=rci0
         dzf(j)=dz0(j)
         dzc(j)=dz0(j)
         flag(j)='open'
      end do

!     radii of fuel nodes
      drf0=(rfo0-rfi0)/dfloat(nf-1)
      do i=1,nf
         if(i.eq.1)then
            rad0(1)=rfi0
         else
            rad0(i)=rad0(i-1)+drf0
         end if
         do j=1,nzz
            drf(i,j)=drf0
         end do
      end do

!     radii of clad nodes
      drc0=(rco0-rci0)/dfloat(nc-1)
      do i=nf+1,nf+nc
         if(i.eq.nf+1)then
            rad0(nf+1)=rci0
         else
            rad0(i)=rad0(i-1)+drc0
         end if
         do j=1,nzz
            drc(i,j)=drc0
         end do
      end do

      do j=1,nzz
      do i=1,nf+nc
         rad(i,j)=rad0(i)
      end do
      end do

!     initial radii of boundaries between nodes
      do i=1,nf+nc-1
         rad_0(i)=(rad0(i) + rad0(i+1))/2.0
      end do

      do i=1,nf
!        XS area
         if(i.eq.1) then
            az0(1)=pi*((rfi0+drf0/2)**2-rfi0**2)
         else if(i.eq.nf) then
            az0(nf)=pi*(rfo0**2-(rfo0-drf0/2.0)**2)
         else
            az0(i)=pi*(rad0(i)+drf0/2.0)**2 - pi*(rad0(i)-drf0/2.0)**2
         end if
      end do

      do i=nf+1,nf+nc
!        axial area of volume
         if(i.eq.nf+1) then
            az0(i)=pi*((rci0+drc0/2.0)**2-rci0**2)
         else if(i.eq.nf+nc) then
            az0(i)=pi*(rco0**2-(rco0-drc0/2.0)**2)
         else
            az0(i)=pi*(rad0(i)+drc0/2.0)**2 - pi*(rad0(i)-drc0/2.0)**2
         end if
      end do

!     set initial energy parameters
      do j=1,nzz
         enth(j)=0.0d0
         fleak(j)=0.0d0
         centh(j)=0.0d0
      end do

      if(n300001.eq.0)then
         vtot=0.
         h0=0.
         do j=1,nzz
            h0=h0+dz0(j)
         end do
         v=pi*(rfo0**2-rfi0**2)*h0
         vtot=vtot+v
         if(vtot.gt.0.)then
            nqvmlt=1
            tqvmlt(1)=0.d0
            qvmlt(1)=v/vtot
         end if
      end if

      fggen=0.d0
      fgrel=0.d0

!     set initial temperature for fuel rod and volumes

      do j=1,nzz
         tcave(j)=tem0
         tfave(j)=tem0
         do i=1,nf+nc+1
            tem(i,j)=tem0
         end do
      end do

!     initial fuel stack height
      h0=0.d0
      do j=1,nzz
         h0=h0+dz0(j)
      end do
!     as-fabricated free volume
      vfree0=pi*(rci0**2-rfi0**2 + rfo0**2)*h0+vgp

      open(unit=702,file='fred.dat')
      write(702,'(a2,67(a13))')'  ', &
     &             'time(s)      ', &
     &             'time(days)   ', &
     &             'temfi(c)     ', &
     &             'temfo(c)     ', &
     &             'temci(c)     ', &
     &             'temco(c)     ', &
     &             'tfave(c)     ', &
     &             'tcave(c)     ', &
     &             'hgap(w/m2/k) ', &
     &             'rfi(mm)      ', &
     &             'rfo(mm)      ', &
     &             'rci(mm)      ', &
     &             'rco(mm)      ', &
     &             'gap(mkm)     ', &
     &             'gapth(mkm)   ', &
     &             'gpres(mPa)   ', &
     &             'pfc(mpa)     ', &
     &             'htc(w/m2-k)  ', &
     &             'fgrpp(%)     ', &
     &             'efs(%)       ', &
     &             'gap_state    '

      inp(1) = dtout
      inp(2) = tend
      inp(3) = nf
      inp(4) = nc
      inp(5) = nzz
      inp(6) = tem0
      write(*,*)'tem0 ', tem0
      return
!     input error
1000  write(*,*)'input error: line number: ', numlin
      stop
      end