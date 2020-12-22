!--------------------------------------------------------------------------
!
!    Copyright (C) 1985-2018  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published
!    by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main
!    directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------
!
! sediment dynamics in water for mercury subroutines (sed4merc_water.f)
!
! revision log :
!
! 30.05.2018    dmc&gr	integration of the 0d module into SHYFEM
! 10.07.2019    cl	3d sinking, debug
!
! notes :
!
! dmc= donata canu; gr= ginevra rosati; cl= célia laurent
! laa= Leslie Aveytua Alcazar
!
!*****************************************************************

        subroutine sed4merc_water(bbottom,dtday,tday,wat_vol,
     +                           wdepth,k,temp,sal,taub,area,
     +                          C,Dssink,Dpsink,Vds,Vdp,
     +                          ds_gm2s, dp_gm2s,vold)


       implicit none

       integer nstate
       parameter (nstate=2)


      real C(nstate)         !sed variables: C(1):Silt, C(2):POM  [mg/L]
      real CD(nstate)        !derivatives CD(1) Dsilt, CD(2) DPOM [g/d]
      integer m              !time indicator for write outputs, debug
      integer k              !node
      logical bbottom        !is on bottom?
      real cold(nstate)      !old state variables
      real dtday,tday        !time step [day], time [day]
      real temp,sal          !temperature [C] and salinity [-]
      real tkref             !reference temperature, 293 [K]
      real loads(nstate)     !atmospheric loadings

      real Sw, POMw           !State variables, silt and POM [mg/L]

      real wdepth, area, wat_vol   !depth [m], area [m2] and volume [m3] of water elements
      real vold            !here the same volnew has to be used!
      real taub                    !bottom stress from subssed.f [Pa]
      real Swm,POMwm               !masses of silt and POM in water [g]
      real Vds, Vdp                !Pd x Stoke's velocity for silt and POM              [m/s]
      real ds_gm2s, dp_gm2s        !deposition Flux for silt and POM [g/m2 s]
c      real dep_gm2s                !total deposition flux (silt+POM) [g/m2 s]
      real Dssink,Dpsink,Dsink     !Sink of silt (Dssink), POM (Dpsink)and sum (Dsink) [g/s]

      integer ipext,ipint,kext     !nodes external and internal numbers
      integer fortfilenum, iter
      logical constant_parameters

c	variables
c     _______________________________________________________
c     assigne old value to variables

      Sw=C(1)
      POMw=C(2)        ![mg/l]
c     ________________________________________________________

       ds_gm2s = Vds*Sw           !Flux: [m s-1]*[g m-3]-> [g m-2 s-1]
       dp_gm2s = Vdp*POMw         !Flux
c
       Dssink = ds_gm2s *area     !Sink of silt [g s-1] = [g m-2 s-1] * [m2] =Dsflux di Ginevra
       Dpsink = dp_gm2s *area      !Sink of POM [g s-1] = Dpflux di Ginevra
       Dsink  = Dssink + Dpsink   ! = Dflux di Ginevra
c
c       dep_gm2s = Dsink/area


c       if (kext .EQ. 2284) then
c       write(487,*) Vdp, Vds, 'sed4MERCw'
c       write(486,*) Dssink, Dpsink, 'sed4MERCw'
c       write (500,*) Sw
c       write (501,*) POMw
c       end if

       Swm   = Sw*wat_vol        ! [g m3]*m3 --> g
       POMwm = POMw *wat_vol  ! masses of soilds in water

        if (POMw .LE. 0.0) then  !if
        write(*,*)'instability - negative POMw in sed4merc_wat kext=',
     +  ipext(k),POMw
        stop
        else if (Sw .LE. 0.0) then
        write(*,*)'instability - negative siltw in sed4merc_wat kext=',
     +  ipext(k)
        stop
        end if


        if (C(2) .LE. 0.0) then  !if
        write(*,*)'instability - negative POMw in sed4merc_wat kext=',
     +  ipext(k),POMw
        stop
        else if (C(1) .LE. 0.0) then
        write(*,*)'instability - negative siltw in sed4merc_wat kext=',
     +  ipext(k)
        stop
        end if

c      if (kext .EQ. 70) then
c      write(440,*) bbottom,dtday,tday,wat_vol,wdepth,temp,sal,taub,area
c      write(441,*) C,Dssink,Dpsink,Vds,Vdp
c      write(442,*) ds_gm2s,dp_gm2s,vold
c      end if

C       _________________________________________________________

        C(1)=Sw
        C(2)=POMw

        !write(*,*)C(1), 'c1 s4m_wat'
        !write(*,*)Sw, 'Sw s4m_wat'
c       __________________________________________________________

c         Dssink=Dssink    ! [g/sec] =Dsflux di Ginevra COMCelia: CANCELLED
c         Dsink=Dpsink     ! [g/sec] =Dflux di Ginevra COMCelia: CANCELLED

        CD(1) = -Dssink *86400.   !g/day
        CD(2) = -Dpsink *86400.   !g/day

c      if (kext .EQ. 70) then
c      write(413,*) CD(1),-Dssink*86400.
c      write(414,*) CD(2),-Dpsink*86400.
c      end if

c       call merc_euler (2,dtday,wat_vol,wat_vol,c,cold,cd)    ! c(i)=( c(i)*vol+dt*cd(i) )/vol
        call merc_euler (2,dtday,wat_vol,c,cold,cd,vold)

        if (C(2) .LE. 0.0) then  !if
        write(*,*) 'POMw<=0',POMw,wdepth,ipext(k),'s4m_wat aft'
c        write(449,*) 'POMw<=0',POMw,wdepth,ipext(k),'s4m_wat aft'
c        stop
        else if (C(1) .LE. 0.0) then
        write(*,*) 'siltw<=0',Sw,wdepth,ipext(k),'s4m_wat aft'
c        write(444,*) 'siltw<=0',Sw,wdepth,ipext(k),'s4m_wat aft'
c        stop
        end if

c       if (ipint(k) .Eq. 70) then
c         write(*,*) 'ipext(3987) =', ipext(3987),'k=',k
c       end if
  
c      write(*,*) 'ipint-3985' , ipint(3985), ipext(330)

      iter=nint(tday*86400.)

      if (MOD (iter,1800) .EQ. 0) then
           kext=ipext(k)
           fortfilenum=-1
           if(kext==3985)then       !ve8 closer
               fortfilenum=250
           elseif(kext==3986) then  !ve8
               fortfilenum=251
           elseif(kext==3988) then  !ve8
               fortfilenum=252
           elseif(kext==4007) then  !ve8
               fortfilenum=253
           elseif(kext==3763) then  !ve7
               fortfilenum=254
           elseif(kext==3765) then  !ve7
               fortfilenum=255
           elseif(kext==3764)then   !ve7
               fortfilenum=256
           elseif(kext==3762) then  !ve7 closer
               fortfilenum=257
           else if(kext==2150)then  !ve6 closer
               fortfilenum=258
           elseif(kext==2009) then  !ve6
               fortfilenum=259
           elseif(kext==2359) then  !ve6
               fortfilenum=260
           elseif(kext==2358) then  !ve6
               fortfilenum=261
           elseif(kext==2341) then  !ve5 closer
               fortfilenum=262
           elseif(kext==2408) then  !ve5
               fortfilenum=263
           elseif(kext==2191)then   !ve5
               fortfilenum=264
           elseif(kext==2192) then  !ve5
               fortfilenum=265
           elseif(kext==2654)then   !ve4 closer
               fortfilenum=266
           elseif(kext==2505) then  !ve4
               fortfilenum=267
           elseif(kext==2655) then  !ve4
               fortfilenum=268
           elseif(kext==2653) then  !ve4
               fortfilenum=269
           else if(kext==1372)then  !ve3 closer
               fortfilenum=270
           elseif(kext==1375) then  !ve3
               fortfilenum=271
           elseif(kext==1331) then  !ve3
               fortfilenum=272
           elseif(kext==1378) then  !ve3
               fortfilenum=273
           elseif(kext==4347) then  !ve3
               fortfilenum=274
             elseif(kext==3216) then  !ve2 closer
               fortfilenum=275
           elseif(kext==3057)then   !ve2
               fortfilenum=276
           elseif(kext==2953) then  !ve2
               fortfilenum=277
           elseif(kext==3217) then  !ve2
               fortfilenum=278
           elseif(kext==2405) then  !ve1
               fortfilenum=279
           elseif(kext==2407)then   !ve1
               fortfilenum=280
           elseif(kext==2284) then  !ve1 closer
               fortfilenum=281
           elseif(kext==2404) then  !ve1
               fortfilenum=282
           endif
           if(fortfilenum.ge.0)then
c               if(fortfilenum==250)
c     +             write(*,*) 'stamp to file 250... at iter=',iter,
c     +             ', tday=', tday
               write(fortfilenum,"(2(i10,','),4(f15.7,','))")
     +         iter,kext,wdepth, Sw, POMw, taub
           endif
         endif

c      if (kext .EQ. 70) then
c      write(561,*) C
c      end if



c      write(*,*) 'iteration', iter
c      write(*,*) 'ipint-3985' , ipint(3985), ipext(330)
c      write(*,*)
c      write(*,*) 'ipint-3763', ipint(3763), ipext(929)
c      write(*,*)
c      write(*,*) 'ipint-2009', ipint(2009), ipext(2135)

        end
c---------------------------------------------------------
c       set atmospheric loading on the surface layer
c--------------------------------------------------------
c********************************************************************
c
c      subroutine load0ds(dt,cds,loads,vol)

c integrate loadings

c      implicit none

c        integer nstate          !total number of state parameters
c        parameter( nstate =     2 )

c      real cds(nstate)      !source term [g]
c      real loads(nstate)      !loading for c [g/ day)on the element]
c      real vol            !volume of box [m**3]
c      real dt

c     integer i

c        write(6,*) cds(1),cds(2),cds(3),dt,'cds before load'
c      do i=1,nstate
c        cds(i) = cds(i) +  loads(i)*dt
c      end do

c        write(6,*) loads(1),loads(2),loads(3),dt,'loads'
c        write(6,*) cds(1),cds(2),cds(3),dt,'cds after load'

c      end

c********************************************************************

c	subroutine  sed4merc_gas_exchange(salin,temp,area,visc,rhow)
c
c 4.05.2017	dmc

c	computing mercury exchange between air and the water column
c	Sorensen et al., 2010 An Improved Global Mercury Model
c	for Air Sea Exchange of Mercury: High Concentrations over the North
c	Atlantic


c	implicit none
c       logical bvis1		!select viscosity calculation bvis1 true
c	save bvis1		!bvis true use model output else use Soerensen

c	auxiliar variables
c         real a,b, p,c,d
c         real e,g,h,m,n,o
c         real v1,v2,v3,v4,v5,v6

c	from the hydrodynamic model
c     	real temp	!water temperature °C
c     	real salin	!water salinity
c       	real uwind10	!wind speed normalised at 10 m above sea surface
c     	real area	!element surface
c     	real visc	 ! water viscosity [cP]
c        real rhow 	! density of the (sea)water  (KG/M**3)
c      	real tempk	!temperature [K]

c       tempk=temp+273.15       !temperature Kelvin

C ======================================================================
C ======================================================================
C
C Compute the density and the dynamic viscosity of water from the temperature
C and the salinity

C compute the dynamic/molecular viscosity
c      VISC0=1.802863d-3 - 6.1086d-5*TEMP + 1.31419d-06*TEMP**2 -
c       &1.35576d-08*TEMP**3 + 2.15123d-06*SALIN + 3.59406d-11*SALIN**2

c        a=0.0001529
c        b=0.000016826
c        p=1.013253
c        c=0.000000083885      !8.3885*(10E-8)
c        d=p*p                 !p**(2)
c        e=0.0024727
c        g= 0.000048429       !4.8429*(10E-5)
c        h= 0.0000047172      !4.7172*(10E-6)
c        m= 0.000000075986    !7.5986*(10E-8)
c        n= 0.0000060574     !6.0574*(10E-6)
c        o= 0.000000002676     !2.676*(10E-9)       
c
c        v1= temp*(0.06144-temp*(0.001451-temp*b))
c        v2=a*p
c        v3=c*d
c        v4=e*salin
c        v5=(n*p-o*d)*temp
c        v6=((temp*g)-temp*(h-temp*m))*salin
c
c        visc=(1.791- v1-v2+v3+v4+ v5+v6)/1000.     ![kg m-1* s-1]

c mpute the water density according to Brydon et al. 1999, J. Geoph. Res.
C 104/C1, 1537-1540, equation 2 with Coefficient of Table 4, without pressure
C component. Ranges TEMP -2 - 40øC, S 0-42, surface water.
C      RHOW=9.20601d-2 + 5.10768d-2*TEMP + 8.05999d-1*SALIN
C     &     -7.40849d-3*TEMP**2 - 3.01036d-3*SALIN*TEMP +
C     %     3.32267d-5*TEMP**3 + 3.21931d-5*SALIN*TEMP**2
C      RHOW=RHOW+1000d0

C compute the water density according to EOS80, Fofonoff 198599,
C J. Geoph. Res. 90/C2, 3332-3342, without pressure component.
c	[kg * m-2]

c     RHOW=999.842594d0 +6.793952d-2*TEMP -9.095290d-3*TEMP**2.
c    &   +1.00168d-4*TEMP**3 -1.120083d-6*TEMP**4 +6.536332d-9*TEMP**5.
c    & +(8.24493d-1 -4.0899d-3*TEMP +7.6438d-5*TEMP**2.
c    &   -8.2467d-7*TEMP**3 +5.3875d-9*TEMP**4.) * SALIN
c    & +(-5.72466d-3 +1.0227d-4*TEMP -1.6546d-6*TEMP**2.) * SALIN**1.5d0
c    & +4.8314d-4*SALIN**2.


c      RHOW=999.842594d0+6.793952d-2*TEMP -9.095290d-3*(TEMP*TEMP)
c     & + 1.00168d-4*TEMP**3.-1.120083d-6*TEMP**4.+6.536332d-9*TEMP**5.
c     & + (8.24493d-1 -4.0899d-3*TEMP +7.6438d-5*(TEMP*TEMP)
c     & - 8.2467d-7*TEMP**3. +5.3875d-9*TEMP**4.) * SALIN
c     & + (-5.72466d-3 +1.0227d-4*TEMP -1.6546d-6*(TEMP*TEMP))
c     & * SALIN**1.5d0 +4.8314d-4*(SALIN*SALIN)
c	if(bvis1)then
c
c	kvis=(VISC/RHOW)*10000	![cm2 s-1]
c	else
c
c	kvis=0.017*exp(-0.025*tempk)
c
c	end if

C ======================================================================
C ======================================================================

c     end ! end of routine
