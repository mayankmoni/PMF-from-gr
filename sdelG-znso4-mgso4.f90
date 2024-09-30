program main
implicit none
character(60) damm
integer :: i, file, Atom
integer, parameter :: file_tot=4
integer, parameter :: file_tot1=4
integer, parameter :: file_tot2=4
integer, parameter :: step=1000
integer, parameter :: divide=5
real(8) :: aT,akb,akbT,zinKT

real(8) :: z(file_tot1,step)
real(8) :: wa1, wa2,wa3,wa4,wa11,wa12,wa31,wa32,wa33,wa34
real(8) :: wg31,wg32,wg33,wg34
real(8) :: freq_1(file_tot1,step), freq_2(file_tot2,step)
real(8) :: freq_3(file_tot1,step), freq_4(file_tot2,step) 
real(8) :: freq_5(file_tot,step),z1(file_tot2,step)
real(8) :: freq_6(file_tot,step),freq_7(file_tot,step)
real(8) :: freq_8(file_tot,step),freq_9(file_tot,step),freq_10(file_tot,step)
real(8) :: freq_11(file_tot,step),freq_12(file_tot,step),freq_13(file_tot,step)

real(8) :: av_1(step), av_2(step),av_3(step)
real(8) :: av_11(step), av_12(step),av_tot(step)
real(8) :: eerr(step), eerr2(step),eerr3(step)
real(8) :: sdrr(step), sdrr2(step),sdrr3(step)
real(8) :: eerr11(step), eerr12(step)
real(8) :: sdrr11(step),sdrr12(step),sdtot(step)

real(8) :: agt_1(file_tot1,step),agt_2(file_tot1,step)
real(8) :: agt_3(file_tot1,step),agt_4(file_tot1,step),agt_5(file_tot1,step)
real(8) :: agt_6(file_tot1,step),agt_7(file_tot1,step),agt_8(file_tot1,step),agt_9(file_tot1,step)
real(8) :: agt_10(file_tot1,step),agt_11(file_tot1,step),agt_12(file_tot1,step),agt_13(file_tot1,step)

real(8) :: av_31(step), av_32(step),av_33(step),av_34(step),av_tot3(step)
real(8) :: ag_31(step), ag_32(step),ag_33(step),ag_34(step)
real(8) :: av_tot4(step)
real(8) :: eerr31(step), eerr32(step),eerr33(step),eerr34(step)
real(8) :: sdrr31(step),sdrr32(step),sdrr33(step),sdrr34(step),sdtot3(step)
real(8) :: sdtot4(step)

open(101,file="ZNSO4-2M-rdf-ZN2tZN2-SO1O2O3O4-SO4-OT-HT-H1-H2-SO4-WAT_1.dat")
open(102,file="MGSO4-2M-rdf-MGtoMG-SO1O2O3O4-SO4-OT-HT-H1-H2-SO4-WAT_1.dat")
open(103,file="ZNSO4-MGSO4-1M1M-rdf-ZN2tZN2-MG-SO1O2O3O4-SO4-OT-HT-H1-H2-SO4-WAT_1.dat")
open(104,file="ZNSO4-MGSO4-1M1M-rdf-MGtMG-ZN2-SO1O2O3O4-SO4-OT-HT-H1-H2-SO4-WAT_1.dat")



open(1001,file="adelG-ZNSO4-2M-rdf-ZN2tZN2-SO1O2O3O4-SO4-OT-HT-H1-H2-SO4-WAT.dat")
open(1002,file="adelG-MGSO4-2M-rdf-MGtoMG-SO1O2O3O4-SO4-OT-HT-H1-H2-SO4-WAT.dat")
open(1003,file="adelG-ZNSO4-MGSO4-1M1M-rdf-ZN2tZN2-MG-SO1O2O3O4-SO4-OT-HT-H1-H2-SO4-WAT.dat")
open(1004,file="adelG-ZNSO4-MGSO4-1M1M-rdf-MGtMG-ZN2-SO1O2O3O4-SO4-OT-HT-H1-H2-SO4-WAT.dat")
!###########################33
955 format (13f8.2)
956 format (14f8.2)

akb=8.314
aT=298.0
akbT=akb*aT*0.001
zinKT=1.0/akbT

do Atom=1,1,1

do file=1, file_tot1

        do i=1, step
                if (file.le.2) then 
                read(100*Atom+file,*) z(file,i), freq_1(file,i),freq_2(file,i),freq_3(file,i),freq_4(file,i),&
                freq_5(file,i),freq_6(file,i),freq_7(file,i),freq_8(file,i),freq_9(file,i),freq_10(file,i),&
                freq_11(file,i),freq_12(file,i)!,freq_13(file,i)     

                agt_1(file,i)=-(akbT*zinKT)*log(freq_1(file,i))
                agt_2(file,i)=-(akbT*zinKT)*log(freq_2(file,i))
                agt_3(file,i)=-(akbT*zinKT)*log(freq_3(file,i))
                agt_4(file,i)=-(akbT*zinKT)*log(freq_4(file,i))
                agt_5(file,i)=-(akbT*zinKT)*log(freq_5(file,i))
                agt_6(file,i)=-(akbT*zinKT)*log(freq_6(file,i))
                agt_7(file,i)=-(akbT*zinKT)*log(freq_7(file,i))
                agt_8(file,i)=-(akbT*zinKT)*log(freq_8(file,i))
                agt_9(file,i)=-(akbT*zinKT)*log(freq_9(file,i))
               agt_10(file,i)=-(akbT*zinKT)*log(freq_10(file,i))
               agt_11(file,i)=-(akbT*zinKT)*log(freq_11(file,i))
               agt_12(file,i)=-(akbT*zinKT)*log(freq_12(file,i))
               !agt_13(file,i)=-(akbT*zinKT)*log(freq_13(file,i))
        else

                read(100*Atom+file,*) z(file,i), freq_1(file,i),freq_2(file,i),freq_3(file,i),freq_4(file,i),&
                freq_5(file,i),freq_6(file,i),freq_7(file,i),freq_8(file,i),freq_9(file,i),freq_10(file,i),&
                freq_11(file,i),freq_12(file,i),freq_13(file,i)

                agt_1(file,i)=-(akbT*zinKT)*log(freq_1(file,i))
                agt_2(file,i)=-(akbT*zinKT)*log(freq_2(file,i))
                agt_3(file,i)=-(akbT*zinKT)*log(freq_3(file,i))
                agt_4(file,i)=-(akbT*zinKT)*log(freq_4(file,i))
                agt_5(file,i)=-(akbT*zinKT)*log(freq_5(file,i))
                agt_6(file,i)=-(akbT*zinKT)*log(freq_6(file,i))
                agt_7(file,i)=-(akbT*zinKT)*log(freq_7(file,i))
                agt_8(file,i)=-(akbT*zinKT)*log(freq_8(file,i))
                agt_9(file,i)=-(akbT*zinKT)*log(freq_9(file,i))
               agt_10(file,i)=-(akbT*zinKT)*log(freq_10(file,i))
               agt_11(file,i)=-(akbT*zinKT)*log(freq_11(file,i))
               agt_12(file,i)=-(akbT*zinKT)*log(freq_12(file,i))
               agt_13(file,i)=-(akbT*zinKT)*log(freq_13(file,i))
        end if
        end do
end do


do i=1,step
!eerr31(i)=0.0
!eerr32(i)=0.0
!eerr33(i)=0.0
!eerr34(i)=0.0
!
!do file=1, file_tot1
!       if(file.le.4) then
!       eerr31(i)=eerr31(i)+(ag_31(i)-(agt_3(file,i)))**2
!       eerr32(i)=eerr32(i)+(ag_32(i)-(agt_4(file,i)))**2
!       eerr33(i)=eerr33(i)+(ag_33(i)-(agt_5(file,i)))**2
!       else
!       eerr34(i)=eerr34(i)+(ag_34(i)-(agt_6(file,i)))**2
!       end if
!end do
!       sdrr31(i)=sqrt(eerr31(i)/12)
!       sdrr32(i)=sqrt(eerr32(i)/12)
!       sdrr33(i)=sqrt(eerr33(i)/12)
!       sdrr34(i)=sqrt(eerr34(i)/12)

           write(1000*atom+1,*) z(1,i),agt_1(1,i),agt_2(1,i),agt_3(1,i),agt_4(1,i),agt_4(1,i),agt_5(1,i),&
                   agt_6(1,i),agt_7(1,i),agt_8(1,i),agt_9(1,i),agt_10(1,i),agt_11(1,i),agt_12(1,i)
!
           write(1000*atom+2,*) z(2,i),agt_1(2,i),agt_2(2,i),agt_3(2,i),agt_4(2,i),agt_4(2,i),agt_5(2,i),&
                   agt_6(2,i),agt_7(2,i),agt_8(2,i),agt_9(2,i),agt_10(2,i),agt_11(2,i),agt_12(2,i)
!
           write(1000*atom+3,*) z(3,i),agt_1(3,i),agt_2(3,i),agt_3(3,i),agt_4(3,i),agt_4(3,i),agt_5(3,i),&
                   agt_6(3,i),agt_7(3,i),agt_8(3,i),agt_9(3,i),agt_10(3,i),agt_11(3,i),agt_12(3,i),agt_13(3,i)

           write(1000*atom+4,*) z(4,i),agt_1(4,i),agt_2(4,i),agt_3(4,i),agt_4(4,i),agt_4(4,i),agt_5(4,i),&
                   agt_6(4,i),agt_7(4,i),agt_8(4,i),agt_9(4,i),agt_10(4,i),agt_11(4,i),agt_12(4,i),agt_13(4,i)


           end do
end do
!########################3333333




end program                        


