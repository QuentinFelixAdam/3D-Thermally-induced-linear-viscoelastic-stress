PROGRAM thermalstress()
IMPLICIT NONE
INTEGER :: nlinesVE,i,j,m
REAL(KIND(1.0D0))::stressvaluesubrout,E0,Einf,nD,TauD,c1VE,c2VE
REAL(KIND(1.0D0))::AlphaT,Tr,Nu,factor
REAL(KIND(1.0D0)),DIMENSION(:),ALLOCATABLE::listVEsubroutlocal,listVEsubrout
REAL(KIND(1.0D0)),DIMENSION(:),ALLOCATABLE::list_Delta_T_t
REAL(KIND(1.0D0)),DIMENSION(:),ALLOCATABLE::list_Delta_T_tau
REAL(KIND(1.0D0)),DIMENSION(:),ALLOCATABLE::list_d_Delta_T_tau
REAL(KIND(1.0D0)),DIMENSION(:),ALLOCATABLE::list_d_Delta_T_t
REAL(KIND(1.0D0)),DIMENSION(:),ALLOCATABLE::list_inverse_aT
REAL(KIND(1.0D0)),DIMENSION(:),ALLOCATABLE::list_integrand1
REAL(KIND(1.0D0)),DIMENSION(:),ALLOCATABLE::list_ksi_prime
REAL(KIND(1.0D0)),DIMENSION(:),ALLOCATABLE::list_integrand2
REAL(KIND(1.0D0)),DIMENSION(:),ALLOCATABLE::list_VE_Stress
REAL(KIND(1.0D0)),DIMENSION(:),ALLOCATABLE::list_ksi_inter
REAL(KIND(1.0D0)),DIMENSION(:),ALLOCATABLE::list_ksi_inter_rev
REAL(KIND(1.0D0)),DIMENSION(:),ALLOCATABLE::list_diff_ksi
REAL(KIND(1.0D0)),DIMENSION(:),ALLOCATABLE::list_diff_ksi_rev
REAL(KIND(1.0D0)),DIMENSION(:),ALLOCATABLE::list_relax

E0 = 30000.0
Einf = 400.0
nD = 0.3
TauD = 10.0**7
c1VE = 30.0
c2VE = 200.0
AlphaT = 0.00002
Tr = 0.0
Nu = 0.25
factor = AlphaT/(1.0-Nu)
nlinesVE = 0

OPEN(51,FILE='Surfacetemp.dat',STATUS='old')
DO
   READ(51,*,END=19)
   nlinesVE = nlinesVE + 1
END DO
19 CLOSE(51)

ALLOCATE(listVEsubroutlocal(nlinesVE))

OPEN(51,FILE='Tsurface.dat',STATUS='old')

DO i = 1,nlinesVE
   READ(51,*) listVEsubroutlocal(i)
END DO

ALLOCATE(list_Delta_T_t(nlinesVE))
ALLOCATE(list_Delta_T_tau(nlinesVE-1))
ALLOCATE(list_d_Delta_T_tau(nlinesVE-1))
ALLOCATE(list_d_Delta_T_t(nlinesVE-1))
ALLOCATE(list_inverse_aT(nlinesVE-1))
ALLOCATE(list_integrand1(nlinesVE-1))
ALLOCATE(list_ksi_prime(nlinesVE-1))
ALLOCATE(list_integrand2(nlinesVE-1))
ALLOCATE(list_VE_Stress(nlinesVE-1))
ALLOCATE(list_ksi_inter(nlinesVE-2))
ALLOCATE(list_ksi_inter_rev(nlinesVE-2))
ALLOCATE(list_diff_ksi(nlinesVE-2))
ALLOCATE(list_diff_ksi_rev(nlinesVE-2))
ALLOCATE(list_relax(nlinesVE-2))


DO i = 1,nlinesVE
   list_Delta_T_t(i) = listVEsubroutlocal(i) - Tr
END DO

DO i = 1,nlinesVE - 1
   list_d_Delta_T_t(i) = (list_Delta_T_t(i+1) - list_Delta_T_t(i))
END DO

DO i = 1,nlinesVE - 1
   list_Delta_T_tau(i) = (list_Delta_T_t(i) + list_Delta_T_t(i+1))*0.5
END DO

DO i = 1,nlinesVE - 1
   list_d_Delta_T_tau(i) = (list_Delta_T_tau(i+1) - list_Delta_T_tau(i))
END DO

DO i = 1,nlinesVE - 1
   list_inverse_aT(i) = 10**((-c1VE*list_Delta_T_tau(i)/(c2VE + list_Delta_T_tau(i))))
   list_inverse_aT(i) = 1.0 / list_inverse_aT(i)
END DO

DO i = 1,nlinesVE - 1
   list_integrand1(i) = 0.5*60.0*(list_inverse_aT(i+1) + list_inverse_aT(i))
END DO

DO i = 2,nlinesVE - 1
   list_ksi_prime(i) = list_ksi_prime(i-1) + list_integrand1(i)
END DO

list_ksi_inter = list_ksi_prime(1:nlinesVE-2)

list_ksi_inter_rev = list_ksi_inter(nlinesVE-2:1:-1)

DO j = 1,nlinesVE-2
   list_diff_ksi(j) = list_ksi_inter(nlinesVE-2) - list_ksi_inter_rev(j)
END DO

list_diff_ksi_rev = list_diff_ksi(nlinesVE-2:1:-1)

DO j = 1,nlinesVE-2
   list_relax(j) = Einf * (1.0 + (list_diff_ksi_rev(j)/TauD)**nD)/ & 
   ((list_diff_ksi_rev(j)/TauD)**nD + Einf/E0)
END DO
   

DO j = 1,nlinesVE-2
   list_integrand2(j) = list_d_Delta_T_t(j) * list_relax(j)
END DO

stressvaluesubrout = SUM(list_integrand2)*factor

OPEN(28,FILE='VEStress.dat',POSITION='append')

WRITE(28,*) stressvaluesubrout

CLOSE(28)

DEALLOCATE(list_ksi_inter)
DEALLOCATE(list_ksi_inter_rev)
DEALLOCATE(list_diff_ksi)
DEALLOCATE(list_diff_ksi_rev)
DEALLOCATE(list_relax)
DEALLOCATE(list_Delta_T_t)
DEALLOCATE(list_Delta_T_tau)
DEALLOCATE(list_d_Delta_T_tau)
DEALLOCATE(list_d_Delta_T_t)
DEALLOCATE(list_inverse_aT)
DEALLOCATE(list_integrand1)
DEALLOCATE(list_ksi_prime)
DEALLOCATE(list_integrand2)
DEALLOCATE(list_VE_Stress)

END PROGRAM thermalstress
