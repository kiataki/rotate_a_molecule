!!! Matheus Kiataki(2017), reecrito em 2019

program rotacoes      
 implicit none

real*8,parameter :: angstobohr = (1.0d0/0.52917724924d0) !gamess value 
!!real*8 :: fatorc=1.889725989D0
real*8, allocatable, dimension(:) :: theta
real*8, allocatable, dimension(:,:) :: coord
real*8, parameter :: pi=3.14159265358979D0  !!!pi=3.141592654
real*8, dimension(3,3) :: Rm 
real, allocatable, dimension(:) :: ne
integer :: na, i, j, op, stat, nopera, unitt
character(9), allocatable, dimension(:) :: atom
character(7) :: flag
character(5) :: group
character(1), allocatable, dimension(:) :: eixo

read(5,'(a,1x,i3)',iostat=stat) flag, na 
call flagerror(flag,'NATOM =')
call sstatus(stat,'flag, na')

allocate(atom(na))
allocate(ne(na))
allocate(coord(3,na))

read(5,'(a,1x,i3)',iostat=stat) flag, nopera 
call flagerror(flag,'NOPER =')
call sstatus(stat,'flag, nopera') 

allocate(eixo(nopera))
allocate(theta(nopera))

read(5,'(a,1x)', advance='no', iostat=stat) flag 
call flagerror(flag,'EIXOS =')
call sstatus(stat,'flag, eixo(op), op=1,nopera') 
read(5,*,iostat=stat) (eixo(op), op=1,nopera)
call sstatus(stat,'flag, eixo(op), op=1,nopera') 

read(5,'(a,1x)', advance='no', iostat=stat) flag
call flagerror(flag,'THETA =')
call sstatus(stat,'flag') 
read(5,*,iostat=stat) (theta(op), op=1,nopera)
call sstatus(stat,'flag, theta(op), op=1,nopera') 

read(5,'(a,1x,a)', iostat=stat) flag, group
call flagerror(flag,'GROUP =')
call sstatus(stat,'flag, group') 

do i=1,na
  read(5,*,iostat=stat) atom(i), ne(i), (coord(j,i), j=1,3)
  call sstatus(stat,'atom(i),ne(i),(coord(j,i), j=1,3)')
end do

open(9,file='antes-da-rotacao-distancia-interatomica.dat')
write(9,*)'distancia inter atomica antes da rotacao em angstrom'
write(9,*)
do i=1,na-1
  do j=i+1,na
    write(9,*) atom(i), atom(j), sqrt( (coord(1,i)-coord(1,j))**2 + (coord(2,i)-coord(2,j))**2 + (coord(3,i)-coord(3,j))**2 ) 
  end do
end do

theta(:)=((theta(:)*pi)/180.0D0)

do op=1,nopera      
  if(eixo(op).eq.'y') then
     Rm(1,1)=dcos(theta(op));  Rm(1,2)=0.0D0; Rm(1,3)=dsin(theta(op))
     Rm(2,1)=0.0D0;            Rm(2,2)=1.0D0; Rm(2,3)=0.0D0 
     Rm(3,1)=-dsin(theta(op)); Rm(3,2)=0.0D0; Rm(3,3)=dcos(theta(op))
  elseif(eixo(op).eq.'x') then
     Rm(1,1)=1.0D0; Rm(1,2)=0.0D0;           Rm(1,3)=0.0D0
     Rm(2,1)=0.0D0; Rm(2,2)=dcos(theta(op)); Rm(2,3)=-dsin(theta(op))
     Rm(3,1)=0.0D0; Rm(3,2)=dsin(theta(op)); Rm(3,3)=dcos(theta(op))   
  elseif(eixo(op).eq.'z') then
     Rm(1,1)=dcos(theta(op)); Rm(1,2)=-dsin(theta(op)); Rm(1,3)=0.0D0
     Rm(2,1)=dsin(theta(op)); Rm(2,2)=dcos(theta(op));  Rm(2,3)=0.0D0
     Rm(3,1)=0.0D0;           Rm(3,2)=0.0D0;            Rm(3,3)=1.0D0   
  else
     stop'Abort, you put invalid axes'
  end if
  coord = matmul(Rm,coord)
end do

open(10,file='coord-giradas-angs.inp')
open(11,file='coord-giradas-bohr.inp')
open(12,file='depois-da-rotacao-distancia-inter-atomica.dat')

do unitt=10,11
  if(unitt.eq.10) write(unitt,'(1x,a,1x,a,1x,a,1x,a,1x,a)') '$CONTRL', 'SCFTYP=RHF', 'COORD=CART', 'UNITS=ANGS', '$END'
  if(unitt.eq.11) write(unitt,'(1x,a,1x,a,1x,a,1x,a,1x,a)') '$CONTRL', 'SCFTYP=RHF', 'COORD=CART', 'UNITS=BOHR', '$END'
  write(unitt,'(1x,a,2x,a)') '$SYSTEM', '$END'
  write(unitt,'(1x,a,5x,a)') '$SCF', '$END'
  write(unitt,'(1x,a,2x,a)') '$STATPT', '$END'
  write(unitt,'(1x,a,3x,a)') '$BASIS', '$END'
  write(unitt,'(1x,a)') '$DATA'
  write(unitt,'(a)')'...' 
  write(unitt,'(a)') group
  if(group.ne.'C1') write(unitt,*)
  do i=1,na
    if(unitt.eq.10) then
      write(unitt,'(a,11x,f4.1,3(f15.9,2x))') atom(i), ne(i), (coord(j,i), j=1,3)
    else
      write(unitt,'(a,11x,f4.1,3(f15.9,2x))') atom(i), ne(i), (coord(j,i)*angstobohr, j=1,3)     
    end if
  end do 
  write(unitt,'(1x,a)') '$END'
end do

write(12,*)'distancia inter atomica depois da rotacao em angstrom'
write(12,*)
do i=1,na-1
  do j=i+1,na
    write(12,*) atom(i), atom(j), sqrt( (coord(1,i)-coord(1,j))**2 + (coord(2,i)-coord(2,j))**2 + (coord(3,i)-coord(3,j))**2 ) 
  end do
end do
!write(*,*) pi
end program rotacoes
!==========================================================================================================
subroutine flagerror(flag2,msgflag)

character(len=*), intent(in) :: flag2
character(len=*), intent(in) :: msgflag

if(flag2.ne.msgflag) then
  write(*,'(4a)')'Abort, ',flag2,' .ne. ',msgflag 
  stop
end if

end subroutine flagerror
!=========================================================================================================
subroutine sstatus(stat,msg)

character(len=*), intent(in) :: msg

if(stat.gt.0)then
  write(*,'(2a)')'Abort, problem reading: ',msg
  stop
elseif(stat.lt.0) then
  write(*,'(2a)')'Abort, end of input file: ',msg
  stop
end if  

end subroutine sstatus

