program Hemholtz
implicit none
! Define variable in main code
real:: dx,dy,gamma,eps,delta,rho,p,M,tf,CFL,dt,time,residual,emax,output
real,dimension(4) :: phi,psi,qr,ql,fe,gn
real,dimension(:),allocatable :: acoustic
real,dimension(:,:),allocatable ::x,y,alpha
real,dimension(:,:,:),allocatable ::q,q0,rhs
integer i,j,ni,nj,pounter,n stage,stage,step,nm
character(20) :: method
! Grid
ni=360
nj=360
dx=2.*3.14159/ni
dy=2.*3.14159/nj
! x=dx/2:dx:2*3.14159−dx/2; y=dy/2:dy:2*3.14159−dy/2;
! Flow parameters
gamma=1.4
eps=3.14159/15
delta=0.05
rho=1.
p=1/gamma
M=0.25
tf=30.
CFL=0.4
method='van Leer';
residual=0.;
allocate(x(ni,nj));x=0.
allocate(y(ni,nj));y=0.
do i=1,ni
do j=1,nj
x(i,j)=dx/2.+(i−1.)*dx;
y(i,j)=dy/2.+(j−1.)*dy;
end do
end do
allocate(q(4,ni+4,nj+4)); q = 0.
allocate(q0(4,ni+4,nj+4)); q0 = 0.
allocate(acoustic(ni*nj)); acoustic = 0.
allocate(rhs(4,ni+4,nj+4))
allocate(alpha(3,3));alpha=0.
pounter=1;
do i=3,ni+2
do j=3,nj+2
q(1,i,j)=rho;
if (j>nj/2) then
q(2,i,j)=rho*M*tanh((3*3.14159−2*y(i−2,j−2))/(2*eps));
else
q(2,i,j)=rho*M*tanh((2*y(i−2,j−2)−3.14159)/(2*eps));
end if
q(3,i,j)=M*delta*sin(x(i−2,j−2));
q(4,i,j)=p/(gamma−1)+0.5*(rho*M**2);
acoustic(pounter)=sqrt(q(1,i,j)**gamma/q(1,i,j))+sqrt((q(2,i,j)**2.+q(3,i,j)**2)/q(1,i,j)**2.)
pounter=pounter+1;
end do
end do
do i=2,pounter−1
emax=max(acoustic(i),acoustic(i−1))
enddo
dt=CFL*dx/emax
nm=anint(tf/dt)
! Debuggggg
!open(unit=15,file='debuggggg.txt',status='unknown')
!
! do j = 3,nj+2
! do i = 3,ni+2
! write(15,'( 6(:,1p,1x,e24.16e3) )')x(i,j), y(i,j), q(1,i,j),q(2,i,j),q(3,i,j),q(4,i,j)
! end do
! end do
! close(15)
! 3−Stage TVD R−K Coeff
n stage=3;
alpha(1,:)=[1., 3./4., 1./3.]
alpha(2,:)=[0., 1./4., 2./3.]
alpha(3,:)=[1., 1./4., 2./3.]
!! Evolve state in time
do step=1,nm
! output=step/nm
write(*,*) step
q0=q;
do stage=1,n stage
! Periodic Boundary Conditions at Ghost Points
do j=3,nj+2
q(:,2,j)=q(:,ni+2,j)
q(:,1,j)=q(:,ni+1,j)
q(:,ni+3,j)=q(:,3,j)
q(:,ni+4,j)=q(:,4,j)
enddo
do i=3,ni+2
q(:,i,2)=q(:,i,nj+2)
q(:,i,1)=q(:,i,nj+1)
q(:,i,nj+3)=q(:,i,3)
q(:,i,nj+4)=q(:,i,4)
enddo
! RHS contribution
rhs = 0.
! write(*,*) 'broke here?'
do j=3,nj+2
do i=2,ni+2
call limiter(q(:,i+1,j)−q(:,i,j),q(:,i,j)−q(:,i−1,j),method,phi)
call limiter(q(:,i+1,j)−q(:,i,j),q(:,i+2,j)−q(:,i+1,j),method,psi)
qL=q(:,i,j) + 0.5*phi*(q(:,i,j)−q(:,i−1,j))
qR=q(:,i+1,j) − 0.5*psi*(q(:,i+2,j)−q(:,i+1,j))
call ausmp(gamma,1.,0.,qL,qR,fE)
rhs(:,i,j)=rhs(:,i,j)+fE*dy
rhs(:,i+1,j)=rhs(:,i+1,j)−fE*dy
enddo
enddo
do j=2,nj+2
do i=3,ni+2
call limiter(q(:,i,j+1)−q(:,i,j),q(:,i,j)−q(:,i,j−1),method,phi)
call limiter(q(:,i,j+1)−q(:,i,j),q(:,i,j+2)−q(:,i,j+1),method,psi)
qL=q(:,i,j) + 0.5*phi*(q(:,i,j)−q(:,i,j−1))
qR=q(:,i,j+1) − 0.5*psi*(q(:,i,j+2)−q(:,i,j+1))
call ausmp(gamma,0.,1.,qL,qR,gN);
rhs(:,i,j)=rhs(:,i,j)+gN*dx
rhs(:,i,j+1)=rhs(:,i,j+1)−gN*dx
enddo
enddo
!Runge−Kutta Stage Update of Interior State
do j=3,nj+2
do i=3,ni+2
q(:,i,j)=alpha(1,stage)*q0(:,i,j)+alpha(2,stage)*q(:,i,j)−alpha(3,stage)*dt*rhs(:,i,j)/(dx*dy)
enddo
enddo
enddo
! monitor residual
time=step*dt;
residual=0.;
do i=3,ni+2
do j=3,nj+2
residual=residual+rhs(1,i,j)**2;
enddo
enddo
residual=sqrt(residual)
! output solution
! intime r{step}=q(3:ni+2,3:nj,1);
! intime u{step}=q(3:ni+2,3:nj,2);
! intime v{step}=q(3:ni+2,3:nj,3);
! intime e{step}=q(3:ni+2,3:nj,4);
!write(*,*) 'broke here???'
! Write the output to a file.
if(mod(step,25)==0) then
open(unit=12,file='kh total.txt',status='unknown',position='append')
do j = 3,nj+2
do i = 3,ni+2
write(12,'( 6(:,1p,1x,e24.16e3) )')x(i−2,j−2), y(i−2,j−2), q(1,i,j),q(2,i,j),q(3,i,j),q(4,i,j)
end do
end do
close(12)
elseif (step==1) then
open(unit=12,file='kh total.txt',status='unknown',position='append')
! write(*,*) 'broke here??????'
do j = 3,nj+2
do i = 3,ni+2
write(12,'( 6(:,1p,1x,e24.16e3) )')x(i−2,j−2), y(i−2,j−2), q(1,i,j),q(2,i,j),q(3,i,j),q(4,i,j)
end do
end do
close(12)
endif
enddo
stop
end
subroutine ausmp(g,nx,ny,left,right,flux)
! Define parameters
implicit none
real :: g,nx,ny,alpha,beta,p left,p right,a left,a right,m,un left,un right,m plus,m minus,p plus,p minus,m face,p face,a face
real,dimension(4)::left,right,flux,state
! Default
alpha=3./16.
beta=1./8.
! Left Side Primitive
state(1)=left(1)
state(2)=left(2)
state(3)=left(3)
state(4)=left(4)
un left=state(2)/state(1)*nx+state(3)/state(1)*ny;
p left=(g−1)*(state(4)−.5*(state(2)/state(1)*state(2)+state(3)/state(1)*state(3)))
if (p left <0) then
write(*,*) 'Negative pressure from left state'
endif
left(4)=left(4)+p left
a left=sqrt(2*(g−1.)/(g+1)*left(4)/left(1))
! Right side Primitive
state(1)=right(1)
state(2)=right(2)
state(3)=right(3)
state(4)=right(4)
un right=state(2)/state(1)*nx+state(3)/state(1)*ny
p right=(g−1)*(state(4)−.5*(state(2)/state(1)*state(2)+state(3)/state(1)*state(3)))
if (p right <0) then
write(*,*) 'Negative pressure from right state'
endif
right(4)=right(4)+p right
a right=sqrt(2*(g−1)/(g+1)*right(4)/right(1))
! Face acoustic speed, Eq 40
a face=min(a left**2/max(a left,abs(un left)),a right**2/max(a right,abs(un right)))
! Mach number and pressure splitting, Eq 19 & 21, Alg A2
m=un left/a face
if (abs(m)>=1) then
m plus=1./2.*(m+abs(m))
p plus=0.5*(1.+sign(1.,m))
else
m plus=1./4.*(m+1.)**2.+beta*(m**2.−1.)**2.
p plus=1./4.*(m+1.)**2.*(2.−m)+alpha*m*(m**2.−1.)**2.
endif
m=un right/a face
if (abs(m)>=1) then
m minus=1./2.*(m−abs(m))
p minus=0.5*(1.−sign(1.,m))
else
m minus=−1./4.*(m−1.)**2.−beta*(m**2.−1.)**2.
p minus=1./4.*(m−1.)**2.*(2.+m)−alpha*m*(m**2.−1.)**2.
endif
m face=m plus+m minus
p face=p plus*p left+p minus*p right
! Flux
flux=a face*(0.5*(m face+abs(m face))*left+0.5*(m face−abs(m face))*right)
flux(2)=flux(2)+p face*nx
flux(3)=flux(3)+p face*ny
return
end
subroutine limiter(a,b,method,phi)
implicit none
real,parameter :: e=1.e−7
real,dimension(4) :: a,b,r,phi
character(20) method
select case (method)
case ('First Order')
phi=0.
case ('Second Order')
phi=1.
case ('MinMod')
r=(a+e)/(b+e)
phi=max(0.,min(1.,r))
case ('Superbee')
r=(a+e)/(b+e)
phi=max(0.,min(2.*r,1.),min(r,2.))
case ('van Leer')
r=(a+e)/(b+e)
phi=(r+abs(r))/(1.+abs(r))
case default
write(*,*)'Unknown limiter'
stop
end select
return
end

