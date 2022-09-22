implicit none
real(16), allocatable:: ar(:,:),ai(:,:),d(:),e(:),e2(:),tau(:,:),zr(:,:),zi(:,:),uno(:,:)
integer:: Nx,Ny,Mx,My,N,ss,p_l,p_r,p_u,p_d,i,j,i_u,i_d,j_l,j_r,k,lx,ly, ierr, b, N_p,N_s,PN,q,M,p
complex(16):: s
real(16):: x,y,pi,t,V,pq,E_d,kT
complex(16)::z
complex*32,allocatable::H1(:,:),H2(:,:),M1(:,:),M2(:,:),HG(:,:),G(:,:),UV(:,:)
real(16),allocatable::UVr(:,:),UVi(:,:),HGr(:,:),HGi(:,:)
complex*32,allocatable::um(:,:),vm(:,:),lm(:),HGo(:,:)
complex*32,allocatable:: HGK(:,:,:,:)
8   format(10f7.2)
9   format(10f12.8)
    open(10, file="gap.txt")

Nx=5; Ny=10
Mx=10; My=5
N=Nx*Ny
allocate(d(0:2*N-1),E(0:2*N-1),e2(0:2*N-1),tau(2,0:2*N-1))
allocate(HG(0:2*N-1, 0:2*N-1),UV(0:2*N-1, 0:2*N-1),HGo(0:2*N-1,0:2*N-1))
allocate(HGi(0:2*N-1,0:2*N-1),HGr(0:2*N-1,0:2*N-1),uno(0:2*N-1,0:2*N-1))
allocate(UVi(0:2*N-1,0:2*N-1),UVr(0:2*N-1,0:2*N-1))
allocate(M1(0:Ny-1,0:Nx-1),M2(0:Ny-1,0:Nx-1))
allocate(H1(0:N-1,0:N-1),H2(0:N-1,0:N-1),G(0:N-1,0:N-1))
allocate(HGK(0:Mx-1,0:My-1,0:2*N-1,0:2*N-1))
allocate(um(0:2*N-1,0:N-1),vm(0:2*N-1,0:N-1),lm(0:N-1))

uno=0.0q0
do i=0,2*N-1
    uno(i,i)=1.0q0
end do


t=1.0q0; V=2.0q0; pq=-0.20q0; kT=0.000000010q0; E_d=1.0q0
pi=4.0q0*atan(1.0q0)
z=(0,1)

do lx=0,Mx-1
do ly=0,My-1
H1=0.0q0;H2=0.0q0;G=0.0q0
do ss=0,N-1
M1=0.0q0; M2=0.0q0
i=ss/Nx; j=mod(ss,Nx)
x=j; y=-i
        i_u=mod(i-1+Ny,Ny); i_d=mod(i+1,Ny); j_l=mod(j-1+Nx,Nx); j_r=mod(j+1,Nx)
        p_l=abs((j-Nx)/Nx); p_r=abs((j+1)/Nx); p_u=abs((i-Ny)/Ny); p_d=abs((i+1)/Ny)

        M1(i,j_l)=-t*exp(p_l*2*pi*z*(-1.0q0*lx/Mx+y/Ny))
        M1(i,j_r)=-t*exp(p_r*2*pi*z*(1.0q0*lx/Mx-y/Ny))
        M1(i_u,j)=-t*exp(2*z*pi*x/N)*exp(p_u*2*pi*z*ly/My)
        M1(i_d,j)=-t*exp(-2*z*pi*x/N)*exp(-p_d*2*pi*z*ly/My)

        M2(i,j_l)=-t*exp(p_l*2*pi*z*(-1.0q0*lx/Mx-y/Ny))
        M2(i,j_r)=-t*exp(p_r*2*pi*z*(1.0q0*lx/Mx+y/Ny))
        M2(i_u,j)=-t*exp(-2*pi*z*x/N)*exp(p_u*2*pi*z*ly/My)
        M2(i_d,j)=-t*exp(2*pi*z*x/N)*exp(-p_d*2*pi*z*ly/My)

        H1(ss,:)=reshape(transpose(M1),(/N/))
        H2(ss,:)=reshape(transpose(M2),(/N/))
        H1(ss,ss)=-pq
        H2(ss,ss)=-pq
        G(ss,ss)=0.10q0
end do
    HG=0.0q0
    HG(0:N-1,0:N-1)=H1
    HG(0:N-1,N:2*N-1)=G
    HG(N:2*N-1,0:N-1)=conjg(G)
    HG(N:2*N-1,N:2*N-1)=-H2

    HGK(lx,ly,:,:)=HG
end do
end do

!do i=0,N-1
!    write(*,9) HG(i,:)
!end do


b=0
N_p=4
M=Mx*My !Aquí hay que estár modificando
PN=N_p*M
do p=0,N_p-1
    UV=0.0q0
    lm=0.0q0
    do lx=0,Mx-1
        do ly=0,My-1
            HGo=HGK(lx,ly,:,:)
            HGr=real(HGo); HGi=imag(HGo)
            call htridi(2*N,2*N,HGr,HGi,D,E,E2,tau)
            UVr=uno
            call tql2(2*N,2*N,D,E,UVr,ierr)
            call htribk(2*N,2*N,HGr,HGi,tau,2*N,UVr,UVi)      !Aquí sacas los nuevos eigenvalores
            UV=UVr+z*UVi
            UV=transpose(UV)
            um=UV(0:2*N-1,0:N-1)
            vm=UV(0:2*N-1,N:2*N-1)

            do i=0,N-1
                do q=0,2*N-1
                    if(abs(D(q))<=E_d) then
                        lm(i)=lm(i)+um(q,i)*conjg(vm(q,i))*tanh(0.5*D(q)/kT)
                    end if
                end do
            end do

            b=b+1
            write(*,8) (1.0*b/PN)*100

        end do
    end do

    do lx=0,Mx-1
        do ly=0,My-1
            H1=0.0q0;H2=0.0q0;G=0.0q0
            do ss=0,N-1
                M1=0.0q0; M2=0.0q0
                i=ss/Nx; j=mod(ss,Nx)
                x=j; y=-i
                i_u=mod(i-1+Ny,Ny); i_d=mod(i+1,Ny); j_l=mod(j-1+Nx,Nx); j_r=mod(j+1,Nx)
                p_l=abs((j-Nx)/Nx); p_r=abs((j+1)/Nx); p_u=abs((i-Ny)/Ny); p_d=abs((i+1)/Ny)
!
                M1(i,j_l)=-t*exp(p_l*2*pi*z*(-1.0q0*lx/Mx+y/Ny))
                M1(i,j_r)=-t*exp(p_r*2*pi*z*(1.0q0*lx/Mx-y/Ny))
                M1(i_u,j)=-t*exp(2*z*pi*x/N)*exp(p_u*2*pi*z*ly/My)
                M1(i_d,j)=-t*exp(-2*z*pi*x/N)*exp(-p_d*2*pi*z*ly/My)

                M2(i,j_l)=-t*exp(p_l*2*pi*z*(-1.0q0*lx/Mx-y/Ny))
                M2(i,j_r)=-t*exp(p_r*2*pi*z*(1.0q0*lx/Mx+y/Ny))
                M2(i_u,j)=-t*exp(-2*pi*z*x/N)*exp(p_u*2*pi*z*ly/My)
                M2(i_d,j)=-t*exp(2*pi*z*x/N)*exp(-p_d*2*pi*z*ly/My)

                H1(ss,:)=reshape(transpose(M1),(/N/))
                H2(ss,:)=reshape(transpose(M2),(/N/))
                H1(ss,ss)=-pq
                H2(ss,ss)=-pq
                G(ss,ss)=(V/(2.0q0*M))*lm(ss)

            end do
            HG=0.0q0
            HG(0:N-1,0:N-1)=H1
            HG(0:N-1,N:2*N-1)=G
            HG(N:2*N-1,0:N-1)=conjg(G)
            HG(N:2*N-1,N:2*N-1)=-H2

            HGK(lx,ly,:,:)=HG
       end do
    end do
end do

do i=0,N-1
    write(10,9) real(G(i,i)),imag(G(i,i))
end do

end program

!==============================================================================
!                       SUBRUTINAS DE DIAGONALIZACION
!==============================================================================
      subroutine htridi(nm,n,ar,ai,d,e,e2,tau)
      integer i,j,k,l,n,ii,nm,jp1
      real(16) ar(nm,n),ai(nm,n),d(n),e(n),e2(n),tau(2,n)
      real(16) f,g,h,fi,gi,hh,si,scale,pythag
      tau(1,n) = 1.0q0
      tau(2,n) = 0.0q0
      do 100 i = 1, n
100   d(i) = ar(i,i)
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0q0
         scale = 0.0q0
         if (l .lt. 1) go to 130
         do 120 k = 1, l
120    scale = scale + abs(ar(i,k)) + abs(ai(i,k))
         if (scale .ne. 0.0q0) go to 140
         tau(1,l) = 1.0q0
         tau(2,l) = 0.0q0
130    e(i) = 0.0q0
         e2(i) = 0.0q0
         go to 290
140    do 150 k = 1, l
            ar(i,k) = ar(i,k) / scale
            ai(i,k) = ai(i,k) / scale
            h = h + ar(i,k) * ar(i,k) + ai(i,k) * ai(i,k)
150    continue
         e2(i) = scale * scale * h
         g = sqrt(h)
         e(i) = scale * g
         f = sqrt(ar(i,l)*ar(i,l)+ai(i,l)*ai(i,l))
         if (f .eq. 0.0q0) go to 160
         tau(1,l) = (ai(i,l) * tau(2,i) - ar(i,l) * tau(1,i)) / f
         si = (ar(i,l) * tau(2,i) + ai(i,l) * tau(1,i)) / f
         h = h + f * g
         g = 1.0q0 + g / f
         ar(i,l) = g * ar(i,l)
         ai(i,l) = g * ai(i,l)
         if (l .eq. 1) go to 270
         go to 170
160    tau(1,l) = -tau(1,i)
         si = tau(2,i)
         ar(i,l) = g
170    f = 0.0q0
         do 240 j = 1, l
            g = 0.0q0
            gi = 0.0q0
            do 180 k = 1, j
               g = g + ar(j,k) * ar(i,k) + ai(j,k) * ai(i,k)
               gi = gi - ar(j,k) * ai(i,k) + ai(j,k) * ar(i,k)
180       continue
            jp1 = j + 1
            if (l .lt. jp1) go to 220
            do 200 k = jp1, l
               g = g + ar(k,j) * ar(i,k) - ai(k,j) * ai(i,k)
               gi = gi - ar(k,j) * ai(i,k) - ai(k,j) * ar(i,k)
200       continue
220       e(j) = g / h
            tau(2,j) = gi / h
            f = f + e(j) * ar(i,j) - tau(2,j) * ai(i,j)
240    continue
         hh = f / (h + h)
         do 260 j = 1, l
            f = ar(i,j)
            g = e(j) - hh * f
            e(j) = g
            fi = -ai(i,j)
            gi = tau(2,j) - hh * fi
            tau(2,j) = -gi
            do 260 k = 1, j
               ar(j,k) = ar(j,k)-f*e(k) - g * ar(i,k) + fi * tau(2,k) + gi * ai(i,k)
               ai(j,k) = ai(j,k)-f*tau(2,k) - g * ai(i,k) - fi * e(k) - gi * ar(i,k)
260    continue
270    do 280 k = 1, l
            ar(i,k) = scale * ar(i,k)
            ai(i,k) = scale * ai(i,k)
280    continue
         tau(2,l) = -si
290    hh = d(i)
         d(i) = ar(i,i)
         ar(i,i) = hh
         ai(i,i) = scale * sqrt(h)
300   continue
      return
      end
!-------------------------------------------------------------------
      subroutine tql2(nm,n,d,e,z,ierr)
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      real(16) d(n),e(n),z(nm,n)
      real(16) c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
      ierr = 0
      if (n .eq. 1) go to 1001
      do 100 i = 2, n
100   e(i-1) = e(i)
      f = 0.0q0
      tst1 = 0.0q0
      e(n) = 0.0q0
      do 240 l = 1, n
         j = 0
         h = abs(d(l)) + abs(e(l))
         if (tst1 .lt. h) tst1 = h
         do 110 m = l, n
            tst2 = tst1 + abs(e(m))
            if (tst2 .eq. tst1) go to 120
110    continue
120    if (m .eq. l) go to 220
130    if (j .eq. 30) go to 1000
         j = j + 1
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0q0 * e(l))
         r = sqrt(p*p+1.0q0*1.0q0)
         d(l) = e(l) / (p + sign(r,p))
         d(l1) = e(l) * (p + sign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
         do 140 i = l2, n
140    d(i) = d(i) - h
145    f = f + h
         p = d(m)
         c = 1.0q0
         c2 = c
         el1 = e(l1)
         s = 0.0q0
         mml = m - l
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = sqrt(p*p+e(i)*e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
180       continue
200    continue
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + abs(e(l))
         if (tst2 .gt. tst1) go to 130
220    d(l) = d(l) + f
240   continue
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
260    continue
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
280    continue
300   continue
      go to 1001
1000  ierr = l
1001  return
      end
!--------------------------------------------------------------------
      subroutine htribk(nm,n,ar,ai,tau,m,zr,zi)
      integer i,j,k,l,m,n,nm
      real(16) ar(nm,n),ai(nm,n),tau(2,n),zr(nm,m),zi(nm,m)
      real(16) h,s,si
      if (m .eq. 0) go to 200
      do 50 k = 1, n
         do 50 j = 1, m
            zi(k,j) = -zr(k,j) * tau(2,k)
            zr(k,j) = zr(k,j) * tau(1,k)
50    continue
      if (n .eq. 1) go to 200
      do 140 i = 2, n
         l = i - 1
         h = ai(i,i)
         if (h .eq. 0.0q0) go to 140
         do 130 j = 1, m
            s = 0.0q0
            si = 0.0q0
            do 110 k = 1, l
               s = s + ar(i,k) * zr(k,j) - ai(i,k) * zi(k,j)
               si = si + ar(i,k) * zi(k,j) + ai(i,k) * zr(k,j)
110       continue
            s = (s / h) / h
            si = (si / h) / h
            do 120 k = 1, l
               zr(k,j) = zr(k,j) - s * ar(i,k) - si * ai(i,k)
               zi(k,j) = zi(k,j) - si * ar(i,k) + s * ai(i,k)
120       continue
130    continue
140    continue
200    return
      end

