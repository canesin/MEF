      subroutine elmq01(e,x,u,p,s,nel,ndm,nst,isw,nin)
c     ***********************************************************
c     Estado plano de deformacao - elemento quadrilatero bilinear
c     ***********************************************************
       integer nel,ndm,nst,isw,i,j,k,l,m,n
       real*8  e(*),x(ndm,*),u(nst),p(nst),s(nst,nst),xj(ndm,2)
       real*8  det,hx(4),hy(4),xji(ndm,2),Nx(4),Ne(4)
       real*8  xi,eta,ri,si,d11,d12,d21,d22,d33
       real*8  my,nu,a,b,c
c ....................................................... 
       goto(100,200) isw

c	   Leitura das constantes físicas:
c      e(1): Modulo de elasticidade, e(2) poietaon
c      se ept=1, e(3): espeetaura
  100  continue
       read(nin,*) e(1),e(2)
       return
c
  200  continue
c      Matriz constitutiva:
       my = e(1)
       nu = e(2)
       a = 1.d0+nu
       b = a*(1.d0-2.d0*nu)
       c = my*(1.d0-nu)/b
       d11 = c
       d12 = my*nu/b
       d21 = d12
       d22 = c
       d33 = my/(2.d0*a)

c     Matriz de rigidez:
c      zera os coeficientes da matriz
       do i = 1, nst
           do j = 1, nst
               s(i,j)=0.d0
           enddo
       enddo

c      Loop para integracao de Gauus em xi e eta
       do i = 1, 2
        do j = 1, 2
            if (i .eq. 1 .and. j .eq. 1) then
                xi = 0.577350269189626
                eta = -0.577350269189626
            else if (i .eq. 1 .and. j .eq. 2) then
                xi = 0.577350269189626
                eta = 0.577350269189626
            else if (i .eq. 2 .and. j .eq. 1) then
                xi = -0.577350269189626
                eta = 0.577350269189626
            else if (i .eq. 2 .and. j .eq. 2) then
                xi = -0.577350269189626
                eta = -0.577350269189626
            endif
c     Derivadas de Ni em relacao a xi :
            Nx(1)  =   (1.d0+eta) / 4.d0
            Nx(2)  = - (1.d0+eta) / 4.d0
            Nx(3)  = - (1.d0-eta) / 4.d0
            Nx(4)  =   (1.d0-eta) / 4.d0
c     Derivadas de Ni em relacao a eta :
            Ne(1)  =  (1.d0+xi) / 4.d0
            Ne(2)  =  (1.d0-xi) / 4.d0
            Ne(3)  = -(1.d0-xi) / 4.d0
            Ne(4)  = -(1.d0+xi) / 4.d0

c     Matriz Jacobiana
c       Produto escalar linha coluna, dot(N,x)
            do k = 1 , 2
                xj(1,k) = 0.d0
                xj(2,k) = 0.d0
                do l = 1 , 4
                    xj(1,k) = xj(1,k) + Nx(l) * x(k,l)
                    xj(2,k) = xj(2,k) + Ne(l) * x(k,l)
                enddo
            enddo
c
c ... Determinante da matriz Jacobiana:  
c
c ......................................................................              
            det = xj(1,1)*xj(2,2)-xj(2,1)*xj(1,2)
c ......................................................................                    
c
c ... Inversa da matriz Jacobiana:  
c
            xji(1,1) =  xj(2,2) / det
            xji(1,2) = -xj(1,2) / det
            xji(2,1) = -xj(2,1) / det
		    xji(2,2) =  xj(1,1) / det
c
c ... Derivadas das funcoes de interpolacao:
            do k = 1, 4
               hx(k) = xji(1,1)*Nx(k) + xji(1,2)*Ne(k)
               hy(k) = xji(2,1)*Nx(k) + xji(2,2)*Ne(k)
            enddo
            do  m = 1, 4
               k = (m-1)*2+1
               do n = 1, 4
                 l = (n-1)*2+1
          s(l,k)= s(l,k) +( hx(n)*d11*hx(m) + hy(n)*d33*hy(m) ) * det
          s(l,k+1)= s(l,k+1)+( hx(n)*d12*hy(m) + hy(n)*d33*hx(m) )*det
          s(l+1,k)=s(l+1,k)+ ( hy(n)*d12*hx(m) + hx(n)*d33*hy(m) )*det
          s(l+1,k+1)=s(l+1,k+1)+( hy(n)*d22*hy(m) + hx(n)*d33*hx(m))*det
                enddo
            enddo
        enddo
       enddo
c
c	produto  p  =  s u :
c
       call lku(s,u,p,nst)
	   return
1000   continue
       print*, '*** Subrotina elmq01: determinante nulo ou negativo do el
      .emento ',nel
       stop
	   end
