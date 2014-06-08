c     *****************************************************************
c     MEF - Programa exemplo do Metodo dos Elementos Finitos
c     programa para a a disciplina de MEF da COPPE eng. civil
c     Programa original pelo prof. Fernando Ribeiro
c     www.coc.ufrj.br/~fernando/downloads/program.zip
c
c       F·bio CÈsar Canesin <fabio.canesin@gmail.com> 28/05/2014
c     *****************************************************************
      program mef

c     npos -> parametro que controla maximo de memoria
c     default npos = 1073741824 bytes -> 1GB reservado para mef.f
        parameter (npos = 1073741824)

c     Torna o array m (de memÛria) uma "vari·vel global"
        common m(npos)
        common /size/ max

c     Cria um real*8 (equivalente ao "double precision" em x64_86) que
c     servira como ponteiro de numeros de ponto flutuante no programa
        real*8 a(1)

c     Alinha as ponteiros de a e m, fazendo com que eles coincidam
        equivalence (m(1), a(1))

c     Popula a vari·vel max, que tem a posiÁ„o m·xima que pode ocupar
c     uma vari·vel, no caso considerando que s„o vari·veis de 2 bytes ?
        max =  npos/2

c     Executa a subrotina de controle
        call contr(m, a)

c     Interrompe a execusao do programa e finaliza mef
      stop
      end program mef

c     *****************************************************************
c     CONTR -> subrotina de controle do programa que gerencia todo o
c     fluxo de execuÁ„o do programa
c     *****************************************************************
      subroutine contr(m,a)
        common /size/ max
        integer m(*)
        real*8  a(1)
        character*80 fname

c     Abertura de arquivos:
c       nin e nout s„o identificadores dos arquivos, quando utilizamos
c       "open(NUM, fname)" NUM deve ser um numero ˙nico para todos os
c       arquivos abertos no programa e fname um array de characteres
c       que contÈm o fname = "path+nome_do_arquivo"
        nin  = 1
        nout = 2
        print*, 'Arquivo de dados:'
        read(*, '(a)') fname
        open(nin, file=fname)
        print*, 'Arquivo de saida:'
        read(*, '(a)') fname
        open(nout, file=fname)

c      Leitura das vari·veis principais:
c       nnode -> n˙mero de nÛs
c       numel -> n˙mero de elementos
c       numat -> n˙mero de tipos de materiais
c       nen   -> n˙mero m·ximo de nÛs por elemento
c       ndf   -> n˙mero de graus de liberdade por nÛ
c       ndm   -> dimens„o do problema
        read(nin,*) nnode, numel, numat, nen, ndf, ndm

c     Esquema de alocaÁ„o de memÛria:
c         1       i1       i2       i3       i4     i5       i6
c     ---------------------------------------------------------
c     |   e   |   ie   |   ix   |   id   |   x   |   f   |
c     ---------------------------------------------------------
c     Essa È a representaÁ„o da disposiÁ„o no array m das vari·veis
c     onde cada i È um ponteiro para a posiÁ„o que comeÁa cada array,
c     esses s„o na ordem:
c       e: a(1) -> constantes fisicas - m·x 10 por tipo de material
        nen1 = nen+1
c       ie: m(i1) -> Tipo de elemento
        i1 = ( 1 + numat*10 - 1)*2 + 1
c       ix: m(i2) -> conectividade + material do elemento
        i2 =  i1 + numat
c       id: m(i3) -> restriÁıes nodais, 0 = livre ; 1 = restringido
        i3 =  i2 + numel*nen1
c       x: a(i4) -> coordenadas nodais
        i4 = (i3 + nnode*ndf)/2 + 1
c       f: a(i5) -> deslocamento presc. se ID(j,i) = 0 ou forÁa nodal
        i5 =  i4 + nnode*ndm
c       a(i6) -> matrix de c·lculo, numero de nÛs x graus de liberdade
        i6 =  i5 + nnode*ndf

c     Executa a subrotina mem respons·vel por
        call mem(i6)

c      Executa a subrotina de leitura dos dados:
        call rdata(a, m(i1), m(i2), m(i3), a(i4), a(i5), nnode, numel,
     &               numat, nen, ndm, ndf, nin)

c      Executa a subrotina que realiza a numeraÁ„o das equacıes
        call numeq(m(i3), nnode, ndf, neq)

c     Alocatpo de mem=ria:
c
c     i4    i5    i6    i7       i8
c     ------------------------------
c      |  x  |  f  |  u  | jdiag  |
c     ------------------------------

        i7 = (i6 + neq - 1)*2 + 1
        i8 = (i7 + neq)/2 + 1
        call mem(i8)

c      Perfil (skyline) da matriz de rigidez:
        call profil(m(i2),m(i3),m(i7),numel,nen,ndf,neq,ncs)

c      Alocatpo de mem=ria:
c
c         i7     i8   i9     i10    i11    i12   i13   i14
c     -----------------------------------------------
c     | jdiag | am |  xl  |  ul  |  fl  | sl  | ld  |
c     -----------------------------------------------
        nst = nen*ndf
        i9  =  i8  + ncs
        i10 =  i9  + nen*ndm
        i11 =  i10 + nst
        i12 =  i11 + nst
        i13 = (i12 + nst*nst - 1)*2 + 1
        i14 = (i13 + nst)/2 + 1
        call mem(i14)

c     Fortas nodais equivalentes:
        call pload(m(i3),a(i5),a(i6),nnode,ndf)

c     Matriz de rigidez e vetor de fortas corrigido:
        call pform(a,m(i1),m(i2),m(i3),a(i4),a(i5),a(i6),m(i7),a(i8),
     &           a(i9),a(i10),a(i11),a(i12),m(i13),numel,ndm,ndf,nen,
     &           nst,2,.true.,.true.)

c     afl = .true.  corrige o vetor de fortas
c     bfl = .true.  monta a matriz de rigidez global
c     isw = c=digo de instrutpo para a rotina de elemento
c
c     Resolutpo do sistema de equat)es:
        call actcol(a(i8),a(i6),m(i7),neq,.true.,.true.)
c
c    afac = .true.  fatoriza a matriz
c    back = .true.  retrosubstitui
c
c
c    Derivadas:
c
c     call pform(.............)
c
c    Impresspo dos resultados:
c
        call wdata(m(i2),m(i3),a(i4),a(i5),a(i6),nnode,numel,ndm,nen,
     &           ndf,nout)
c
      return
      end

c     *****************************************************************
c     *****************************************************************
      subroutine rdata(e,ie,ix,id,x,f,nnode,numel,numat,nen,ndm,ndf,nin)
        integer nnode, numel, numat, ndm, nen, ndf, ie(*)
        integer ix(nen+1, *), id(ndf, *), iaux(6)
        real*8 e(10, *), x(ndm, *), f(ndf, *), aux(6), dum(1)

        do 100 i = 1, numat
            read(nin,*) ma,iel
            call elmlib(e(1,ma),dum,dum,dum,dum,1,iel,1,1,1,nin)
            ie(ma) = iel
 100    continue
        do 200 i = 1, nnode
            read(nin,*) k, (x(j,k),j=1,ndm)
 200    continue
        do 300 i = 1, numel
            read(nin,*) k, (ix(j,k),j=1,nen+1)
 300    continue
 400    continue
        read (nin,*) k, (iaux(j), j = 1, ndf)
        if (k .le. 0) goto 500
        do 410 j = 1, ndf
            id(j,k) = iaux(j)
 410    continue
        goto 400
 500    continue
        read (nin,*) k, (aux(j), j = 1, ndf)
        if (k .le. 0) goto 600
        do 510 j = 1, ndf
            f(j,k) = aux(j)
 510    continue
        goto 500
 600    continue
      return
      end

c     *****************************************************************
c     *****************************************************************
      subroutine numeq(id,nnode,ndf,neq)
        integer nnode,ndf,neq,id(ndf,*)

        neq = 0
        do 100 n = 1, nnode
        do 100 i = 1, ndf
            j = id(i,n)
            if (j .eq. 0) then
                neq = neq + 1
                id(i,n) = neq
            else
                id(i,n) = 0
            endif
 100    continue
      return
      end

c     *****************************************************************
c     *****************************************************************
      subroutine profil(ix,id,jdiag,numel,nen,ndf,neq,ncs)
        integer numel,nen,ndf,neq,ncs,ix(nen+1,*),id(ndf,*),jdiag(*)

        call mzero (jdiag,1,neq)

c     alturas de coluna:
        do 400 nel = 1, numel
        do 400 i = 1, nen
            noi = ix(i,nel)
            if (noi .eq. 0) goto 400
            do 300 k = 1, ndf
                kk = id(k, noi)
                if (kk .eq. 0) goto 300
                do 200 j = i, nen
                    noj = ix(j, nel)
                    if (noj .eq. 0) goto 200
                    do 100 l = 1, ndf
                        ll = id(l, noj)
                        if (ll .eq. 0) goto 100
                        m = max0(kk, ll)
                        jdiag(m) = max0(jdiag(m), iabs(kk-ll))
 100                continue
 200            continue
 300        continue
 400    continue

c     ponteiros da diagonal:
        ncs = 1
        jdiag(1) = 1
        if (neq .eq. 1) return
        do 500 i = 2, neq
            jdiag(i) = jdiag(i) + jdiag(i-1) + 1
 500    continue
        ncs = jdiag(neq)
      return
      end

c     *****************************************************************
c     *****************************************************************
      subroutine pload (id,f,u,nnode,ndf)
        integer nnode, ndf, id(ndf,*)
        real*8 f(ndf,*), u(*)

        do 100 i = 1, nnode
        do 100 j = 1, ndf
            k = id(j,i)
            if (k .gt. 0) u(k) = f(j,i)
 100    continue
      return
      end

c     *****************************************************************
c     *****************************************************************
      subroutine pform(e,ie,ix,id,x,f,u,jdiag,am,xl,ul,fl,sl,ld,numel,
     &                 ndm,ndf,nen,nst,isw,afl,bfl)
        integer numel,ndm,ndf,nen,nst,isw
        integer ie(*),ix(nen+1,*),id(ndf,*),jdiag(*),ld(ndf,*)
        real*8  e(10, *), x(ndm, *), f(ndf, *), u(*), am(*),
     &          xl(ndm, *), ul(ndf, *), fl(*), sl(*)
        logical afl, bfl

c    loop nos elementos:
c
        do 700 nel = 1, numel

c        loop nos n=s do elemento:
        do 600 i = 1, nen
            no = ix(i,nel)
            if (no .gt. 0) goto 300

c            zera os vetores locais se no  0:
            do 100 j = 1, ndm
                xl(j,i) = 0.d0
 100        continue
            do 200 j = 1, ndf
                ul(j,i) = 0.d0
                ld(j,i) = 0
 200        continue
            goto 600

c            forma os vetores locais se no > 0:
 300        continue
            do 400 j = 1, ndm
                xl(j,i) = x(j,no)
 400        continue
            do 500 j = 1, ndf
c                numeratpo global das equat)es do elemento:
                k = id(j,no)
                ld(j,i) = k
c                deslocamentos prescritos:
                ul(j,i) = 0.d0
                if (k .eq. 0) ul(j,i) = f(j,no)
 500        continue
 600    continue
        ma  = ix(nen+1,nel)
        iel = ie(ma)

c        biblioteca de elementos:
        call elmlib(e(1,ma),xl,ul,fl,sl,nel,iel,ndm,nst,isw,1)

c        vetores globais:
        call addstf(am,u,jdiag,fl,sl,ld,nst,afl,bfl)
 700    continue
      return
      end

c     *****************************************************************
c     *****************************************************************
      subroutine addstf(a,b,jdiag,p,s,ld,nst,afl,bfl)
        integer jdiag(*),ld(*),nst
        real*8 a(*),b(*),p(*),s(nst,*)
        logical afl,bfl

        do 200 j = 1, nst
            k = ld(j)
            if (k .eq. 0) goto 200

c        corretpo do vetor de fortas:
            if (afl) b(k) = b(k) - p(j)
            if (.not. bfl)  goto 200

c        matriz de rigidez:
            l = jdiag(k) - k
            do 100 i = 1, nst
                m = ld(i)
                if (m .gt. k .or. m .eq. 0) goto 100
                m = l + m
                a(m) = a(m) + s(i,j)
 100        continue
 200    continue
      return
      end

c     *****************************************************************
c     *****************************************************************
      subroutine elmlib(e,xl,ul,fl,sl,nel,iel,ndm,nst,isw,nin)
        integer nel,iel,ndm,nst,isw
        real*8 e(*),xl(*),ul(*),fl(*),sl(*)

        goto (100,200) iel
        print*, ' TIPO DE ELEMENTO INEXISTENTE, ',iel,'ELEMENTO: ',nel
        stop
 100    call elmt01(e,xl,ul,fl,sl,nel,ndm,nst,isw,nin)
      return
 200    call elmt02(e,xl,ul,fl,sl,nel,ndm,nst,isw,nin)
      return
c  300    call elmt01(e,xl,ul,fl,sl,nel,ndm,nst,isw,nin)
      return
      end

c     *****************************************************************
c     Estado plano de deformacao - elemento triangular linear
c     *****************************************************************
      subroutine elmt01(e,x,u,p,s,nel,ndm,nst,isw,nin)
        integer nel,ndm,nst,isw
        real*8  e(*),x(ndm,*),u(nst),p(nst),s(nst,nst)
        real*8  det,xj11,xj12,xj21,xj22,hx(3),hy(3)
        real*8  xji11,xji12,xji21,xji22,d11,d12,d21,d22,d33
        real*8  my,nu,a,b,c

        goto(100,200) isw

c    leitura das constantes ffsicas:
 100    continue
        read(nin,*) e(1),e(2)
      return

c    matriz de rigidez:
 200    continue

c     Matriz Jacobiana:
        xj11 = x(1,1)-x(1,3)
        xj12 = x(2,1)-x(2,3)
        xj21 = x(1,2)-x(1,3)
        xj22 = x(2,2)-x(2,3)
        det  = xj11*xj22-xj12*xj21
        if (det .le. 0.d0) goto 1000
c
c    Inversa da matriz Jacobiana:
c
        xji11 =  xj22/det
        xji12 = -xj12/det
        xji21 = -xj21/det
        xji22 =  xj11/det
c
c     Derivadas das funcoes de interpolacao:
c
        hx(1) =  xji11
        hx(2) =  xji12
        hx(3) = -xji11-xji12
        hy(1) =  xji21
        hy(2) =  xji22
        hy(3) = -xji21-xji22
c
c     Matriz constitutiva:
c
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
c
c     Matriz de rigidez:
c
        wt = 0.5d0*det
        do 220 j = 1, 3
            k = (j-1)*2+1
            do 210 i = 1, 3
                l = (i-1)*2+1
                s(l,k)   = ( hx(i)*d11*hx(j) + hy(i)*d33*hy(j) ) * wt
                s(l,k+1) = ( hx(i)*d12*hy(j) + hy(i)*d33*hx(j) ) * wt
                s(l+1,k) = ( hy(i)*d21*hx(j) + hx(i)*d33*hy(j) ) * wt
                s(l+1,k+1) = ( hy(i)*d22*hy(j) + hx(i)*d33*hx(j) ) * wt

 210        continue
 220    continue

c    produto  p  =  s u :
        call lku(s,u,p,nst)
      return
1000    continue
        print*, '*** Subrotina ELMT01: determinante nulo ou negativo
     & do elemento ',nel
      stop
      end

c     *****************************************************************
c     ELMT02 - Estado plano de tensao - elemento triangular linear
c     *****************************************************************
      subroutine elmt02(e,x,u,p,s,nel,ndm,nst,isw,nin)
        integer nel,ndm,nst,isw
        real*8  e(*),x(ndm,*),u(nst),p(nst),s(nst,nst)
        real*8  det,xj11,xj12,xj21,xj22,hx(3),hy(3)
        real*8  xji11,xji12,xji21,xji22,d11,d12,d21,d22,d33
        real*8  my,nu,thic,a

        goto(100,200) isw

c    leitura das constantes ffsicas:
 100    continue
        read(nin,*) e(1),e(2),e(3)
      return

c    matriz de rigidez:
 200    continue

c     Matriz Jacobiana:
        xj11 = x(1,1)-x(1,3)
        xj12 = x(2,1)-x(2,3)
        xj21 = x(1,2)-x(1,3)
        xj22 = x(2,2)-x(2,3)
        det  = xj11*xj22-xj12*xj21
        if (det .le. 0.d0) goto 1000
c
c    Inversa da matriz Jacobiana:
c
        xji11 =  xj22/det
        xji12 = -xj12/det
        xji21 = -xj21/det
        xji22 =  xj11/det
c
c     Derivadas das funcoes de interpolacao:
c
        hx(1) =  xji11
        hx(2) =  xji12
        hx(3) = -xji11-xji12
        hy(1) =  xji21
        hy(2) =  xji22
        hy(3) = -xji21-xji22
c
c     Matriz constitutiva:
c
        my   = e(1)
        nu   = e(2)
        thic = e(3)
        a = my/(1.d0-nu*nu)
        d11 = a
        d12 = a*nu
        d21 = d12
        d22 = a
        d33 = my/(2.d0*(1.d0+nu))
c
c     Matriz de rigidez:
c
        wt = 0.5d0*det*thic
        do 220 j = 1, 3
            k = (j-1)*2+1
            do 210 i = 1, 3
                l = (i-1)*2+1
                s(l,k)     = ( hx(i)*d11*hx(j) + hy(i)*d33*hy(j) ) * wt
                s(l,k+1)   = ( hx(i)*d12*hy(j) + hy(i)*d33*hx(j) ) * wt
                s(l+1,k)   = ( hy(i)*d21*hx(j) + hx(i)*d33*hy(j) ) * wt
                s(l+1,k+1) = ( hy(i)*d22*hy(j) + hx(i)*d33*hx(j) ) * wt
 210        continue
 220    continue

c    produto  p  =  s u :
        call lku(s,u,p,nst)
      return
1000    continue
        print*, '*** Subrotina ELMT02: determinante nulo ou negativo do
     & elemento ',nel
      stop
      end

c     ****************************************************************
c       Programa para resolucao de sistemas de equacoes algebricas
c       lineares por eliminacao de gauss com decomposicao LtDL, vﬂlido
c       somente para matrizes simetricas. Armazenamento skyline.
c
c       - Parametros de entrada :
c
c        a       = Coeficientes da matriz, armazenados por altura
c               efetiva de coluna.
c        b      = Vetor independente.
c       jdiag = Vetor apontador  do armazenamento.
c       neq   = Numero de equacoes.
c       afac  = Flag. Fatoriza a matriz se afac =.true.
c       back = Flag. Retrosubstitui se back = .true.
c
c       - Parametros de saida :
c
c       a     = Coeficientes da matriz triangularizada.
c       b     = Vetor solucao.
c       jdiag = Inalterado.
c       neq   = Inalterado.
c       afac  = Inalterado.
c       back  = Inalterado.
c     ****************************************************************
      subroutine actcol(a, b, jdiag, neq, afac, back)
        implicit real*8 (a-h, o-z)
        common/engys/ aengy
        dimension a(*), b(*), jdiag(*)
        logical afac, back

c.... Fatorizatpo da matriz e redutpo do vetor independente:
        aengy = 0.0d0
        jr = 0
        do 600 j = 1,neq
            jd = jdiag(j)
            jh = jd - jr
            is = j - jh + 2
            if(jh-2) 600,300,100
 100        if(.not.afac) goto 500
            ie = j - 1
            k = jr + 2
            id = jdiag(is - 1)

c....        Reduz as equat)es, exceto os termos da diagonal:
            do 200 i = is, ie
                ir = id
                id = jdiag(i)
                ih = min0(id-ir-1,i-is+1)
                if(ih.gt.0) a(k) = a(k) - dot(a(k-ih),a(id-ih),ih)
 200        k = k + 1

c....        Reduz os termos da diagonal:
 300        if(.not.afac) goto 500
            ir = jr+1
            ie = jd - 1
            k = j - jd
            do 400 i = ir, ie
                id = jdiag(k+i)
                if(a(id).eq.0.0d0) goto 400
                d = a(i)
                a(i) = a(i)/a(id)
                a(jd) = a(jd) - d*a(i)
 400        continue

c....        Reduz o vetor independente:
 500          if(back) b(j) = b(j) - dot(a(jr+1),b(is-1),jh-1)
 600        jr = jd
            if(.not.back) return

c.... Divide pelos pivots da diagonal:
            do 700 i = 1,neq
                id = jdiag(i)
                if(a(id).ne.0.0d0) b(i) = b(i)/a(id)
 700            aengy = aengy + b(i)*b(i)*a(id)

c.... Retrosubstitui:
            j = neq
            jd = jdiag(j)
 800    d = b(j)
        j = j - 1
        if(j.le.0) return
        jr = jdiag(j)
        if(jd-jr.le.1) goto 1000
        is = j - jd + jr + 2
        k = jr - is + 1
        do 900 i = is,j
 900    b(i) = b(i) - a(i+k)*d
1000    jd = jr
        goto 800
      end

c     ****************************************************************
c     ****************************************************************
      subroutine lku(s,u,p,nst)
        integer nst
        real*8 s(nst,nst),u(nst),p(nst)
        do 200 i = 1, nst
            p(i) = 0.d0
            do 100 j = 1, nst
                p(i) = p(i) + s(i,j) * u(j)
 100        continue
 200    continue
      return
      end

c     ****************************************************************
c     ****************************************************************
      real*8 function dot(a,b,n)
        integer n
        real*8 a(*), b(*)

        dot = 0.d0
        do 100 i = 1, n
            dot = dot + a(i) * b(i)
 100    continue
      return
      end

c     ****************************************************************
c     ****************************************************************
      subroutine mem(npos)
        common /size/ max
        integer npos
        if ( (npos-1) .gt. max ) then
            print*, ' Mem=ria insuficiente !'
        stop
        endif
      return
      end

c     ****************************************************************
c     ****************************************************************
      subroutine mzero(m,i1,i2)
        integer i1,i2,m(*)
        do 100 i = 1, i2-i1 + 1
            m(i) = 0
 100    continue
      return
      end

c     ****************************************************************
c     ****************************************************************
      subroutine  azero(a,i1,i2)
        integer i1,i2
        real*8 a(*)
        do 100 i = 1,  i2 - i1 + 1
            a(i) = 0.d0
 100    continue
      return
      end

c     ****************************************************************
c     ****************************************************************
      subroutine wdata(ix,id,x,f,u,nnode,numel,ndm,nen,ndf,nout)
        integer nnode,numel,ndm,nen,ndf,ix(nen+1,*),id(ndf,*)
        real*8  x(ndm,*),f(ndf,*),u(*),aux(6)
        write(nout,'(a,i5)') 'coor ', nnode
        do 100 i = 1, nnode
            write(nout,'(i10,3e15.5)') i,(x(j,i),j=1,2),0.0
 100    continue
        write(nout,'(a,i10)') 'elem ', numel
        do 200 i = 1, numel
            write(nout,'(10i10)') i,3,(ix(j,i),j=1,nen)
 200    continue
        write(nout,'(a,i2)') 'nosc ',ndf
        do 300 i = 1, nnode
            do 310 j = 1, ndf
                aux(j) = f(j,i)
                k = id(j,i)
                if(k .gt. 0) aux(j) = u(k)
 310        continue
            write (nout,'(i10,6e15.5e3)') i,(aux(j),j=1,ndf)
 300    continue
        write(nout,'(a,i2)') 'nvec '
        do 400 i = 1, nnode
            do 410 j = 1, ndf
                aux(j) = f(j,i)
                k = id(j,i)
                if(k .gt. 0) aux(j) = u(k)
 410        continue
            write (nout,'(i10,6e15.5e3)') i,(aux(j),j=1,ndf),0.0
 400    continue
        write (nout,'(a)') 'end '
      return
      end
