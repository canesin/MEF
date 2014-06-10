C     *****************************************************************
C     MEF - Programa exemplo do Metodo dos Elementos Finitos
C     programa para a a disciplina de MEF da COPPE eng. civil
C     Programa original pelo prof. Fernando Ribeiro
C     www.coc.ufrj.br/~fernando/downloads/program.zip
C
C       Fabio Cesar Canesin <fabio.canesin@gmail.com> 28/05/2014
C     *****************************************************************
      program mef

C     npos -> parametro que controla maximo de memoria
        parameter (npos = 268435456)

C     Torna o array m (de memoria) uma "variavel global"
        common m(npos)
        common /size/ max

C     Cria um real*8 (equivalente ao "double precision" em x64_86) que
C     servira como ponteiro de numeros de ponto flutuante no programa
        real*8 a(1)

C     Alinha as ponteiros de a e m, fazendo com que eles coincidam
        equivalence (m(1), a(1))

C     Popula a variavel max, que tem a posicao maxima que pode ocupar
C     uma variavel, no caso considerando que sao variaveis de 2 bytes ?
        max =  npos/2

C     Executa a subrotina de controle
        call contr(m, a)

C     Interrompe a execusao do programa e finaliza mef
      stop
      end program mef

C     *****************************************************************
C     CONTR -> subrotina de controle do programa que gerencia todo o
C     fluxo de execucao do programa
C     *****************************************************************
      subroutine contr(m,a)
        common /size/ max
        integer m(*)
        real*8  a(1)
        character*80 fname

C     Abertura de arquivos:
C       nin e nout sao identificadores dos arquivos, quando utilizamos
C       "open(NUM, fname)" NUM deve ser um numero unico para todos os
C       arquivos abertos no programa e fname um array de characteres
C       que contem o fname = "path+nome_do_arquivo"
        nin  = 1
        nout = 2
        print*, 'Arquivo de dados:'
        read(*, '(a)') fname
        open(nin, file=fname)
        print*, 'Arquivo de saida:'
        read(*, '(a)') fname
        open(nout, file=fname)

C     Leitura das variaveis principais:
C       O programa le em ordem a primeira linha como:
C       nnode -> numero de nos
C       numel -> numero de elementos
C       numat -> numero de tipos de materiais
C       nen   -> numero maximo de nos por elemento
C       ndf   -> numero de graus de liberdade por no
C       ndm   -> dimensao do problema
        read(nin,*) nnode, numel, numat, nen, ndf, ndm

C     Esquema de alocacao de memoria:
C         1       i1       i2       i3       i4     i5       i6
C     ---------------------------------------------------------
C     |   e   |   ie   |   ix   |   id   |   x   |   f   |
C     ---------------------------------------------------------
C     Essa e a representacao da disposicao no array m das variaveis
C     onde cada i e um ponteiro para a posicao que comeca cada array,
C     esses sao na ordem:
C       e: a(1) -> constantes fisicas - max 10 por tipo de material
        nen1 = nen+1
C       ie: m(i1) -> Tipo de elemento
        i1 = ( 1 + numat*10 - 1)*2 + 1
C       ix: m(i2) -> conectividade + material do elemento
        i2 =  i1 + numat
C       id: m(i3) -> restricoes nodais, 0 = livre ; 1 = restringido
        i3 =  i2 + numel*nen1
C       x: a(i4) -> coordenadas nodais
        i4 = (i3 + nnode*ndf)/2 + 1
C       f: a(i5) -> deslocamento presc. se ID(j,i) = 0 ou forca nodal
        i5 =  i4 + nnode*ndm
C       a(i6) -> vetor incognitas, numero de nos x graus de liberdade
        i6 =  i5 + nnode*ndf

C     Executa a subrotina mem responsavel por verificar disponibilidade
C     de memoria, sempre chamada antes de alguma operacao que vai alocar
C     mais posicoes em m

        call mem(i6)

C     Executa a subrotina de leitura dos dados:
C       rdata(e: a(1) -> constantes fisicas,
C             ie: m(i1) -> Tipo de elemento,
C             ix: m(i2) -> conectividade + material do elemento,
C             id: m(i3) -> restricoes nodais,
C             x: a(i4) -> coordenadas nodais,
C             f: a(i5) -> deslocamento prescitos,
C             nnode: numero de nos,
C             numel: numero de elementos,
C             numat: numero de tipos de materiais,
C             nen:   numero maximo de nos por elemento,
C             ndm:   numero de dimensoes do problema,
C             ndf:   numero de graus de liberdade por no,
C             nin:   identificador do arquivo de entrada de dados)
        call rdata(a, m(i1), m(i2), m(i3), a(i4), a(i5), nnode, numel,
     &               numat, nen, ndm, ndf, nin)

C     Executa a subrotina que realiza a numeracao das equacoes:
C       rdata(id: m(i3) -> restricoes nodais,
C             nnode: numero de nos,
C             ndf:   numero de graus de liberdade por no,
C             neq:   numero de equacoes)
        call numeq(m(i3), nnode, ndf, neq)

C     Alocacao de memoria:
C        i4    i5    i6      i7      i8
C      -----------------------------
C      |  x  |  f  |  u  |  jdiag  |
C      -----------------------------
C       jdiag: m(i7) -> vetor de posicoes dos elementos diagonais da
C                       matriz de rigidez no format skyline
        i7 = (i6 + neq - 1)*2 + 1
        i8 = (i7 + neq)/2 + 1
        call mem(i8)

C     Gerar perfil skyline da matriz de rigidez:
C       profil(ix: m(i2) -> conectividade + material do elemento,
C              id: m(i3) -> restricoes nodais,
C              jdiag: m(i7) -> vetor de posicoes dos elementos diagonais
C                              da matriz de rigidez no format skyline
C              numel: numero de elementos,
C              nen:   numero maximo de nos por elemento,
C              ndf:   numero de graus de liberdade por no,
C              neq:   numero de equacoes,
C              ncs:   ???)
        call profil(m(i2),m(i3),m(i7),numel,nen,ndf,neq,ncs)

C     Alocacao de memoria:
C         i7    i8    i9    i10     i11   i12   i13   i14
C     -----------------------------------------------
C     | jdiag | am |  xl  |  ul  |  fl  | sl  | ld  |
C     -----------------------------------------------
        nst = nen*ndf
        i9  =  i8  + ncs
        i10 =  i9  + nen*ndm
        i11 =  i10 + nst
        i12 =  i11 + nst
        i13 = (i12 + nst*nst - 1)*2 + 1
        i14 = (i13 + nst)/2 + 1
        call mem(i14)

C     Calculo das forcas nodais equivalentes:
        call pload(m(i3),a(i5),a(i6),nnode,ndf)

C     Matriz de rigidez e vetor de fortas corrigido:
        call pform(a,m(i1),m(i2),m(i3),a(i4),a(i5),a(i6),m(i7),a(i8),
     &           a(i9),a(i10),a(i11),a(i12),m(i13),numel,ndm,ndf,nen,
     &           nst,2,.true.,.true.)

C     afl = .true.  corrige o vetor de fortas
C     bfl = .true.  monta a matriz de rigidez global
C     isw = c=digo de instrutpo para a rotina de elemento
c
C     Resolucao do sistema de equacoes:
        call actcol(a(i8),a(i6),m(i7),neq,.true.,.true.)
c
C     afac = .true.  fatoriza a matriz
C     back = .true.  retrosubstitui
c
C     Derivadas:
C     call pform(.............)
c
C     Impressao dos resultados:
        call wdata(m(i2),m(i3),a(i4),a(i5),a(i6),nnode,numel,ndm,nen,
     &           ndf,nout)
c
      return
      end

C     *****************************************************************
C     Subrotina para leitura de dados do arquivo .dat
C     *****************************************************************
      subroutine rdata(e,ie,ix,id,x,f,nnode,numel,numat,nen,ndm,ndf,nin)
        integer nnode, numel, numat, ndm, nen, ndf, ie(*)
        integer ix(nen+1, *), id(ndf, *), iaux(6)
        real*8 e(10, *), x(ndm, *), f(ndf, *), aux(6), dum(1)

C       Loop de 1 a numat (numero de materiais, 1 linha do .dat)
        do 100 i = 1, numat
C           Segunda linha .dat: ma -> material, iel -> tipo de elemento
            read(nin,*) ma,iel
C           Biblioteca de elem. para ler propriedades materiais:
C              elmlib(e:  constantes fisicas,
C                     xl: coordenadas locais ??,
C                     ul: vetor incognitas local ??,
C                     fl: vetor forcas local ??,
C                     sl:,
C                     nel: numero do elemento ??,
C                     iel: tipo do elemento,
C                     ndm: numero de dimensoes,
C                     nst: ,
C                     isw: se 1 ler propriedades, se 2 calculo de Kij,
C                     nin: id arquivo de entrada)
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

C     *****************************************************************
C     Subrotina que numera as equacoes
C     *****************************************************************
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

C     *****************************************************************
C     Subrotina que cria o perfil Skyline da matriz de rigidez
C     *****************************************************************
      subroutine profil(ix,id,jdiag,numel,nen,ndf,neq,ncs)
        integer numel,nen,ndf,neq,ncs,ix(nen+1,*),id(ndf,*),jdiag(*)

        call mzero (jdiag,1,neq)

C     alturas de coluna:
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

C     Ponteiros da diagonal:
        ncs = 1
        jdiag(1) = 1
        if (neq .eq. 1) return
        do 500 i = 2, neq
            jdiag(i) = jdiag(i) + jdiag(i-1) + 1
 500    continue
        ncs = jdiag(neq)
      return
      end

C     *****************************************************************
C     *****************************************************************
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

C     *****************************************************************
C     *****************************************************************
      subroutine pform(e,ie,ix,id,x,f,u,jdiag,am,xl,ul,fl,sl,ld,numel,
     &                 ndm,ndf,nen,nst,isw,afl,bfl)
        integer numel,ndm,ndf,nen,nst,isw
        integer ie(*),ix(nen+1,*),id(ndf,*),jdiag(*),ld(ndf,*)
        real*8  e(10, *), x(ndm, *), f(ndf, *), u(*), am(*),
     &          xl(ndm, *), ul(ndf, *), fl(*), sl(*)
        logical afl, bfl

C     loop nos elementos:
        do 700 nel = 1, numel

C       loop nos do elemento:
        do 600 i = 1, nen
            no = ix(i,nel)
            if (no .gt. 0) goto 300

C           Zera os vetores locais se no 0:
            do 100 j = 1, ndm
                xl(j,i) = 0.d0
 100        continue
            do 200 j = 1, ndf
                ul(j,i) = 0.d0
                ld(j,i) = 0
 200        continue
            goto 600

C            forma os vetores locais se no > 0:
 300        continue
            do 400 j = 1, ndm
                xl(j,i) = x(j,no)
 400        continue
            do 500 j = 1, ndf
C                numeratpo global das equat)es do elemento:
                k = id(j,no)
                ld(j,i) = k
C                deslocamentos prescritos:
                ul(j,i) = 0.d0
                if (k .eq. 0) ul(j,i) = f(j,no)
 500        continue
 600    continue
        ma  = ix(nen+1,nel)
        iel = ie(ma)

C        biblioteca de elementos:
        call elmlib(e(1,ma),xl,ul,fl,sl,nel,iel,ndm,nst,isw,1)

C        vetores globais:
        call addstf(am,u,jdiag,fl,sl,ld,nst,afl,bfl)
 700    continue
      return
      end

C     *****************************************************************
C     *****************************************************************
      subroutine addstf(a,b,jdiag,p,s,ld,nst,afl,bfl)
        integer jdiag(*),ld(*),nst
        real*8 a(*),b(*),p(*),s(nst,*)
        logical afl,bfl

        do 200 j = 1, nst
            k = ld(j)
            if (k .eq. 0) goto 200

C           Correcao do vetor de forcas:
            if (afl) b(k) = b(k) - p(j)
            if (.not. bfl)  goto 200

C           Matriz de rigidez:
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

C     *****************************************************************
C     Subrotina emlib - biblioteca de elementos
C     *****************************************************************
      subroutine elmlib(e,xl,ul,fl,sl,nel,iel,ndm,nst,isw,nin)
        integer nel,iel,ndm,nst,isw
        real*8 e(*),xl(*),ul(*),fl(*),sl(*)

        goto (100,200,300,400) iel
        print*, ' TIPO DE ELEMENTO INEXISTENTE: ',iel,'ELEMENTO: ',nel
        stop
 100    call elmt01(e,xl,ul,fl,sl,nel,ndm,nst,isw,nin)
      return
 200    call elmt02(e,xl,ul,fl,sl,nel,ndm,nst,isw,nin)
      return
 300    call elmq01(e,xl,ul,fl,sl,nel,ndm,nst,isw,nin,0)
      return
 400    call elmq01(e,xl,ul,fl,sl,nel,ndm,nst,isw,nin,1)
      return
      end

C     *****************************************************************
C     Estado plano de deformacao - elemento triangular linear
C     *****************************************************************
      subroutine elmt01(e,x,u,p,s,nel,ndm,nst,isw,nin)
        integer nel,ndm,nst,isw
        real*8  e(*),x(ndm,*),u(nst),p(nst),s(nst,nst)
        real*8  det,xj11,xj12,xj21,xj22,hx(3),hy(3)
        real*8  xji11,xji12,xji21,xji22,d11,d12,d21,d22,d33
        real*8  my,nu,a,b,c

        goto(100,200) isw

C     leitura das constantes ffsicas:
 100    continue
        read(nin,*) e(1),e(2)
      return

C     Matriz de rigidez:
 200    continue

C     Matriz Jacobiana:
        xj11 = x(1,1)-x(1,3)
        xj12 = x(2,1)-x(2,3)
        xj21 = x(1,2)-x(1,3)
        xj22 = x(2,2)-x(2,3)
        det  = xj11*xj22-xj12*xj21
        if (det .le. 0.d0) goto 1000

C     Inversa da matriz Jacobiana:
        xji11 =  xj22/det
        xji12 = -xj12/det
        xji21 = -xj21/det
        xji22 =  xj11/det

C     Derivadas das funcoes de interpolacao:
        hx(1) =  xji11
        hx(2) =  xji12
        hx(3) = -xji11-xji12
        hy(1) =  xji21
        hy(2) =  xji22
        hy(3) = -xji21-xji22

C     Matriz constitutiva:
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

C     Matriz de rigidez:
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

C       Produto  p = S u :
        call lku(s,u,p,nst)
      return
1000    continue
        print*, "*** Subrotina ELMT01: determinante nulo ou negativo
     & do elemento ",nel
      stop
      end

C     *****************************************************************
C     ELMT02 - Estado plano de tensao - elemento triangular linear
C     *****************************************************************
      subroutine elmt02(e,x,u,p,s,nel,ndm,nst,isw,nin)
        integer nel,ndm,nst,isw
        real*8  e(*),x(ndm,*),u(nst),p(nst),s(nst,nst)
        real*8  det,xj11,xj12,xj21,xj22,hx(3),hy(3)
        real*8  xji11,xji12,xji21,xji22,d11,d12,d21,d22,d33
        real*8  my,nu,thic,a

        goto(100,200) isw

C     leitura das constantes fisicas:
 100    continue
        read(nin,*) e(1),e(2),e(3)
      return

C     Matriz de rigidez:
 200    continue

C     Matriz Jacobiana:
        xj11 = x(1,1)-x(1,3)
        xj12 = x(2,1)-x(2,3)
        xj21 = x(1,2)-x(1,3)
        xj22 = x(2,2)-x(2,3)
        det  = xj11*xj22-xj12*xj21
        if (det .le. 0.d0) goto 1000

C     Inversa da matriz Jacobiana:
        xji11 =  xj22/det
        xji12 = -xj12/det
        xji21 = -xj21/det
        xji22 =  xj11/det

C     Derivadas das funcoes de interpolacao:
        hx(1) =  xji11
        hx(2) =  xji12
        hx(3) = -xji11-xji12
        hy(1) =  xji21
        hy(2) =  xji22
        hy(3) = -xji21-xji22

C     Matriz constitutiva:
        my   = e(1)
        nu   = e(2)
        thic = e(3)
        a = my/(1.d0-nu*nu)
        d11 = a
        d12 = a*nu
        d21 = d12
        d22 = a
        d33 = my/(2.d0*(1.d0+nu))

C     Matriz de rigidez:
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

C       Produto  p = S u :
        call lku(s,u,p,nst)
      return
1000    continue
        print*, '*** Subrotina ELMT02: determinante nulo ou negativo do
     & elemento ',nel
      stop
      end

      subroutine elmq01(e,x,u,p,s,nel,ndm,nst,isw,nin,ept)
C     ***********************************************************
C     Elemento quadrilatero bilinear, EPD e EPT
C     ***********************************************************
       integer nel,ndm,nst,isw,i,j,k,l,m,n,ept
       real*8  e(*),x(ndm,*),u(nst),p(nst),s(nst,nst),xj(ndm,2)
       real*8  det,hx(4),hy(4),xji(ndm,2),Nx(4),Ne(4)
       real*8  xi,eta,ri,si,d11,d12,d21,d22,d33
       real*8  my,nu,lam,lamb,mu

       goto(100,200) isw

C     Leitura das constantes fiicas:
C      e(1): Modulo de elasticidade, e(2) poisson
C      se ept=1, e(3): espessura
  100  continue
       if (ept .eq. 1) then
          read(nin,*) e(1),e(2),e(3)
       else
          read(nin,*) e(1),e(2)
       endif
       return

  200  continue
C      Matriz constitutiva:
       my = e(1)
       nu = e(2)
       lam = nu*my/((1.d0+nu)*(1-2.d0*nu))
       mu = my/(2d0*(1.d0+nu))
       if (ept .eq. 1) then
          lamb = 2.d0*(lam*mu)/(lam+2.d0*mu)
          lam = lmab
       else
          continue
       endif

       d11 = lam + 2d0*mu
       d12 = lam
       d33 = mu
       d21 = d12
       d22 = d11

C      Matriz de rigidez:
C      Zera os coeficientes da matriz:
       do i = 1, nst
           do j = 1, nst
               s(i,j)=0.d0
           enddo
       enddo

C      Loop para integracao de Gauus em xi e eta
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
C     Derivadas de Ni em relacao a xi :
            Nx(1)  =   (1.d0+eta) / 4.d0
            Nx(2)  = - (1.d0+eta) / 4.d0
            Nx(3)  = - (1.d0-eta) / 4.d0
            Nx(4)  =   (1.d0-eta) / 4.d0
C     Derivadas de Ni em relacao a eta :
            Ne(1)  =  (1.d0+xi) / 4.d0
            Ne(2)  =  (1.d0-xi) / 4.d0
            Ne(3)  = -(1.d0-xi) / 4.d0
            Ne(4)  = -(1.d0+xi) / 4.d0

C     Matriz Jacobiana
C       Produto escalar linha coluna, dot(N,x)
            do k = 1 , 2
                xj(1,k) = 0.d0
                xj(2,k) = 0.d0
                do l = 1 , 4
                    xj(1,k) = xj(1,k) + Nx(l) * x(k,l)
                    xj(2,k) = xj(2,k) + Ne(l) * x(k,l)
                enddo
            enddo
c
C     Determinante da matriz Jacobiana:  
            det = xj(1,1)*xj(2,2)-xj(2,1)*xj(1,2)
c
C     Inversa da matriz Jacobiana:  
c
            xji(1,1) =  xj(2,2) / det
            xji(1,2) = -xj(1,2) / det
            xji(2,1) = -xj(2,1) / det
		    xji(2,2) =  xj(1,1) / det
c
C     Derivadas das funcoes de interpolacao:
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
       print*, '*** Subrotina elmq01: determinante nulo ou negativo do
     & elemento ',nel
       stop
	   end
     
C     ****************************************************************
C       Programa para resolucao de sistemas de equacoes algebricas
C       lineares por eliminacao de gauss com decomposicao LtDL, valido
C       somente para matrizes simetricas. Armazenamento skyline.
c
C       - Parametros de entrada :
c
C        a       = Coeficientes da matriz, armazenados por altura
C               efetiva de coluna.
C        b      = Vetor independente.
C       jdiag = Vetor apontador  do armazenamento.
C       neq   = Numero de equacoes.
C       afac  = Flag. Fatoriza a matriz se afac =.true.
C       back = Flag. Retrosubstitui se back = .true.
c
C       - Parametros de saida :
c
C       a     = Coeficientes da matriz triangularizada.
C       b     = Vetor solucao.
C       jdiag = Inalterado.
C       neq   = Inalterado.
C       afac  = Inalterado.
C       back  = Inalterado.
C     ****************************************************************
      subroutine actcol(a, b, jdiag, neq, afac, back)
        implicit real*8 (a-h, o-z)
        common/engys/ aengy
        dimension a(*), b(*), jdiag(*)
        logical afac, back

c.... Fatorizacao da matriz e reducao do vetor independente:
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

c....        Reduz as equacoes, exceto os termos da diagonal:
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

C     ****************************************************************
C     ****************************************************************
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

C     ****************************************************************
C     Funcao produto escalar de dois vetores
C     ****************************************************************
      real*8 function dot(a,b,n)
        integer n
        real*8 a(*), b(*)

        dot = 0.d0
        do 100 i = 1, n
            dot = dot + a(i) * b(i)
 100    continue
      return
      end

C     ****************************************************************
C     Subrotina para verificar a disponibilidade de memoria
C     ****************************************************************
      subroutine mem(npos)
        common /size/ max
        integer npos
        if ( (npos-1) .gt. max ) then
            print*, ' Memoria insuficiente !'
        stop
        endif
      return
      end

C     ****************************************************************
C     Subrotina para zerar um vetor de real de i2-i1 elementos
C     ****************************************************************
      subroutine mzero(m,i1,i2)
        integer i1,i2,m(*)
        do 100 i = 1, i2-i1 + 1
            m(i) = 0
 100    continue
      return
      end

C     ****************************************************************
C     Subrotina para zerar um vetor de real de i2-i1 elementos
C     ****************************************************************
      subroutine  azero(a,i1,i2)
        integer i1,i2
        real*8 a(*)
        do 100 i = 1,  i2 - i1 + 1
            a(i) = 0.d0
 100    continue
      return
      end

C     ****************************************************************
C     Subrotina para escrever os resultados
C     ****************************************************************
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
