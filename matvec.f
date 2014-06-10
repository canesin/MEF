C     *****************************************************************
C     MATVEC produto matriz vertor c/ armazenamento skyline
C     *****************************************************************
      subroutine matvec(jdiag,a,vector,neq)
          implicit real*8 (a-h,o-z)
          integer i,j,clnefet,neq
	  real*8 vector(neq), aux(neq)
	  dimension a(*), jdiag(*)
	  do j=1, neq
		aux(j)=0
	  enddo
	  if(neq .ge. 1) then
		aux(1)= a(1)*vector(1)
	  endif
	  if(neq .ge. 2) then
		do j=2,neq
			clnefet= jdiag(j) - jdiag(j-1)
			aux(j)= aux(j) + a(jdiag(j)) * vector(j)
			if (clnefet .ge. 2) then
				do i=1, (clnefet-1)
				aux(j-i)= aux(j-i) + a(jdiag(j)-i) * vector(j)
				aux(j)= aux(j) + a(jdiag(j)-i) * vector(j-i)
			    enddo
		          endif
		enddo
	  endif
	  do j=1, neq
		vector(j)=aux(j)
	  enddo
      end

