
C      ***ABHIJIT MUKHERJEE ***
C      *** BUBBLE ***
C     ****************************
	SUBROUTINE SETUP3D
C	****************************
	IMPLICIT NONE
	INTEGER I,J,K
	INCLUDE "VAR3D.DEC"
	
c	***L&M&N ARE THE TOTAL NO. OF CV'S***
C	***L1&M1&N1 ARE TOTAL NO. OF GRID POINTS***
C	***L2&M2&N2 ARE TOTAL NO. OF INTERIOR GRID POINTS***
	L1=L+1
	M1=M+1
	L2=L1+1
	M2=M1+1
	N1=N+1
	N2=N1+1
!
	DX=LEN/L
	DY=HEI/M
	DZ=WID/N
C	PRINT *, "UO, LO, TO, REL, PRL, BOL, DX, DY, DZ"
C	PRINT *, UO, LO, TO, REL, PRL, BOL, DX, DY, DZ
C	***X DIRECTION GRID POINTS***
!	
C
	DO 110 I=1, L2
	XX(I)=(I-2)*DX
110	CONTINUE
C	XX(L2)=XX(L1)+DX/2.
C	***Y DIRECTION GRID POINTS***
!	
C	YY(1)=-DY/2.
C	YY(2)=0.
	DO 120 J=1, M2
	YY(J)=(J-2)*DY
120	CONTINUE
C	YY(M2)=YY(M1)+DY/2.
!
C	***Z DIRECTION GRID POINTS***
!	
C	ZZ(1)=-DZ/2.
C	ZZ(2)=0.
	DO 130 K=1, N2
	ZZ(K)=(K-2)*DZ
130	CONTINUE
C	ZZ(N2)=ZZ(N1)+DY/2.	
!

C	***X DIRECTION STAGGERED GRID***
!
	DO 210 I=1, L1
	XU(I)=-DX/2.+(I-1)*DX
210	CONTINUE
!
C	***Y DIRECTION STAGGERED GRID***
!
	DO 220 J=1,	M1
	YV(J)=-DY/2.+(J-1)*DY
220	CONTINUE
!
C	***Z DIRECTION STAGGERED GRID***
!
	DO 230 K=1,	N1
	ZW(K)=-DZ/2.+(K-1)*DZ
230	CONTINUE
!
	RETURN
	END
C	********************************************




	
