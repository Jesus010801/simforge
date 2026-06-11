      PROGRAM AperR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxnat=1000000)
      parameter(maxn=100000)
      real*8 coorx,coory,coorz,coorxn,cooryn,coorzn
      real*8 boxx,boxy,boxz,Areanm2
      character*72 lin1,lin2,lin3
      character*79 line2
      character*1 cad,cad2,cadn
      integer numaa,numab,numac,num3,i,j,k,Aperes
      character*4 fin
      character*50 name,file


!______________________________________________________________________
!                  ARCHIVOS USADOS 

      open(11,file='area_2.dat')
      rewind 11

      open(12,file='areaAng2.dat')
      rewind 12

!______________________________________________________________________
!                  Lectura de Polcalcina 

      read(11,17)Areanm2
      Aperes=Areanm2*100
      write(12,18)Aperes



17    format(f5.3)
18    format(i4)

      return
      end
